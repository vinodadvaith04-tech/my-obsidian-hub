#!/usr/bin/env python3
"""
miniML batch detection driven by the ephys master table
=======================================================

Adapted from the supervisor's notebook "batch analysis_SdK.ipynb" (same miniML
settings, same pretrained model call), with these changes:

  * Files are chosen from the master table (sheet "files", use_for_miniML == "Y")
    and located by NAME anywhere under --data-root, so the folder layout does not
    matter (the supervisor's glob only looked in one flat folder).
  * Each recording is run once per direction listed in the master column
    "miniML_directions" ("negative", "positive" or "negative,positive").
    Both-direction runs make it possible to quantify inward AND outward events
    in the same trace. The pretrained model was trained on negative events;
    positive runs invert the trace inside miniML (training_direction='negative').
  * Every output row (per file, per event, per cell/holding) carries new_name,
    cell key, holding_command_mV, holding_flag and direction, so any detected
    event can be traced back to its holding.
  * Empty recordings give NaN amplitudes (not 0) and frequency 0.
  * Results are keyed by a hash of the detection settings + model file, so
    changing a setting or the model re-runs the files instead of silently
    skipping stale pickles. Unchanged settings are skipped (resumable).
  * No machine-specific paths: everything comes from the command line.
  * Never writes into the data folders and never edits the master table
    (--update-master writes a COPY).

Usage (from a terminal, in the environment where miniML/TensorFlow work):

  python miniml_batch_master.py \
      --master  "/path/Ephys_master_table_DRAFT.xlsx" \
      --data-root "/path/Ephys pilot & practice" \
      --results "/path/miniML_results" \
      --miniml-core "/path/miniML-main/core" \
      --model "/path/miniML-main/models/GC_lstm_model.h5" \
      --dry-run            # first: only checks the master table against the files

Outputs in --results:
  per_file/<stem>__<direction>__<hash>.{summary.json,events.csv,pickle}
  miniML_file_summary.csv     one row per recording x direction
  miniML_events_all.csv       one row per detected event
  miniML_cell_holding_summary.csv   per cell x holding x direction
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

# miniML's core still calls np.trapz, which NumPy >= 2 renamed to np.trapezoid.
if not hasattr(np, "trapz"):
    np.trapz = np.trapezoid

# --------------------------------------------------------------------------
# Detection settings = the supervisor's notebook settings (unchanged)
# --------------------------------------------------------------------------
PARAMS = dict(
    window_size=600,
    model_threshold=0.5,
    batch_size=512,
    peak_w=5,
    rel_prom_cutoff=0.25,
    convolve_win=20,
    gradient_convolve_win=40,
    channel=0,
    scaling=1.0,
    unit="pA",
)
# The supervisor's QC'd recordings had bad segments overwritten with 10000.
# Our recordings have no such cuts, so this is off. If set, those samples are
# excluded from the valid time (NOTE: they would still sit inside miniML's
# min-max scaled windows, so cut segments should rather be handled by splitting).
SENTINEL = None

META_COLS = ["new_name", "session", "date_header", "T.Code", "sex", "slice", "cell",
             "internal", "psc_content", "drug", "holding_command_mV", "holding_flag",
             "ljp_in_note", "duration_s"]


# --------------------------------------------------------------------------
# master table / file lookup
# --------------------------------------------------------------------------
def default_directions(row) -> str:
    """Event directions to run when the master table has no 'miniML_directions'.

    KCl internal: Cl- currents are INWARD at -60 mV command, so IPSCs are negative too.
    KGlu / Cs-gluconate at intermediate or depolarised holdings: IPSCs are outward.
    """
    c = str(row.get("psc_content", "") or "")
    if c.startswith("EPSC"):
        return "negative"
    if c.startswith("IPSC") or "KCl" in c:
        return "negative"
    return "negative,positive"


def load_master(path: Path) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="files",
                       dtype={"date_header": str, "counter": str, "cell": str, "slice": str})
    need = {"new_name", "use_for_miniML", "holding_command_mV", "holding_flag"}
    missing = need - set(df.columns)
    if missing:
        sys.exit(f"Master table lacks columns: {sorted(missing)}")
    df = df[df["use_for_miniML"].astype(str).str.upper() == "Y"].copy()
    if "miniML_directions" not in df.columns:
        df["miniML_directions"] = ""
    df["miniML_directions"] = df.apply(
        lambda r: (str(r["miniML_directions"]).strip()
                   if str(r["miniML_directions"]).strip() not in ("", "nan")
                   else default_directions(r)), axis=1)
    stem = df["new_name"].str.replace(r"\.abf$", "", regex=True)
    df["recording_name"] = stem.str.rsplit("_", n=1).str[0]   # = cell key (supervisor's rsplit)
    df["recording_index"] = stem.str.rsplit("_", n=1).str[1]  # = per-day counter
    return df.reset_index(drop=True)


def index_abf(root: Path):
    idx, dups = {}, collections.defaultdict(list)
    for p in root.rglob("*.abf"):
        if p.name in idx:
            dups[p.name].append(p)
        else:
            idx[p.name] = p
    return idx, dups


def params_hash(params: dict, model_path: Path) -> str:
    sig = dict(params, model=model_path.name, model_bytes=model_path.stat().st_size)
    return hashlib.sha1(json.dumps(sig, sort_keys=True).encode()).hexdigest()[:8]


# --------------------------------------------------------------------------
# detection (needs TensorFlow + miniML; imported lazily so --dry-run works anywhere)
# --------------------------------------------------------------------------
def detect_one(trace, direction: str, model, P: dict):
    """Run miniML on one MiniTrace in one direction. Returns the EventDetection object."""
    from miniML import EventDetection
    det = EventDetection(
        data=trace, model=model, model_threshold=P["model_threshold"],
        window_size=P["window_size"], batch_size=P["batch_size"],
        event_direction=direction, verbose=0)
    det.detect_events(eval=True, peak_w=P["peak_w"], rel_prom_cutoff=P["rel_prom_cutoff"],
                      convolve_win=P["convolve_win"], gradient_convolve_win=P["gradient_convolve_win"])
    return det


EVENT_COLS = ["event_idx", "peak_time_s", "start_time_s", "interevent_interval_s", "amplitude",
              "amplitude_abs", "score", "charge", "risetime_ms", "slope", "halfdecay_ms", "halfwidth_ms"]


def events_table(det) -> pd.DataFrame:
    """One row per detected event (empty frame if none)."""
    ev = getattr(det, "events", None)
    if ev is None or len(ev) == 0 or not hasattr(det, "event_stats"):
        return pd.DataFrame(columns=EVENT_COLS)
    es = det.event_stats
    amp = np.asarray(es.amplitudes, dtype=float)
    n = len(amp)
    cols = dict(
        peak_time_s=det.event_peak_times, start_time_s=det.event_start_times,
        interevent_interval_s=det.interevent_intervals, score=es.event_scores,
        charge=es.charges, risetime_ms=np.asarray(es.risetimes, float) * 1000,
        slope=es.slopes, halfdecay_ms=np.asarray(es.halfdecays, float) * 1000,
        halfwidth_ms=np.asarray(es.halfwidths, float) * 1000)
    out = {"event_idx": np.arange(n), "amplitude": amp, "amplitude_abs": np.abs(amp)}
    for k, v in cols.items():
        v = np.asarray(v, dtype=float).ravel()
        if v.size != n:
            raise ValueError(f"event field {k} has {v.size} values for {n} events")
        out[k] = v
    return pd.DataFrame(out)[EVENT_COLS]


def summarise(ev: pd.DataFrame, trace, direction: str) -> dict:
    sig = np.asarray(trace.data)
    total_s = sig.size * trace.sampling
    invalid_s = float(np.sum(sig == SENTINEL) * trace.sampling) if SENTINEL is not None else 0.0
    valid_s = total_s - invalid_s
    n = len(ev)
    nan = np.nan
    return dict(
        direction=direction, unit=trace.y_unit, sampling_interval_s=trace.sampling,
        recording_time_s=total_s, cut_time_s=invalid_s, valid_time_s=valid_s,
        n_events=n, frequency_Hz=(n / valid_s if valid_s > 0 else nan),
        amplitude_abs_mean=(ev["amplitude_abs"].mean() if n else nan),
        amplitude_abs_median=(ev["amplitude_abs"].median() if n else nan),
        charge_mean=(ev["charge"].abs().mean() if n else nan),
        risetime_ms_mean=(ev["risetime_ms"].mean() if n else nan),
        halfdecay_ms_mean=(ev["halfdecay_ms"].mean() if n else nan),
        halfwidth_ms_mean=(ev["halfwidth_ms"].mean() if n else nan))


def process(row, path: Path, direction: str, model, P: dict, outdir: Path, h: str) -> str:
    from miniML import MiniTrace
    stem = f"{path.stem}__{direction}__{h}"
    done = outdir / f"{stem}.summary.json"
    if done.exists():
        return "skipped (already done with these settings)"
    trace = MiniTrace.from_axon_file(filename=str(path), channel=P["channel"],
                                     scaling=P["scaling"], unit=P["unit"])
    det = detect_one(trace, direction, model, P)
    ev = events_table(det)
    ev.to_csv(outdir / f"{stem}.events.csv", index=False)
    if len(ev):
        det.save_to_pickle(filename=str(outdir / f"{stem}.pickle"),
                           include_prediction=False, include_data=False)
    s = summarise(ev, trace, direction)
    s.update(new_name=row["new_name"], params_hash=h, status="ok" if len(ev) else "ok_no_events",
             processed_at=time.strftime("%Y-%m-%d %H:%M:%S"))
    done.write_text(json.dumps(s, default=float))
    return s["status"]


# --------------------------------------------------------------------------
# aggregation (works without TensorFlow, from the per-file outputs)
# --------------------------------------------------------------------------
def aggregate(master: pd.DataFrame, results: Path, h: str):
    pf = results / "per_file"
    sums = [json.loads(p.read_text()) for p in sorted(pf.glob(f"*__{h}.summary.json"))]
    if not sums:
        print("Nothing to aggregate for these settings.")
        return None, None, None
    meta = [c for c in META_COLS + ["recording_name", "recording_index", "miniML_directions"]
            if c in master.columns]
    s = pd.DataFrame(sums).merge(master[meta], on="new_name", how="left")
    lead = ["new_name", "recording_name", "recording_index", "direction"]
    s = s[lead + [c for c in meta if c not in lead] +
          [c for c in s.columns if c not in lead + meta]]
    s.to_csv(results / "miniML_file_summary.csv", index=False)

    parts = []
    for p in sorted(pf.glob(f"*__{h}.events.csv")):
        nm, direction = p.name.split("__")[0:2]
        e = pd.read_csv(p)
        if len(e):
            e.insert(0, "direction", direction)
            e.insert(0, "new_name", nm + ".abf")
            parts.append(e)
    ev = (pd.concat(parts, ignore_index=True) if parts
          else pd.DataFrame(columns=["new_name", "direction"] + EVENT_COLS))
    em = master[["new_name", "recording_name", "holding_command_mV", "holding_flag",
                 "internal", "psc_content", "sex", "T.Code", "date_header", "cell"]]
    ev = ev.merge(em, on="new_name", how="left")
    ev.to_csv(results / "miniML_events_all.csv", index=False)

    g = s.groupby(["recording_name", "holding_command_mV", "direction"], dropna=False)
    cell = g.agg(n_files=("new_name", "nunique"), n_events=("n_events", "sum"),
                 valid_time_s=("valid_time_s", "sum"),
                 session=("session", "first"), date_header=("date_header", "first"),
                 sex=("sex", "first"), T_code=("T.Code", "first"), internal=("internal", "first"),
                 psc_content=("psc_content", "first"), holding_flag=("holding_flag", "first")).reset_index()
    cell["frequency_Hz"] = cell["n_events"] / cell["valid_time_s"]
    if len(ev):
        a = (ev.groupby(["recording_name", "holding_command_mV", "direction"], dropna=False)
               ["amplitude_abs"].agg(amplitude_abs_mean="mean", amplitude_abs_median="median").reset_index())
        cell = cell.merge(a, on=["recording_name", "holding_command_mV", "direction"], how="left")
    cell.to_csv(results / "miniML_cell_holding_summary.csv", index=False)
    return s, ev, cell


def write_master_copy(master_path: Path, s, cell, results: Path):
    out = master_path.with_name(master_path.stem + "_with_detection.xlsx")
    import shutil
    shutil.copy(master_path, out)
    with pd.ExcelWriter(out, engine="openpyxl", mode="a", if_sheet_exists="replace") as xw:
        s.to_excel(xw, sheet_name="detection_per_file", index=False)
        cell.to_excel(xw, sheet_name="detection_per_cell_holding", index=False)
    print("Master copy with detection sheets:", out)


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", required=True, type=Path)
    ap.add_argument("--data-root", required=True, type=Path)
    ap.add_argument("--results", required=True, type=Path)
    ap.add_argument("--miniml-core", type=Path, help="folder containing miniML.py (miniML-main/core)")
    ap.add_argument("--model", type=Path, help="path to GC_lstm_model.h5 (or a transfer-learned model)")
    ap.add_argument("--only", help="process only recordings whose name contains this text (e.g. a date or cell key)")
    ap.add_argument("--limit", type=int, help="process at most N recordings (testing)")
    ap.add_argument("--dry-run", action="store_true", help="check master vs files, do not run miniML")
    ap.add_argument("--aggregate-only", action="store_true", help="rebuild the summary CSVs from existing per-file results")
    ap.add_argument("--update-master", action="store_true", help="also write <master>_with_detection.xlsx (a copy)")
    a = ap.parse_args()

    master = load_master(a.master)
    if a.only:
        master = master[master["new_name"].str.contains(a.only, regex=False)].reset_index(drop=True)
    if a.limit:
        master = master.head(a.limit)
    idx, dups = index_abf(a.data_root)
    found = master["new_name"].isin(idx)
    print(f"Master rows to analyse: {len(master)} | found on disk: {int(found.sum())} | missing: {int((~found).sum())}")
    for nm in master.loc[~found, "new_name"].tolist()[:20]:
        print("  MISSING:", nm)
    for nm, ps in dups.items():
        if nm in set(master["new_name"]):
            sys.exit(f"Duplicate file name on disk (must be unique): {nm} -> {ps}")
    nohold = master["holding_command_mV"].isna()
    if nohold.any():
        print(f"WARNING: {int(nohold.sum())} recordings have no holding_command_mV (flag: "
              f"{master.loc[nohold, 'holding_flag'].unique().tolist()}). They will run, but cannot be analysed by holding.")
    runs = sum(len(str(d).split(",")) for d in master.loc[found, "miniML_directions"])
    print(f"Detection runs (recording x direction): {runs}")
    print(master.loc[found, "miniML_directions"].value_counts().to_string())
    unknown_in_master = {n for n in idx if n[:8].isdigit() and "_Adv00_" in n} - set(
        pd.read_excel(a.master, sheet_name="files")["new_name"].dropna())
    if unknown_in_master:
        print(f"NOTE: {len(unknown_in_master)} correctly-named files on disk are not in the master table.")
    if a.dry_run:
        return
    if not a.model or not a.miniml_core:
        sys.exit("--model and --miniml-core are required unless --dry-run")

    results = a.results
    (results / "per_file").mkdir(parents=True, exist_ok=True)
    h = params_hash(PARAMS, a.model)
    print("Settings hash:", h, "| params:", PARAMS)
    (results / f"settings_{h}.json").write_text(json.dumps(
        dict(PARAMS, model=a.model.name, model_bytes=a.model.stat().st_size), indent=2))

    if not a.aggregate_only:
        sys.path.insert(0, str(a.miniml_core))
        import tensorflow as tf
        model = tf.keras.models.load_model(str(a.model))
        todo = [(r, idx[r["new_name"]], d.strip()) for _, r in master[found].iterrows()
                for d in str(r["miniML_directions"]).split(",")]
        for i, (r, path, d) in enumerate(todo, 1):
            t0 = time.time()
            try:
                st = process(r, path, d, model, PARAMS, results / "per_file", h)
            except Exception as e:  # keep going; record the failure
                st = f"ERROR {type(e).__name__}: {e}"
                (results / "per_file" / f"{path.stem}__{d}__{h}.error.txt").write_text(traceback.format_exc())
            print(f"[{i}/{len(todo)}] {path.name} [{d}] -> {st} ({time.time() - t0:.0f}s)", flush=True)

    s, ev, cell = aggregate(master, results, h)
    if s is not None and a.update_master:
        write_master_copy(a.master, s, cell, results)
    print("Done. Results in", results)


if __name__ == "__main__":
    main()
