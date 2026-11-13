#!/usr/bin/env python3
"""
miniML batch detection driven by the ONE master workbook (Ephys_master.xlsx)
==========================================================================

Adapted from the supervisor's notebook "batch analysis_SdK.ipynb" (same miniML settings, same pretrained
model call). What this version does:

  * Recordings come from the master sheet "files" (use_for_miniML == "Y") and are located on disk by
    path_on_mac (or by file name anywhere under --data-root), so the folder layout does not matter.
  * Each recording is run once per direction in files.miniML_directions ("negative", "positive" or
    "negative,positive"). Positive runs invert the trace inside miniML (pretrained model = negative events).
    Where the column is empty, a rule by holding is used (see directions_from_holding).
  * Results are ADDED to the master workbook, no separate CSV files:
        miniML_per_file          one row per recording x direction x settings
        miniML_events            one row per detected event
        miniML_per_cell_holding  rebuilt from the two sheets above
        runs                     one row per settings hash (settings, model, dates, counts)
    Recordings already in miniML_per_file with the same settings are skipped (use --redo to run again).
    Changing a setting or the model (e.g. a transfer-learned model) changes the settings hash, so the
    new results are added next to the old ones instead of replacing them.
  * While running, each finished recording is cached in <master folder>/_work/miniML_per_file/ so a
    stopped run can be resumed; the master workbook is written at the end (and every --save-every runs).
  * Every output row carries new_name, cell key, holding_command_mV, holding_flag and direction.
  * The master is backed up to _backups/ before each write. Close it in Excel first.

Usage (terminal, in an environment where miniML/TensorFlow work):

  python scripts/miniml_batch_master.py \
      --miniml-core "/path/miniML-main/core" \
      --model "/path/miniML-main/models/GC_lstm_model.h5" \
      --dry-run            # first: only checks the master against the files on disk
  (--master and --data-root default to the Ephys folder that holds scripts/)
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ephys_master as em

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
# Direction rule by COMMAND holding (set with Advaith, 12 Nov 2026):
#   KCl internal                      -> negative only (Cl- currents are inward; downward events = IPSCs)
#   holding >= 0 mV (0, 5, 10, 15)    -> positive only (outward IPSCs)
#   -50 ... -35 mV                    -> both directions (EPSCs inward + IPSCs outward expected)
#   <= -60 mV                         -> negative only (EPSCs)
#   -55 mV (gap)                      -> negative only (decided with Advaith, 12 Nov: few/no IPSCs expected)
#   -30 mV (gap)                      -> both directions (decided with Advaith, 12 Nov)
#   any other value in a gap          -> both directions (nothing lost); see GAP_OVERRIDES
BOTH_RANGE = (-50.0, -35.0)
NEG_ONLY_AT_OR_BELOW = -60.0
POS_ONLY_AT_OR_ABOVE = 0.0
GAP_RULE = "both"                   # fallback for any gap value not listed below
GAP_OVERRIDES = {-55.0: "negative", -30.0: "negative,positive"}


def directions_from_holding(holding, internal: str = "") -> str:
    if str(internal) == "KCl":
        return "negative"
    try:
        h = float(holding)
    except (TypeError, ValueError):
        return "negative,positive"          # unknown holding: be safe, run both
    if np.isnan(h):
        return "negative,positive"
    if h >= POS_ONLY_AT_OR_ABOVE:
        return "positive"
    if BOTH_RANGE[0] <= h <= BOTH_RANGE[1]:
        return "negative,positive"
    if h <= NEG_ONLY_AT_OR_BELOW:
        return "negative"
    if h in GAP_OVERRIDES:
        return GAP_OVERRIDES[h]
    if GAP_RULE == "both":
        return "negative,positive"
    return "negative" if h < np.mean([BOTH_RANGE[0], NEG_ONLY_AT_OR_BELOW]) + 1e-9 else "positive"


def default_directions(row) -> str:
    """Used only when the master table has no 'miniML_directions' value for a row."""
    return directions_from_holding(row.get("holding_command_mV"), row.get("internal", ""))


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


def process(row, path: Path, direction: str, model, P: dict, outdir: Path, h: str, use_cache: bool = True) -> str:
    from miniML import MiniTrace
    stem = f"{Path(row['new_name']).stem}__{direction}__{h}"
    done = outdir / f"{stem}.summary.json"
    if use_cache and done.exists():
        return "skipped (cached from an earlier, unfinished run)"
    trace = MiniTrace.from_axon_file(filename=str(path), channel=P["channel"],
                                     scaling=P["scaling"], unit=P["unit"])
    det = detect_one(trace, direction, model, P)
    ev = events_table(det)
    ev.to_csv(outdir / f"{stem}.events.csv", index=False)
    s = summarise(ev, trace, direction)
    sig = np.asarray(trace.data, dtype=float)
    s["noise_hf_pA"] = float(1.4826 * np.median(np.abs(np.diff(sig) - np.median(np.diff(sig)))) / np.sqrt(2))   # robust sample-to-sample noise SD
    s["noise_mad_pA"] = float(1.4826 * np.median(np.abs(sig - np.median(sig))))                                  # robust SD incl. slow fluctuation
    s["hold_current_median_pA"] = float(np.median(sig))
    s.update(new_name=row["new_name"], params_hash=h, status="ok" if len(ev) else "ok_no_events",
             processed_at=time.strftime("%Y-%m-%d %H:%M:%S"))
    done.write_text(json.dumps(s, default=float))
    return s["status"]


# --------------------------------------------------------------------------
# add the cached per-file results to the master workbook
# --------------------------------------------------------------------------
def to_master(master_path: Path, files: pd.DataFrame, cache: Path, h: str, settings: dict):
    sums = [json.loads(p.read_text()) for p in sorted(cache.glob(f"*__{h}.summary.json"))]
    if not sums:
        print("Nothing new to add to the master.")
        return
    fm = em.meta_table(files)
    new_pf = em.with_meta_per_file(pd.DataFrame(sums), fm)
    parts = []
    for p in sorted(cache.glob(f"*__{h}.events.csv")):
        nm, direction = p.name.split("__")[0:2]
        e = pd.read_csv(p)
        if len(e):
            e.insert(0, "direction", direction)
            e.insert(0, "new_name", nm + ".abf")
            parts.append(e)
    new_ev = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["new_name", "direction"] + EVENT_COLS)
    new_ev.insert(2, "params_hash", h)
    new_ev = em.with_meta_events(new_ev, fm)
    new_ev.insert(3, "params_hash", new_ev.pop("params_hash"))

    old_pf = em.read_sheet(master_path, "miniML_per_file")
    old_ev = em.read_sheet(master_path, "miniML_events")
    pf = em.merge_frames(old_pf, new_pf, em.KEYS["miniML_per_file"], sort=["date_header", "new_name", "direction"])
    ev = em.merge_frames(old_ev, new_ev, em.KEYS["miniML_events"], sort=["new_name", "direction", "peak_time_s"])
    cell = em.cell_holding_table(pf, ev)
    runs = em.read_sheet(master_path, "runs")
    row = pd.DataFrame([dict(tool="miniML", params_hash=h, model=settings.get("model"), settings=json.dumps(settings),
                             first_processed=pf.loc[pf.params_hash == h, "processed_at"].min(),
                             last_processed=pf.loc[pf.params_hash == h, "processed_at"].max(),
                             n_recording_directions=int((pf.params_hash == h).sum()),
                             n_events=int((ev.params_hash == h).sum()),
                             notes=settings.get("notes") or (
                                 runs.loc[(runs.tool == "miniML") & (runs.params_hash == h), "notes"].iloc[0]
                                 if len(runs) and ((runs.tool == "miniML") & (runs.params_hash == h)).any() else ""))])
    runs = em.merge_frames(runs, row, ["tool", "params_hash"])
    em.write_sheets(master_path, {"miniML_per_file": pf, "miniML_events": ev,
                                  "miniML_per_cell_holding": cell, "runs": runs})
    print(f"Master updated: miniML_per_file {len(pf)} rows (+{len(new_pf)} new/updated), "
          f"miniML_events {len(ev)}, miniML_per_cell_holding {len(cell)}")


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", type=Path, default=em.default_master())
    ap.add_argument("--data-root", type=Path, default=None, help="folder with the .abf files (default: the master's folder)")
    ap.add_argument("--work", type=Path, default=None, help="cache folder for unfinished runs (default <master folder>/_work/miniML_per_file)")
    ap.add_argument("--miniml-core", type=Path, help="folder containing miniML.py (miniML-main/core)")
    ap.add_argument("--model", type=Path, help="path to GC_lstm_model.h5 (or a transfer-learned model)")
    ap.add_argument("--only", help="process only recordings whose name contains this text (e.g. a date or cell key)")
    ap.add_argument("--names", help="comma-separated new_name values to process (testing a few recordings)")
    ap.add_argument("--limit", type=int, help="process at most N recordings (testing)")
    ap.add_argument("--max-hours", type=float, help="stop cleanly after this many hours (re-run the same command to resume)")
    ap.add_argument("--redo", action="store_true", help="run again even if the master already has this recording + direction + settings")
    ap.add_argument("--save-every", type=int, default=0, help="also write to the master every N finished runs (0 = only at the end)")
    ap.add_argument("--dry-run", action="store_true", help="check master vs files, do not run miniML")
    ap.add_argument("--aggregate-only", action="store_true", help="add the cached results to the master without running anything")
    ap.add_argument("--notes", default="", help="free text stored in the runs sheet for this settings hash")
    a = ap.parse_args()

    mp = a.master
    if not mp.exists():
        sys.exit(f"Master workbook not found: {mp}")
    root = a.data_root or mp.parent
    files_all = em.read_sheet(mp, "files")
    need = {"new_name", "use_for_miniML", "holding_command_mV", "holding_flag"}
    if need - set(files_all.columns):
        sys.exit(f"files sheet lacks columns: {sorted(need - set(files_all.columns))}")
    master = files_all[files_all["use_for_miniML"].astype(str).str.upper() == "Y"].copy()
    if "miniML_directions" not in master.columns:
        master["miniML_directions"] = ""
    master["miniML_directions"] = master.apply(
        lambda r: (str(r["miniML_directions"]).strip() if str(r["miniML_directions"]).strip() not in ("", "nan")
                   else default_directions(r)), axis=1)
    stem = master["new_name"].astype(str).str.replace(r"\.abf$", "", regex=True)
    master["recording_index"] = stem.str.rsplit("_", n=1).str[1]
    if a.only:
        master = master[master["new_name"].str.contains(a.only, regex=False)]
    if a.names:
        wanted = [n.strip() for n in a.names.split(",") if n.strip()]
        master = master[master["new_name"].isin(wanted)]
    master = (master.assign(_ci=pd.to_numeric(master["recording_index"], errors="coerce"))
                    .sort_values(["date_header", "_ci"], ascending=False)
                    .drop(columns="_ci").reset_index(drop=True))
    if a.limit:
        master = master.head(a.limit)

    idx, dups = em.index_abf(root)
    paths = {r["new_name"]: em.locate(r, root, idx) for _, r in master.iterrows()}
    found = master["new_name"].map(lambda n: paths[n] is not None)
    print(f"Master rows to analyse: {len(master)} | found on disk: {int(found.sum())} | missing: {int((~found).sum())}")
    for nm in master.loc[~found, "new_name"].tolist()[:20]:
        print("  MISSING:", nm)
    nohold = master["holding_command_mV"].isna()
    if nohold.any():
        print(f"WARNING: {int(nohold.sum())} recordings have no holding_command_mV. They will run, but cannot be analysed by holding.")
    print(master.loc[found, "miniML_directions"].value_counts().to_string())
    if a.dry_run:
        return
    if not a.aggregate_only and (not a.model or not a.miniml_core):
        sys.exit("--model and --miniml-core are required (except with --dry-run / --aggregate-only)")

    work = a.work or (mp.parent / "_work" / "miniML_per_file")
    work.mkdir(parents=True, exist_ok=True)
    if a.model:
        h = params_hash(PARAMS, a.model)
        settings = dict(PARAMS, model=a.model.name, model_bytes=a.model.stat().st_size, notes=a.notes)
        (work / f"settings_{h}.json").write_text(json.dumps(settings, indent=2))
    else:   # --aggregate-only without a model: use the newest cached settings
        js = sorted(work.glob("settings_*.json"), key=lambda p: p.stat().st_mtime)
        if not js:
            sys.exit("--aggregate-only needs --model (or an earlier cached run)")
        h = js[-1].stem.split("_")[1]
        settings = json.loads(js[-1].read_text())
    print("Settings hash:", h, "| params:", PARAMS)

    if not a.aggregate_only:
        done_pf = em.read_sheet(mp, "miniML_per_file")
        have = set()
        if len(done_pf) and not a.redo:
            have = set(zip(done_pf.loc[done_pf.params_hash == h, "new_name"], done_pf.loc[done_pf.params_hash == h, "direction"]))
        todo = [(r, paths[r["new_name"]], d.strip()) for _, r in master[found].iterrows()
                for d in str(r["miniML_directions"]).split(",") if (r["new_name"], d.strip()) not in have]
        print(f"Detection runs to do: {len(todo)} (already in the master with these settings: "
              f"{sum(len(str(d).split(',')) for d in master.loc[found, 'miniML_directions']) - len(todo)})")
        if todo:
            sys.path.insert(0, str(a.miniml_core))
            import tensorflow as tf
            model = tf.keras.models.load_model(str(a.model))
            t_start = time.time()
            for i, (r, path, d) in enumerate(todo, 1):
                if (work / "STOP").exists():
                    print("STOP file found -> stopping cleanly after the last finished run. Delete the STOP file and re-run to resume.", flush=True)
                    break
                if a.max_hours and (time.time() - t_start) > a.max_hours * 3600:
                    print(f"--max-hours {a.max_hours} reached -> stopping cleanly; re-run the same command to resume.", flush=True)
                    break
                t0 = time.time()
                try:
                    st = process(r, path, d, model, PARAMS, work, h, use_cache=not a.redo)
                except Exception as e:  # keep going; record the failure
                    st = f"ERROR {type(e).__name__}: {e}"
                    (work / f"{Path(path).stem}__{d}__{h}.error.txt").write_text(traceback.format_exc())
                print(f"[{i}/{len(todo)}] {Path(path).name} [{d}] -> {st} ({time.time() - t0:.0f}s)", flush=True)
                if a.save_every and i % a.save_every == 0:
                    to_master(mp, files_all, work, h, settings)

    # only recordings of this selection that finished with these settings
    to_master(mp, files_all, work, h, settings)
    print("Done.")


if __name__ == "__main__":
    main()
