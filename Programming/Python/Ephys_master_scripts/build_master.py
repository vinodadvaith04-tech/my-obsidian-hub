"""
build_master.py - ONE-OFF migration (13 Nov 2026): put everything analysed so far into one workbook.

Reads the earlier separate files and writes Ephys_master.xlsx. Never changes or deletes its inputs.

Usage (from the Ephys folder):   python3 scripts/build_master.py [--folder "<Ephys folder>"]

Inputs (relative to the Ephys folder):
  Ephys_master_table_DRAFT_plus_slice4.xlsx            sheets files, cells, rename_log
  Cs_pilot_package/results_older_Cs_miniML/*.csv       miniML, all recordings before slice 4 (both directions for Cs)
  Cs_pilot_package/results_slice4_miniML/*.csv         miniML, slice 4 (both directions)
  oldDetector_results/*.csv, Cs_pilot_package/results_slice4_second_detector/*.csv   second detector
  Cs_pilot_package/Cs_report_tables_and_figures/cs_qc_manual.csv
"""
import argparse, json, sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ephys_master as em

ap = argparse.ArgumentParser()
ap.add_argument("--folder", type=Path, default=Path(__file__).resolve().parent.parent)
ap.add_argument("--out", type=Path, default=None)
a = ap.parse_args()
R = a.folder
OUT = a.out or (R / em.MASTER_NAME)
if OUT.exists():
    sys.exit(f"{OUT} already exists - not overwriting. Use the other scripts to add to it.")

SRC = R / "Ephys_master_table_DRAFT_plus_slice4.xlsx"
files = pd.read_excel(SRC, sheet_name="files")
cells = pd.read_excel(SRC, sheet_name="cells")
rename_log = pd.read_excel(SRC, sheet_name="rename_log")
print("files", files.shape, "cells", cells.shape, "rename_log", rename_log.shape)

# ---------------------------------------------------------------- files: where is each recording now?
idx, dups = em.index_abf(R)
rel = {}
for i, r in files.iterrows():
    p = None
    for k in ("new_name", "old_name"):
        n = r.get(k)
        if isinstance(n, str) and n in idx:
            p = idx[n]; break
    rel[i] = p.relative_to(R).as_posix() if p else np.nan
files.insert(files.columns.get_loc("new_folder") + 1, "path_on_mac", pd.Series(rel))
print("recordings found on disk:", files["path_on_mac"].notna().sum(), "of", len(files),
      "| duplicates by name:", len(dups))
files["analysis_note"] = np.nan   # free text for you

# ---------------------------------------------------------------- cells: add the slice-4 cells
have = set(zip(cells["date"].astype(int), cells["cell_new"].astype(int)))
CELL_FLAGS = {1: "short: only 0 and +10 mV, no negative holdings",
              2: "most complete and stable cell; no -70 mV; cap check missing in the files",
              3: "unstable early (large holding current, high noise at -80 and -50 mV), stable from -40 mV",
              4: "cell lost after two files at -70 mV",
              5: "leaky (about -1.8 nA at -80 mV); file-to-holding mapping ambiguous (files 56, 59, 60, 62)",
              6: "unstable early (high noise at -70 and -50 mV), stable from -30 mV"}
new_cells = []
s4 = files[(files.session == "Cs internal pilot") & (files.date_header == 20261112)]
for c, g in s4.groupby("cell"):
    c = int(c)
    if (20261112, c) in have:
        continue
    gf = g[g.use_for_miniML == "Y"].sort_values("counter")
    new_cells.append({"date": 20261112, "session": "Cs internal pilot", "slice": 4, "cell_new": c,
                      "cell_old_folder": float(c), "T.Code": "NT", "sex": "unk", "internal": "Cs-gluconate",
                      "internal_source": "stated", "n_files": len(g), "n_sPSC_gapfree": len(gf),
                      "holding_commands (gap-free, time order)": ", ".join(f"{h:g}" for h in gf.holding_command_mV),
                      "flags": "12 Nov note: " + CELL_FLAGS.get(c, "")})
cells = pd.concat([cells, pd.DataFrame(new_cells)], ignore_index=True)
print("cells now", len(cells), "(added", len(new_cells), ")")

# ---------------------------------------------------------------- miniML
def load(d, f):
    return pd.read_csv(R / d / f, low_memory=False)

pf_old = load("Cs_pilot_package/results_older_Cs_miniML", "miniML_file_summary.csv")
pf_s4 = load("Cs_pilot_package/results_slice4_miniML", "miniML_file_summary.csv")
ev_old = load("Cs_pilot_package/results_older_Cs_miniML", "miniML_events_all.csv")
ev_s4 = load("Cs_pilot_package/results_slice4_miniML", "miniML_events_all.csv")
per_file = pd.concat([pf_old, pf_s4], ignore_index=True).drop_duplicates(["new_name", "direction", "params_hash"], keep="last")
events = pd.concat([ev_old, ev_s4], ignore_index=True).drop_duplicates()
# The Cs-only result files carry meta columns only for the Cs recordings (they were merged with a Cs subset
# of the master). Rebuild the meta columns of ALL rows from the `files` sheet, so every row is complete.
META = ["session", "date_header", "T.Code", "sex", "slice", "cell", "internal", "psc_content", "drug",
        "holding_command_mV", "holding_flag", "ljp_in_note", "duration_s", "miniML_directions"]
fm = files.dropna(subset=["new_name"]).drop_duplicates("new_name").set_index("new_name")[META].copy()
stem = fm.index.to_series().str.replace(r"\.abf$", "", regex=True)
fm.insert(0, "recording_index", stem.str.rsplit("_", n=1).str[1].values)
fm.insert(0, "recording_name", stem.str.rsplit("_", n=1).str[0].values)
print("rows with missing metadata before the rebuild: per_file", int(per_file.recording_name.isna().sum()),
      "| events", int(events.recording_name.isna().sum()))
drop = ["recording_name", "recording_index", "event_type_assumed"] + META
res_cols = [c for c in per_file.columns if c not in drop + ["new_name", "direction"]]
per_file = per_file[["new_name", "direction"] + res_cols].merge(fm, left_on="new_name", right_index=True, how="left")
per_file.insert(4, "event_type_assumed", [em.assumed_event_type(d, i) for d, i in zip(per_file.direction, per_file.internal)])
lead = ["new_name", "recording_name", "recording_index", "direction", "event_type_assumed"] + META
per_file = per_file[lead + [c for c in per_file.columns if c not in lead]]
ev_meta = ["recording_name", "holding_command_mV", "holding_flag", "internal", "psc_content", "sex", "T.Code", "date_header", "cell"]
ev_res = [c for c in events.columns if c not in ev_meta + ["event_type_assumed"]]
events = events[ev_res].merge(fm[ev_meta], left_on="new_name", right_index=True, how="left")
events.insert(2, "event_type_assumed", [em.assumed_event_type(d, i) for d, i in zip(events.direction, events.internal)])
events = events[[c for c in events.columns if c not in ev_meta] + ev_meta]
assert per_file.recording_name.notna().all() and events.recording_name.notna().all()
# every event belongs to one recording x direction x settings (a single settings hash so far)
hashes = per_file["params_hash"].unique().tolist()
assert len(hashes) == 1, hashes
events.insert(3, "params_hash", hashes[0])
per_file = per_file.sort_values(["date_header", "new_name", "direction"]).reset_index(drop=True)
events = events.sort_values(["new_name", "direction", "peak_time_s"]).reset_index(drop=True)
# consistency: events per recording x direction must equal n_events in the per-file table
cnt = events.groupby(["new_name", "direction"]).size().rename("n")
chk = per_file.set_index(["new_name", "direction"]).n_events.to_frame().join(cnt).fillna({"n": 0})
bad = chk[chk.n_events != chk.n]
assert len(bad) == 0, bad.head()
cell_hold = em.cell_holding_table(per_file, events)
print("miniML: per_file", per_file.shape, "events", events.shape, "cell_holding", cell_hold.shape)

# check against the earlier separate summaries
old_ch = load("miniML_results", "miniML_cell_holding_summary.csv")
m = old_ch.merge(cell_hold, on=["recording_name", "holding_command_mV", "direction"], suffixes=("_old", "_new"))
diff_ev = m[m.n_events_old != m.n_events_new]
print("cell-holding rows of the first run:", len(old_ch), "| matched:", len(m),
      "| rows with a different event count (extra directions added later):", len(diff_ev))

# ---------------------------------------------------------------- second detector
ad_f = pd.concat([load("oldDetector_results", "oldDetector_file_summary.csv"),
                  load("Cs_pilot_package/results_slice4_second_detector", "oldDetector_file_summary.csv")]
                 ).drop_duplicates("new_name", keep="last").reset_index(drop=True)
ad_e = pd.concat([load("oldDetector_results", "oldDetector_events_all.csv"),
                  load("Cs_pilot_package/results_slice4_second_detector", "oldDetector_events_all.csv")]
                 ).drop_duplicates().reset_index(drop=True)
print("second detector: per_file", ad_f.shape, "events", ad_e.shape)

# ---------------------------------------------------------------- manual QC (Cs pilot notes)
qc = load("Cs_pilot_package/Cs_report_tables_and_figures", "cs_qc_manual.csv")
qc.insert(0, "session", "Cs internal pilot")

# ---------------------------------------------------------------- runs
settings = json.loads((R / "miniML_results" / f"settings_{hashes[0]}.json").read_text())
runs = pd.DataFrame([
    dict(tool="miniML", params_hash=hashes[0], model=settings.get("model"), settings=json.dumps(settings),
         first_processed=per_file["processed_at"].min(), last_processed=per_file["processed_at"].max(),
         n_recording_directions=len(per_file), n_events=len(events),
         notes="pretrained model, supervisor's notebook settings; positive runs invert the trace. "
               "Direction per recording is in files.miniML_directions (all Cs recordings: both)"),
    dict(tool="alt_det (second detector)", params_hash="alt_v1", model="matched filter + fit (psc_detection.py)",
         settings="", first_processed=np.nan, last_processed=np.nan,
         n_recording_directions=len(ad_f), n_events=len(ad_e),
         notes="cross-check of miniML; both polarities always; not for sharing with Sam")])

# ---------------------------------------------------------------- README
README = pd.DataFrame([
    ("What this is", "", "The ONE master workbook: every recording, every analysis result so far. Scripts read from it and add to it; filter any sheet in Excel to choose what to analyse."),
    ("Sheet", "runs", "One row per detection setup (settings hash, model, dates, counts). A new model or new settings = a new row; old results stay."),
    ("Sheet", "files", "One row per recording (.abf). Identity, session, protocol, internal, drug, holding (command, uncorrected for the junction potential), and use_for_miniML / miniML_directions that steer the detection. path_on_mac = where the file is now (relative to the Ephys folder). Edit holdings here by hand if needed."),
    ("Sheet", "cells", "One row per cell: files, gap-free files, holdings in time order, flags."),
    ("Sheet", "rename_log", "Old path -> new path of the 842 renamed files (the 86 slice-4 files and 9 ignored files were not renamed on the Mac)."),
    ("Sheet", "manual_qc", "Manual Ra / Rm / Cm / pipette checks read from the lab notes (Cs pilot so far)."),
    ("Sheet", "miniML_per_file", "miniML result per recording x direction x settings: events, frequency, amplitude, kinetics, noise, holding current."),
    ("Sheet", "miniML_events", "One row per detected miniML event, with cell and holding. params_hash links it to its settings."),
    ("Sheet", "miniML_per_cell_holding", "miniML results pooled per cell x holding x direction x settings (rebuilt from the two sheets above)."),
    ("Sheet", "alt_det_per_file", "Second detector, per recording (inward and outward). Cross-check only; do not share with Sam."),
    ("Sheet", "alt_det_events", "Second detector, one row per event. Cross-check only; do not share with Sam."),
    ("Keys", "", "files: new_name (old_name for ignored files). miniML_per_file: new_name + direction + params_hash. miniML_events: new_name + direction + params_hash + event_idx. alt_det_per_file: new_name."),
    ("Directions", "", "negative = inward (downward) events, positive = outward (upward) events. event_type_assumed (EPSC/IPSC) is an assumption from the sign of the event, not a measurement."),
    ("Holdings", "", "All holdings are COMMAND holdings (no junction-potential correction). The 4 Nov note shows a -16 mV offset (-60 read as about -76 mV); not checked for other days."),
    ("Workflow: new data", "", "1) put the .abf files anywhere under the Ephys folder. 2) python scripts/add_new_recordings.py - adds rows to files (status NEW). 3) fill in session, internal, holding_command_mV, use_for_miniML and miniML_directions in files (Excel). 4) python scripts/miniml_batch_master.py ... - runs the detection and adds the results to the sheets above."),
    ("Backups", "", "Every script that writes here first copies the workbook to _backups/ (last 5 kept). Close the workbook in Excel before running a script."),
    ("Built", "", "13 Nov 2026 from Ephys_master_table_DRAFT_plus_slice4.xlsx and the separate miniML / second-detector / QC tables (kept in _archive_before_master/)."),
], columns=["topic", "name", "description"])

em.write_sheets(OUT, {"README": README, "runs": runs, "files": files, "cells": cells, "rename_log": rename_log,
                      "manual_qc": qc, "miniML_per_file": per_file, "miniML_events": events,
                      "miniML_per_cell_holding": cell_hold, "alt_det_per_file": ad_f, "alt_det_events": ad_e},
                backup=False)
print("written", OUT)
