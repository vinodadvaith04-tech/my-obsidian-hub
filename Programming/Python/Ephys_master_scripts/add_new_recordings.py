"""
add_new_recordings.py - bring new .abf recordings into the master workbook (sheet files).

What it does
  * scans the data root (default: the folder that holds Ephys_master.xlsx) for .abf files
  * a file counts as KNOWN if its name matches new_name or old_name in the files sheet, or its relative path
    matches path_on_mac; everything else is NEW
  * every NEW file is appended as a row with what the ABF header says (date, time, protocol, mode, duration,
    sweeps, header holding) and status "NEW - fill in". Nothing is guessed: sex, slice, cell, internal,
    holding_command_mV stay empty and use_for_miniML = N until you fill them in (in Excel).
  * existing rows are never changed. A backup of the workbook goes to _backups/ first.

After you filled in a row (holding_command_mV, internal, use_for_miniML = Y, miniML_directions) run
  python scripts/miniml_batch_master.py        (it only runs what is not yet in the master)

Usage:  python scripts/add_new_recordings.py [--data-root FOLDER] [--exclude TEXT ...] [--dry-run]
"""
import argparse, sys
from pathlib import Path
import pandas as pd, pyabf

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ephys_master as em

SKIP_DIRS = ("_backups", "_work", "_archive", "scripts")

ap = argparse.ArgumentParser()
ap.add_argument("--master", type=Path, default=em.default_master())
ap.add_argument("--data-root", type=Path, default=None)
ap.add_argument("--exclude", action="append", default=[], help="skip paths containing this text (repeatable)")
ap.add_argument("--dry-run", action="store_true")
a = ap.parse_args()
root = a.data_root or a.master.parent

files = em.read_sheet(a.master, "files")
known_names = set(files.new_name.dropna()) | set(files.old_name.dropna())
known_paths = set(files.path_on_mac.dropna().astype(str)) if "path_on_mac" in files else set()

new = []
for p in sorted(root.rglob("*.abf")):
    rel = p.relative_to(root).as_posix()
    if p.name.startswith("._") or not p.is_file() or any(part.startswith(SKIP_DIRS) for part in p.relative_to(root).parts[:-1]):
        continue
    if any(x in rel for x in a.exclude):
        continue
    if p.name in known_names or rel in known_paths:
        continue
    new.append(p)

print(f"{len(files)} recordings in the master, {len(new)} new .abf files found under {root}")
rows = []
for p in new:
    rel = p.relative_to(root).as_posix()
    try:
        ab = pyabf.ABF(str(p), loadData=False)
        proto = ab.protocol or ""
        gap = "GapFree" in proto
        ptype = "sPSC gap-free" if gap else ("cap check" if "Cap" in proto else ("spike frequency" if "spike" in proto else "other"))
        try:
            hh = float(ab.holdingCommand[0])
        except Exception:
            hh = None
        rows.append({"new_name": p.name, "old_name": p.name, "old_folder": p.parent.relative_to(root).as_posix(),
                     "new_folder": p.parent.relative_to(root).as_posix(), "path_on_mac": rel,
                     "status": "NEW - fill in", "date_header": int(ab.abfDateTime.strftime("%Y%m%d")),
                     "start_time": ab.abfDateTime.strftime("%H:%M:%S"), "protocol": proto, "protocol_type": ptype,
                     "mode": "VC" if ab.adcUnits[0] == "pA" else "CC", "duration_s": round(ab.dataLengthSec, 1),
                     "sweeps": ab.sweepCount, "holding_header_mV": hh, "use_for_miniML": "N",
                     "analysis_note": "new recording: sex, slice, cell, internal, holding_command_mV not filled in yet"})
    except Exception as e:
        rows.append({"new_name": p.name, "old_name": p.name, "path_on_mac": rel, "status": "NEW - header unreadable",
                     "use_for_miniML": "N", "analysis_note": f"pyabf failed: {e}"})
    print(" NEW", rel)

if not rows:
    print("nothing to add")
    sys.exit(0)
if a.dry_run:
    print("dry run, master not changed")
    sys.exit(0)
out = em.merge_frames(files, pd.DataFrame(rows), em.KEYS["files"])
em.write_sheets(a.master, {"files": out})
print(f"added {len(rows)} rows to sheet files (now {len(out)}). Fill in the empty columns in Excel, then run miniml_batch_master.py")
