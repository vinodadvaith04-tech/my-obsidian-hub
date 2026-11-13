"""
run_alt_detector.py - the second detector (matched filter + fit, psc_detection.py) on the recordings of the
master workbook. A cross-check of miniML, both polarities always. Not for sharing with Sam.

Reads  : Ephys_master.xlsx, sheet files (use_for_miniML == Y); psc_detection.py (same folder as this script)
Writes : sheets alt_det_per_file and alt_det_events of the same workbook (added, never replaced unless --redo)
Cache  : <master folder>/_work/alt_det_per_file/ (so a stopped run can be resumed)

Usage:  python scripts/run_alt_detector.py [--data-root <folder with .abf>] [--only TEXT] [--max-hours H] [--redo]
Stop cleanly by creating a file named STOP in the cache folder.
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np, pandas as pd, pyabf

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ephys_master as em

ap = argparse.ArgumentParser()
ap.add_argument("--master", type=Path, default=em.default_master())
ap.add_argument("--data-root", type=Path, default=None)
ap.add_argument("--psc-dir", type=Path, default=Path(__file__).resolve().parent, help="folder holding psc_detection.py")
ap.add_argument("--only", help="only recordings whose name contains this text")
ap.add_argument("--max-hours", type=float)
ap.add_argument("--redo", action="store_true")
ap.add_argument("--dry-run", action="store_true")
a = ap.parse_args()
mp = a.master
root = a.data_root or mp.parent
files = em.read_sheet(mp, "files")
m = files[files.use_for_miniML.astype(str).str.upper() == "Y"].copy()
if a.only:
    m = m[m.new_name.str.contains(a.only, regex=False)]
have = set()
if not a.redo:
    have = set(em.read_sheet(mp, "alt_det_per_file").get("new_name", pd.Series(dtype=str)).dropna())
m = m[~m.new_name.isin(have)]
m["_c"] = pd.to_numeric(m["counter"], errors="coerce")
m = m.sort_values(["date_header", "_c"], ascending=False)
idx, _ = em.index_abf(root)
print(f"recordings to run: {len(m)} (already in the master: {len(have)})")
if a.dry_run or not len(m):
    sys.exit(0)
sys.path.insert(0, str(a.psc_dir))
import psc_detection as P

cache = mp.parent / "_work" / "alt_det_per_file"
cache.mkdir(parents=True, exist_ok=True)
t0 = time.time()
for i, (_, r) in enumerate(m.iterrows(), 1):
    if (cache / "STOP").exists():
        print("STOP file found"); break
    if a.max_hours and time.time() - t0 > a.max_hours * 3600:
        print("max-hours reached"); break
    f = cache / (r.new_name[:-4] + ".old.json")
    if f.exists() and not a.redo:
        continue
    p = em.locate(r, root, idx)
    if p is None:
        print("MISSING", r.new_name); continue
    ts = time.time()
    ab = pyabf.ABF(str(p)); ab.setSweep(0); I = np.asarray(ab.sweepY, float); fs = float(ab.dataRate)
    ok, rej, di = P.detect(I, fs)
    ev = [dict(new_name=r.new_name, polarity=e["polarity"], peak_s=e["peak_s"], amp_pA=e["amp_pA"],
               rise1090_ms=e["rise1090_ms"], tau_decay_ms=e["td"], halfwidth_ms=e["halfwidth_ms"],
               snr=e["snr"], isolated=bool(e["isolated"])) for e in ok]
    dur = di["duration_s"] - di["masked_s"]
    summ = dict(new_name=r.new_name, holding_command_mV=r.holding_command_mV, internal=r.internal,
                duration_s=di["duration_s"], analysed_s=dur, noise_lp_pA=di["noise_lp_pA"],
                limit_inward_pA=di["limit_inward_pA"], limit_outward_pA=di["limit_outward_pA"],
                n_inward=sum(e["polarity"] == "inward" for e in ev), n_outward=sum(e["polarity"] == "outward" for e in ev))
    summ["freq_inward_Hz"] = summ["n_inward"] / dur
    summ["freq_outward_Hz"] = summ["n_outward"] / dur
    pd.DataFrame(ev).to_csv(cache / (r.new_name[:-4] + ".old.events.csv"), index=False)
    f.write_text(json.dumps(summ, default=float))
    print(f"[{i}/{len(m)}] {r.new_name} {time.time() - ts:.0f}s", flush=True)

names = set(m.new_name)
S = [json.loads(x.read_text()) for x in cache.glob("*.old.json") if x.name[:-9] + ".abf" in names]
if not S:
    sys.exit("nothing new to add")
E = [pd.read_csv(x) for x in cache.glob("*.old.events.csv") if x.name[:-15] + ".abf" in names and x.stat().st_size > 5]
new_f = pd.DataFrame(S)
new_e = pd.concat(E, ignore_index=True) if E else pd.DataFrame(columns=["new_name", "polarity", "peak_s", "amp_pA"])
old_f = em.read_sheet(mp, "alt_det_per_file"); old_e = em.read_sheet(mp, "alt_det_events")
fo = em.merge_frames(old_f, new_f, em.KEYS["alt_det_per_file"], sort=["new_name"])
eo = em.merge_frames(old_e, new_e, em.KEYS["alt_det_events"], sort=["new_name", "peak_s"])
runs = em.read_sheet(mp, "runs")
row = pd.DataFrame([dict(tool="alt_det (second detector)", params_hash="alt_v1", n_recording_directions=len(fo), n_events=len(eo))])
runs = em.merge_frames(runs, row, ["tool", "params_hash"]) if len(runs) else row
em.write_sheets(mp, {"alt_det_per_file": fo, "alt_det_events": eo, "runs": runs})
print(f"master updated: alt_det_per_file {len(fo)} rows, alt_det_events {len(eo)} rows")
