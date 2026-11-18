"""
pick_audit_windows.py - choose audit windows from the recordings of one session, to be added with build_audit.py.

Picks a spread of recordings (low / middle / high event rate, different days) and, in each, three 3-s windows that do not overlap:
  dense         the window with the most events (both detectors together)  -> bursts, overlapping events
  disagreement  the window where miniML and the second detector disagree most
  random        a random window
Writes a csv (file, start_s, why_chosen, stratum) for:   python scripts/build_audit.py --windows-from <csv> --batch "<name>"

Usage:  python scripts/pick_audit_windows.py --session "Helena KCl IPSC" [--n-recordings 6] [--seed 1] [--out windows.csv]
Needs miniML and second-detector results for the session in the master.
"""
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ephys_master as em

WIN, TOL = 3.0, 0.010


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", type=Path, default=em.default_master())
    ap.add_argument("--session", required=True)
    ap.add_argument("--n-recordings", type=int, default=6)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    files = em.read_sheet(a.master, "files")
    runs = em.read_sheet(a.master, "runs")
    h = runs[runs.tool == "miniML"].sort_values("last_processed").params_hash.iloc[-1]
    pf = em.read_sheet(a.master, "miniML_per_file"); pf = pf[(pf.params_hash == h) & (pf.session == a.session)]
    me = em.read_sheet(a.master, "miniML_events"); me = me[(me.params_hash == h) & me.new_name.isin(pf.new_name)]
    oe = em.read_sheet(a.master, "alt_det_events"); oe = oe[oe.new_name.isin(pf.new_name)]
    done_alt = set(em.read_sheet(a.master, "alt_det_per_file").new_name)
    pf = pf[pf.new_name.isin(done_alt)].drop_duplicates("new_name").copy()
    if not len(pf): sys.exit("no recordings with both detectors for this session")
    # spread: sort by miniML rate, take evenly spaced recordings (keeps a mix of low and high rates), at most one per cell
    pf["cellkey"] = pf.new_name.str.replace(r"_\d{3}\.abf$", "", regex=True)
    pf = pf.sort_values("frequency_Hz").drop_duplicates("cellkey")
    idx = np.unique(np.round(np.linspace(0, len(pf) - 1, min(a.n_recordings, len(pf)))).astype(int))
    chosen = pf.iloc[idx]
    rows = []
    for r in chosen.itertuples():
        dur = float(r.duration_s)
        m = me[me.new_name == r.new_name]["peak_time_s"].values
        direction = set(me[me.new_name == r.new_name]["direction"]) or {"negative"}
        pol = "inward" if "negative" in direction else "outward"
        o = oe[(oe.new_name == r.new_name) & (oe.polarity == pol)]["peak_s"].values
        starts = np.arange(0, dur - WIN, 0.5)
        def count(x, s): return int(((x >= s) & (x < s + WIN)).sum()) if len(x) else 0
        def unm(x, y):   # events of x with no event of y within TOL
            return np.array([t for t in x if len(y) == 0 or np.min(np.abs(y - t)) > TOL])
        mu, ou = unm(m, o), unm(o, m)
        dense = np.array([count(m, s) + count(o, s) for s in starts])
        dis = np.array([count(mu, s) + count(ou, s) for s in starts])
        taken = []
        def pick(score, why):
            order = np.argsort(-score, kind="stable")
            for k in order:
                s = starts[k]
                if all(abs(s - t) >= WIN for t in taken):
                    taken.append(s); rows.append(dict(file=r.new_name, start_s=round(float(s), 3), why_chosen=why,
                                                      stratum=f"{r.internal} | {r.new_name.split('_')[0]} | rate {r.frequency_Hz:.2f} Hz")); return
        pick(dense + rng.uniform(0, 0.01, len(dense)), "dense")
        pick(dis + rng.uniform(0, 0.01, len(dis)), "disagreement")
        pick(rng.uniform(0, 1, len(starts)), "random")
    out = pd.DataFrame(rows)
    path = a.out or (a.master.parent / "_work" / "audit_windows_pick.csv")
    path.parent.mkdir(exist_ok=True); out.to_csv(path, index=False)
    print(chosen[["new_name", "frequency_Hz"]].to_string()); print(out.to_string()); print("written", path)


if __name__ == "__main__":
    main()
