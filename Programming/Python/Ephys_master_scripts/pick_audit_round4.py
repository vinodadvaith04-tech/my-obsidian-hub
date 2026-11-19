"""
pick_audit_round4.py - choose the windows of audit round 4: the weakest recall spots (Cs at intermediate holdings), K-gluconate false events, and a fresh random test set.

Round 1 (188 windows) showed:  miniML is precise (96 %) but misses ~20 % of real events, mostly IPSC-like events
(blunter rise, long decay) that the model scores near 0, in Cs- and K-gluconate recordings; KCl recordings are fine.
So round 2 spends the effort where new labels help the model most:

  Cs, 17-18 Nov (slices 5-6, main experiment)  all holdings, half "second detector sees events miniML does not", half random
  Cs, earlier days (4-12 Nov)                  same, windows not in round 1
  K-gluconate at <= -55 mV                     (EPSCs; round 1 recall only ~0.45-0.6), half missed-shape, half random
  false-positive hunt                          non-KCl windows where miniML marks low-score events nobody else sees
  KCl                                          a few random windows as a stability check

Windows never overlap each other or a round-1 window of the same recording (>= 3 s apart); at most 2 per recording.
Writes a csv (file, start_s, why_chosen, stratum) for:
    python scripts/build_audit.py --windows-from <csv> --batch "4: round 2"
Usage:  python scripts/pick_audit_round2.py [--n-windows 120] [--seed 2] [--out audit_round4_windows.csv]
Needs: pandas numpy openpyxl
"""
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ephys_master as em

WIN, TOL, STEP = 3.0, 0.010, 1.0
MIN_DUR = 30.0
WHY = {"missed": "disagreement", "fp": "disagreement", "random": "random"}   # keeps the labels the tool/scoring already know


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", type=Path, default=em.default_master())
    ap.add_argument("--audit-dir", type=Path, default=None)
    ap.add_argument("--n-windows", type=int, default=80)
    ap.add_argument("--seed", type=int, default=4)
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    audit_dir = a.audit_dir or (a.master.parent / "audit")

    runs = em.read_sheet(a.master, "runs")
    h = runs[runs.tool == "miniML"].sort_values("last_processed").params_hash.iloc[-1]
    pf = em.read_sheet(a.master, "miniML_per_file"); pf = pf[pf.params_hash == h].drop_duplicates("new_name")
    me = em.read_sheet(a.master, "miniML_events"); me = me[me.params_hash == h]
    oe = em.read_sheet(a.master, "alt_det_events")
    pf = pf[pf.new_name.isin(set(em.read_sheet(a.master, "alt_det_per_file").new_name)) & (pf.duration_s >= MIN_DUR)].copy()
    pf["day"] = pf.new_name.str[:8]
    pf["allowed"] = [{("down" if d.strip() == "negative" else "up") for d in str(x).split(",")} for x in pf.miniML_directions]
    me = me.assign(pol=np.where(me.direction == "negative", "down", "up"))
    oe = oe.assign(pol=np.where(oe.polarity == "inward", "down", "up"))
    me_by = {k: v for k, v in me.groupby("new_name")}; oe_by = {k: v for k, v in oe.groupby("new_name")}

    done = {}
    xl = audit_dir / "audit_labels.xlsx"
    if xl.exists():
        w = pd.read_excel(xl, sheet_name="audit_windows")
        for f, s in zip(w.file, w.start_s): done.setdefault(f, []).append(float(s))

    def tables(r):
        m = me_by.get(r.new_name); o = oe_by.get(r.new_name)
        m = m[m.pol.isin(r.allowed)] if m is not None else pd.DataFrame(columns=me.columns)
        o = o[o.pol.isin(r.allowed)] if o is not None else pd.DataFrame(columns=oe.columns)
        mt, ot = m.peak_time_s.values, o.peak_s.values
        o_only = np.array([t for t in ot if len(mt) == 0 or np.min(np.abs(mt - t)) > TOL])
        m_only = m[np.array([len(ot) == 0 or np.min(np.abs(ot - t)) > TOL for t in mt], dtype=bool)]
        m_low = m_only[m_only.score < 0.8].peak_time_s.values
        return mt, ot, o_only, m_low

    def cnt(x, s): return int(((x >= s) & (x < s + WIN)).sum()) if len(x) else 0

    cands = []   # one row per candidate window
    for r in pf.itertuples():
        mt, ot, o_only, m_low = tables(r)
        for s in np.arange(0, r.duration_s - WIN, STEP):
            if any(abs(s - t) < WIN for t in done.get(r.new_name, [])): continue
            cands.append(dict(file=r.new_name, start_s=float(s), internal=r.internal, holding=r.holding_command_mV, day=r.day,
                              cell=r.new_name.rsplit("_", 1)[0], n_m=cnt(mt, s), n_o=cnt(ot, s), n_missed=cnt(o_only, s), n_lowfp=cnt(m_low, s)))
    c = pd.DataFrame(cands)
    print(len(c), "candidate windows in", c.file.nunique(), "recordings")

    taken = {}   # file -> list of starts
    rows = []

    def free(f, s):
        return len(taken.get(f, [])) < 2 and all(abs(s - t) >= WIN for t in taken.get(f, []))

    def take(df, n, kind, score=None, tag=""):
        if n <= 0 or df.empty: return 0
        df = df.copy()
        df["_s"] = (df[score] if score else 0) + rng.uniform(0, 0.5, len(df))
        got = 0
        # spread over cells: round-robin so one cell does not take everything
        by_cell = {k: g.sort_values("_s", ascending=False) for k, g in df.groupby("cell")}
        order = list(by_cell); rng.shuffle(order)
        idx = {k: 0 for k in order}
        while got < n and order:
            for k in list(order):
                g = by_cell[k]
                while idx[k] < len(g) and not free(g.iloc[idx[k]].file, g.iloc[idx[k]].start_s): idx[k] += 1
                if idx[k] >= len(g): order.remove(k); continue
                r = g.iloc[idx[k]]; idx[k] += 1
                taken.setdefault(r.file, []).append(r.start_s)
                rows.append(dict(file=r.file, start_s=round(r.start_s, 3), why_chosen=WHY[kind],
                                 stratum=f"{tag} | {r.internal} | {r.holding:g} mV | {kind}"))
                got += 1
                if got >= n: break
        return got

    N = a.n_windows   # 80: 28 weak-spot Cs + 10 Cs disagreement + 16 K-gluconate false-event hunt + 26 random
    audited_cells = {f.rsplit("_", 1)[0] for f in done}
    c["fresh"] = ~c.cell.isin(audited_cells)
    c["act"] = c.n_m + c.n_o
    def pool(df):
        f = df[df.fresh]
        return f if f.cell.nunique() >= 4 else df
    cs = c[c.internal == "Cs-gluconate"]
    # (A) holdings where round 1-3 showed the weakest recall: inward events at -50..-20 mV, outward at -30 mV, and +15 mV
    for hold, n in ((-50, 5), (-40, 5), (-30, 6), (-20, 4), (0, 2), (10, 2), (15, 4)):
        p = cs[(cs.holding == hold) & (cs.act >= 3)]
        take(p, n, "missed", "act", f"Weak spot {hold:g} mV")
    # (B) Cs windows where the second detector marks events that miniML does not
    take(pool(cs[cs.n_missed >= 2]), 10, "missed", "n_missed", "Detectors disagree")
    # (C) K-gluconate: miniML-only low-score events (38 % of these were false in round 3)
    kg = c[c.internal == "KGlu"]
    take(pool(kg[kg.n_lowfp > 0]), 16, "fp", "n_lowfp", "K-gluconate FP hunt")
    # (D) random windows = unbiased fresh test set
    for internal, n in (("Cs-gluconate", 14), ("KGlu", 8), ("KCl", 4)):
        take(pool(c[c.internal == internal]), n, "random", None, "Random")

    out = pd.DataFrame(rows)
    # order so that any first N windows are a mix (interleave strata)
    out["_g"] = out.groupby("stratum").cumcount()
    out = out.sort_values(["_g", "stratum"], kind="stable").drop(columns="_g").reset_index(drop=True)
    path = a.out or (a.master.parent / "audit" / "audit_round4_windows.csv")
    path.parent.mkdir(exist_ok=True); out.to_csv(path, index=False)
    print(len(out), "windows,", out.file.nunique(), "recordings")
    print(out.stratum.str.split(" \\| ").str[0].value_counts().to_string())
    print(out.stratum.str.split(" \\| ").str[-1].value_counts().to_string())
    print("written", path)


if __name__ == "__main__":
    main()
