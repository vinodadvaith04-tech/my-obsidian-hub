"""
pick_audit_round2.py - choose the windows of audit round 2 (aimed at what round 1 showed miniML misses).

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
Usage:  python scripts/pick_audit_round2.py [--n-windows 120] [--seed 2] [--out audit_round2_windows.csv]
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
    ap.add_argument("--n-windows", type=int, default=120)
    ap.add_argument("--seed", type=int, default=2)
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

    N = a.n_windows
    plan = [("Cs 17-18 Nov", c[(c.internal == "Cs-gluconate") & c.day.isin(["20261117", "20261118"])], round(N * 0.34)),
            ("Cs earlier", c[(c.internal == "Cs-gluconate") & ~c.day.isin(["20261117", "20261118"])], round(N * 0.25))]
    for tag, pool, n in plan:
        hs = sorted(pool.holding.unique())
        w = np.array([np.sqrt((pool.holding == x).sum()) for x in hs]); share = np.maximum(2, np.round(n * w / w.sum())).astype(int)
        for hv, k in zip(hs, share):
            p = pool[pool.holding == hv]
            take(p[p.n_missed > 0], int(np.ceil(k / 2)), "missed", "n_missed", tag)
            take(p, k // 2, "random", None, tag)
    kg = c[(c.internal == "KGlu") & (c.holding <= -55)]
    take(kg[kg.n_missed > 0], round(N * 0.085), "missed", "n_missed", "K-gluconate negative")
    take(kg, round(N * 0.085), "random", None, "K-gluconate negative")
    nk = c[c.internal != "KCl"]
    take(nk[nk.n_lowfp > 0], round(N * 0.17), "fp", "n_lowfp", "FP hunt")
    take(c[c.internal == "KCl"], round(N * 0.07), "random", None, "KCl")

    out = pd.DataFrame(rows)
    # order so that any first N windows are a mix (interleave strata)
    out["_g"] = out.groupby("stratum").cumcount()
    out = out.sort_values(["_g", "stratum"], kind="stable").drop(columns="_g").reset_index(drop=True)
    path = a.out or (a.master.parent / "audit" / "audit_round2_windows.csv")
    path.parent.mkdir(exist_ok=True); out.to_csv(path, index=False)
    print(len(out), "windows,", out.file.nunique(), "recordings")
    print(out.stratum.str.split(" \\| ").str[0].value_counts().to_string())
    print(out.stratum.str.split(" \\| ").str[-1].value_counts().to_string())
    print("written", path)


if __name__ == "__main__":
    main()
