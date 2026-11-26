"""
audit_correction_factors.py - correction factors from the RANDOM audit windows: how many real events there are per miniML marker,
per internal x holding class x direction, with cell-clustered bootstrap CIs.

   corrected events = miniML markers x factor          (factor = truth events / miniML markers in the audited random windows)
   factor = precision / recall   (precision < 1 pulls it down, recall < 1 pushes it up)

Usage: python audit_correction_factors.py --audit-xlsx audit_labels.xlsx --out out_dir [--n-boot 4000]
Truth is a lower bound (accepted miniML + accepted second-detector markers + your clicks). Skipped windows are excluded.
"""
import argparse, re
from pathlib import Path
import numpy as np, pandas as pd

rng = np.random.default_rng(11)
hclass = lambda h: "pos (>=0)" if h >= 0 else ("neg (<=-55)" if h <= -55 else "int (-50..-30)")


def ids(s):
    return [] if (s is None or (isinstance(s, float) and np.isnan(s)) or str(s).strip() == "") else [x for x in re.split(r"[;,]", str(s)) if x.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit-xlsx", type=Path, required=True); ap.add_argument("--out", type=Path, required=True); ap.add_argument("--n-boot", type=int, default=4000)
    a = ap.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    x = a.audit_xlsx
    lab = pd.read_excel(x, sheet_name="audit_labels", dtype=str).fillna("")
    win = pd.read_excel(x, sheet_name="audit_windows"); key = pd.read_excel(x, sheet_name="audit_events_key"); tr = pd.read_excel(x, sheet_name="audit_truth_events")
    f = ["miniML_false_ids", "old_false_ids", "miniML_missed", "old_missed"]
    done = set(lab[(lab[f] != "").any(axis=1)].window)
    win = win[win.window.isin(done) & (win.why_chosen == "random")].copy()
    win["cell"] = win.file.str.replace(r"_\d+\.abf$", "", regex=True); win["hclass"] = win.holding_mV.map(hclass); win["stratum"] = win.internal + " | " + win.hclass
    shown = win.set_index("window").directions_shown.astype(str); ok = {"down": "negative", "up": "positive"}
    # per window x polarity: miniML markers (in directions miniML was run on) and truth events
    mk = key[key.window.isin(win.window) & key.id.astype(str).str.startswith("M")]
    mk = mk[[ok.get(p, "") in shown[w] for w, p in zip(mk.window, mk.polarity)]]
    t = tr[tr.window.isin(win.window)]; t = t[[ok.get(p, "") in shown[w] for w, p in zip(t.window, t.polarity)]]
    rows = []
    for pol in ["down", "up"]:
        d = win[["window", "cell", "stratum", "internal"]].copy()
        d["markers"] = d.window.map(mk[mk.polarity == pol].groupby("window").size()).fillna(0)
        d["truth"] = d.window.map(t[t.polarity == pol].groupby("window").size()).fillna(0)
        d["tp"] = d.window.map(t[(t.polarity == pol) & t.found_by_miniML].groupby("window").size()).fillna(0)
        d["polarity"] = pol
        d = d[[ok[pol] in shown[w] for w in d.window]]          # windows where this direction was run
        rows.append(d)
    d = pd.concat(rows)
    out = []
    for (st, pol), g in list(d.groupby(["stratum", "polarity"])) + [(("ALL Cs-gluconate", pol), gg) for pol, gg in d[d.internal == "Cs-gluconate"].groupby("polarity")] \
            + [(("ALL K-Glu", pol), gg) for pol, gg in d[d.internal == "KGlu"].groupby("polarity")]:
        if g.markers.sum() == 0: continue
        cells = g.cell.unique(); grp = {c: gg for c, gg in g.groupby("cell")}
        def ratio(df): return df.truth.sum() / df.markers.sum() if df.markers.sum() else np.nan
        bs = np.array([ratio(pd.concat([grp[c] for c in rng.choice(cells, len(cells))])) for _ in range(a.n_boot)])
        out.append(dict(stratum=st, direction=pol, windows=len(g), cells=len(cells), miniML_markers=int(g.markers.sum()), truth_events=int(g.truth.sum()),
                        factor=ratio(g), lo=np.nanpercentile(bs, 2.5), hi=np.nanpercentile(bs, 97.5), recall=g.tp.sum() / g.truth.sum() if g.truth.sum() else np.nan))
    r = pd.DataFrame(out); r.to_csv(a.out / "correction_factors.csv", index=False)
    pd.set_option("display.width", 220); print(r.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
