"""
audit_recall_by_stratum.py - recall / precision of miniML (and the second detector) from the RANDOM audit windows,
per internal x holding class x polarity, with cell-clustered bootstrap CIs.

Usage: python audit_recall_by_stratum.py --audit-xlsx audit_labels.xlsx --out out_dir [--batch round6] [--n-boot 5000]

Truth = events in audit_truth_events (accepted miniML markers + accepted second-detector markers + your clicks): a lower bound.
Only windows you finished (not "skip") count; only windows with why_chosen == "random" are treated as unbiased.
Unit for the bootstrap = the cell (file stem without the recording counter).
"""
import argparse, re
from pathlib import Path
import numpy as np, pandas as pd

rng = np.random.default_rng(7)


def hclass(h): return "pos (>=0)" if h >= 0 else ("neg (<=-55)" if h <= -55 else "int (-50..-30)")


def ids(s):
    return [] if (s is None or (isinstance(s, float) and np.isnan(s)) or str(s).strip() == "") else [x for x in re.split(r"[;,]", str(s)) if x.strip()]


def load(xl, batch=None):
    lab = pd.read_excel(xl, sheet_name="audit_labels", dtype=str).fillna("")
    win = pd.read_excel(xl, sheet_name="audit_windows")
    tr = pd.read_excel(xl, sheet_name="audit_truth_events")
    f = ["miniML_false_ids", "old_false_ids", "miniML_missed", "old_missed"]
    done = set(lab[(lab[f] != "").any(axis=1)].window)
    win = win[win.window.isin(done) & (win.why_chosen == "random")].copy()
    if batch: win = win[win.batch == batch]
    win["cell"] = win.file.str.replace(r"_\d+\.abf$", "", regex=True)
    win["hclass"] = win.holding_mV.map(hclass)
    win["stratum"] = win.internal + " | " + win.hclass
    lab = lab.set_index("window")
    win["fp_m"] = win.window.map(lambda w: len(ids(lab.loc[w, "miniML_false_ids"])))
    win["un_m"] = win.window.map(lambda w: len(ids(lab.loc[w, "unsure_ids"])))
    win["fp_o"] = win.window.map(lambda w: len(ids(lab.loc[w, "old_false_ids"])))
    tr = tr[tr.window.isin(win.window)].copy()
    # only events in a direction miniML was actually run on in that window (negative holding -> inward=down, positive -> outward=up)
    shown = win.set_index("window").directions_shown.astype(str)
    pol_ok = {"down": "negative", "up": "positive"}
    keep = [pol_ok.get(p, "") in shown.get(w, "") for w, p in zip(tr.window, tr.polarity)]
    print(f"events outside the directions miniML was run on (not counted as misses): {len(tr) - sum(keep)}")
    tr = tr[keep].copy()
    tr["amp_abs"] = tr.amp_pA.abs()
    return win, tr


def boot(df, stat, n_boot):
    """cluster bootstrap over cells; df has column cell; stat(df)->float"""
    cells = df.cell.unique(); groups = {c: g for c, g in df.groupby("cell")}
    out = []
    for _ in range(n_boot):
        pick = rng.choice(cells, len(cells), replace=True)
        out.append(stat(pd.concat([groups[c] for c in pick])))
    out = np.array(out, float)
    return np.nanpercentile(out, 2.5), np.nanpercentile(out, 97.5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit-xlsx", type=Path, required=True); ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--batch", default=None); ap.add_argument("--n-boot", type=int, default=3000)
    a = ap.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    win, tr = load(a.audit_xlsx, a.batch)
    tr = tr.drop(columns=["internal"]).merge(win[["window", "cell", "stratum", "internal", "hclass"]], on="window")
    rows = []
    for pol in ["all", "down", "up"]:
        t = tr if pol == "all" else tr[tr.polarity == pol]
        for st, g in list(t.groupby("stratum")) + [("ALL", t)]:
            if st == "ALL": w = win
            else: w = win[win.stratum == st]
            rec = lambda d: d.found_by_miniML.mean() if len(d) else np.nan
            reco = lambda d: d.found_by_old.mean() if len(d) else np.nan
            lo, hi = boot(g, rec, a.n_boot) if len(g) > 5 and g.cell.nunique() > 1 else (np.nan, np.nan)
            rows.append(dict(polarity=pol, stratum=st, windows=len(w), cells=w.cell.nunique(), truth_events=len(g),
                             truth_per_window=len(g) / max(len(w), 1), miniML_found=int(g.found_by_miniML.sum()),
                             miniML_recall=rec(g), recall_lo=lo, recall_hi=hi, old_recall=reco(g)))
    res = pd.DataFrame(rows)
    # precision (miniML): markers in random windows
    p = win.groupby("stratum").agg(windows=("window", "size"), n_miniML=("n_miniML", "sum"), fp=("fp_m", "sum"), unsure=("un_m", "sum"), n_old=("n_old", "sum"), fp_old=("fp_o", "sum"))
    p["miniML_precision"] = 1 - p.fp / (p.n_miniML - p.unsure).clip(lower=1)
    p["old_precision"] = 1 - p.fp_old / p.n_old.clip(lower=1)
    # recall by amplitude (all strata pooled and by internal)
    bins = [0, 10, 15, 20, 30, 1e9]; lab_ = ["<10", "10-15", "15-20", "20-30", ">30"]
    tr["amp_bin"] = pd.cut(tr.amp_abs, bins, labels=lab_)
    amp = tr.groupby(["internal", "amp_bin"], observed=True).agg(truth=("window", "size"), found=("found_by_miniML", "sum")).reset_index()
    amp["recall"] = amp.found / amp.truth
    suffix = f"_{a.batch}" if a.batch else "_all_random"
    res.to_csv(a.out / f"recall_by_stratum{suffix}.csv", index=False)
    p.to_csv(a.out / f"precision_by_stratum{suffix}.csv")
    amp.to_csv(a.out / f"recall_by_amplitude{suffix}.csv", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print(res[res.polarity == "all"].round(3).to_string(index=False))
    print(res[res.polarity != "all"].round(3).to_string(index=False))
    print(p.round(3).to_string()); print(amp.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
