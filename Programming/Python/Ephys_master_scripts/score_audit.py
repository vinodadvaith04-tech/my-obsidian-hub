"""
score_audit.py - score the filled-in audit_labels.xlsx (one workbook, sheets audit_labels / audit_windows / audit_events_key) for miniML and the old detector.

Usage:  python scripts/score_audit.py [--audit <audit_folder>]    (default: the audit folder next to scripts/)
Only windows where you filled in at least one of the four count/id columns are scored
(an empty row = not audited yet; to say "all fine" write 0 in miniML_missed and old_missed).

Definitions (per detector, over the audited windows):
  false positive (FP)  = marker id listed in *_false_ids
  unsure               = id listed in unsure_ids (left out of precision)
  true positive (TP)   = every other marker
  missed               = real events without a marker (count in *_missed)
  precision = TP/(TP+FP)    recall = TP/(TP+missed)    F1
Also written: audit_score.txt, audit_events_labelled.xlsx (every marker with TP/FP/unsure; usable as training labels).
"""
import argparse
from pathlib import Path
import pandas as pd, numpy as np


def ids(s):
    if pd.isna(s) or str(s).strip() == "": return []
    return [x.strip() for x in str(s).replace(",", ";").split(";") if x.strip()]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--audit", type=Path, default=Path(__file__).resolve().parent.parent / "audit")
    a = ap.parse_args()
    xl = a.audit / "audit_labels.xlsx"
    lab = pd.read_excel(xl, sheet_name="audit_labels", dtype=str).fillna("")
    keep = ["window", "page", "batch", "file", "start_s", "miniML_false_ids", "miniML_missed", "old_false_ids", "old_missed", "missed_times_s", "unsure_ids", "note"]
    lab = lab[[c for c in keep if c in lab.columns]]
    lab["miniML_missed"] = lab["miniML_missed"].astype(str).str.replace(r"\.0$", "", regex=True)
    lab["old_missed"] = lab["old_missed"].astype(str).str.replace(r"\.0$", "", regex=True)
    win = pd.read_excel(xl, sheet_name="audit_windows")
    key = pd.read_excel(xl, sheet_name="audit_events_key")
    done = lab[(lab[["miniML_false_ids", "old_false_ids", "miniML_missed", "old_missed"]].notna()
                & (lab[["miniML_false_ids", "old_false_ids", "miniML_missed", "old_missed"]] != "")).any(axis=1)].copy()
    print(f"{len(done)} of {len(lab)} windows audited")
    fp = set(); un = set()
    for r in done.itertuples():
        for i in ids(r.miniML_false_ids): fp.add((r.window, "M" + i.lstrip("MmOo")))
        for i in ids(r.old_false_ids): fp.add((r.window, "O" + i.lstrip("MmOo")))
        for i in ids(r.unsure_ids): un.add((r.window, i.upper() if i[:1].upper() in "MO" else None))
    k = key[key["window"].isin(done["window"])].copy()
    k["label"] = ["FP" if (w, i) in fp else ("unsure" if (w, i) in un else "TP") for w, i in zip(k["window"], k["id"])]
    k = k.merge(win[["window", "batch", "stratum", "internal", "holding_mV"]], on="window", how="left")
    k["amp_abs"] = k["amp_pA"].abs()
    k.to_excel(a.audit / "audit_events_labelled.xlsx", index=False)
    miss = {"miniML": done.set_index("window")["miniML_missed"].apply(lambda x: 0 if pd.isna(x) or x == "" else int(float(x))),
            "old detector": done.set_index("window")["old_missed"].apply(lambda x: 0 if pd.isna(x) or x == "" else int(float(x)))}
    out = []

    def line(name, kk, mm):
        tp = (kk["label"] == "TP").sum(); f = (kk["label"] == "FP").sum()
        p = tp / (tp + f) if tp + f else np.nan
        if mm is None:
            return f"{name:48s} n_marked={len(kk):4d}  TP={tp:4d} FP={f:4d}  precision={p:.3f}  (recall not split by amplitude/internal-free)"
        r = tp / (tp + mm) if tp + mm else np.nan
        f1 = 2 * p * r / (p + r) if p and r else np.nan
        return f"{name:48s} n_marked={len(kk):4d}  TP={tp:4d} FP={f:4d} missed={mm:4d}  precision={p:.3f} recall={r:.3f} F1={f1:.3f}"
    out.append(f"AUDIT SCORE  ({len(done)} windows)\n")
    for det in ("miniML", "old detector"):
        kk = k[k["detector"] == det]; mm = int(miss[det].sum())
        out.append(line(det + " - overall", kk, mm))
        for s, g in kk.groupby("internal"):
            w = done.merge(win[["window", "internal"]], on="window")
            wi = w[w["internal"] == s]["window"]
            out.append(line(f"   {det} - {s}", g, int(miss[det].reindex(wi).sum())))
        for lo, hi, nm in [(0, 10, "<10 pA"), (10, 20, "10-20 pA"), (20, 1e9, ">20 pA")]:
            g = kk[(kk["amp_abs"] >= lo) & (kk["amp_abs"] < hi)]
            out.append(line(f"   {det} - amplitude {nm}", g, None))
        for ag, nm in [(True, "also found by the other detector"), (False, "found only by this detector")]:
            g = kk[kk["other_detector_agrees"] == ag]
            tp = (g["label"] == "TP").sum(); f = (g["label"] == "FP").sum()
            out.append(f"   {det} - {nm}: n={len(g)}  precision={tp / (tp + f) if tp + f else float('nan'):.3f}")
        out.append("")
    txt = "\n".join(out)
    (a.audit / "audit_score.txt").write_text(txt, encoding="utf-8")
    print(txt)


if __name__ == "__main__":
    main()
