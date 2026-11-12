"""
build_audit.py - make the hand-audit pages that compare miniML and the old detector on the same 3 s windows.

Usage:
  python build_audit.py --miniml-dir <miniML_results> --old-dir <old_det> --data-root <folder with .abf> --out <audit_folder> [--n-windows 144] [--seed 1]

Output (in --out):
  page_01.png ...            4 windows per page; each window has 2 panels (top: miniML, bottom: old detector) on the SAME trace
  audit_windows.csv          what every window is (file, start, holding, internal, stratum, why chosen)
  audit_events_key.csv       every marker: id (M3 / O5), time, amplitude, score or SNR, whether the other detector also found it (+-10 ms)
  audit_labels.csv           THE FILE YOU FILL IN (one row per window)
Windows are ordered so that ANY first N pages are already a stratified sample (do as many pages as you have time for).
Only the direction(s) allowed by the command holding are shown for BOTH detectors (same rule as the miniML run).
"""
import argparse, math
from pathlib import Path
import numpy as np, pandas as pd, pyabf
from scipy.signal import butter, sosfiltfilt
from scipy.ndimage import median_filter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WIN = 3.0
TOL = 0.010


def hold_class(h):
    if h >= 0: return "positive holding (>= 0 mV)"
    if h >= -50: return "intermediate (-50 to -30 mV)"
    return "negative holding (<= -55 mV)"


def stratum_of(r):
    return f"{r['internal']} | {r['drug_tag']} | {hold_class(r['holding_command_mV'])}"


def drug_tag(r):
    d = str(r.get("drug", "") or "").lower()
    if "gabazine" in d or "sr95531" in d: return "gabazine"
    if "dnqx" in d or "ap5" in d: return "DNQX+AP5"
    if "ttx" in d: return "TTX"
    return "no blocker"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--miniml-dir", type=Path, required=True)
    ap.add_argument("--old-dir", type=Path, required=True)
    ap.add_argument("--data-root", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--n-windows", type=int, default=144)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    ms = pd.read_csv(a.miniml_dir / "miniML_file_summary.csv")
    me = pd.read_csv(a.miniml_dir / "miniML_events_all.csv")
    oe = pd.read_csv(a.old_dir / "oldDetector_events_all.csv")
    # one row per recording
    rec = ms.drop_duplicates("new_name").copy()
    rec["drug_tag"] = rec.apply(drug_tag, axis=1)
    rec["stratum"] = rec.apply(stratum_of, axis=1)
    rec["dirs"] = rec["miniML_directions"]
    rec = rec[rec["duration_s"] > 10]
    me["peak"] = me["peak_time_s"]
    me["pol"] = np.where(me["direction"] == "negative", "down", "up")
    oe["pol"] = np.where(oe["polarity"] == "inward", "down", "up")
    # keep only directions allowed by the holding rule, for both detectors
    allowed = {r.new_name: set(("down" if d == "negative" else "up") for d in str(r.dirs).split(",")) for r in rec.itertuples()}
    oe = oe[[p in allowed.get(f, set()) for f, p in zip(oe["new_name"], oe["pol"])]].copy()
    me_by = {k: v for k, v in me.groupby("new_name")}
    oe_by = {k: v for k, v in oe.groupby("new_name")}

    # disagreement events (found by only one detector, +-10 ms)
    def unmatched(f):
        m = me_by.get(f); o = oe_by.get(f)
        mt = m["peak"].values if m is not None else np.array([]); ot = o["peak_s"].values if o is not None else np.array([])
        out = []
        for t in mt:
            if len(ot) == 0 or np.min(np.abs(ot - t)) > TOL: out.append(t)
        for t in ot:
            if len(mt) == 0 or np.min(np.abs(mt - t)) > TOL: out.append(t)
        return np.array(out)

    # allocate windows to strata (square-root of recordings, at least 4 each)
    strata = rec.groupby("stratum").size()
    w = np.sqrt(strata); alloc = np.maximum(4, np.round(a.n_windows * w / w.sum())).astype(int)
    wins = {}
    for s, n in alloc.items():
        r_s = rec[rec["stratum"] == s]
        lst = []
        n_dis = n // 2
        # disagreement windows
        pool = []
        for r in r_s.itertuples():
            for t in unmatched(r.new_name): pool.append((r.new_name, t, r.duration_s))
        rng.shuffle(pool)
        for f, t, dur in pool[:n_dis]:
            st = float(np.clip(t - WIN * rng.uniform(0.2, 0.8), 0, max(0.1, dur - WIN)))
            lst.append((f, st, "disagreement"))
        while len(lst) < n:
            r = r_s.iloc[rng.integers(len(r_s))]
            st = float(rng.uniform(0, max(0.1, r.duration_s - WIN)))
            lst.append((r.new_name, st, "random"))
        rng.shuffle(lst)
        wins[s] = lst
    # round-robin interleave so any prefix is stratified
    order = []
    keys = list(wins)
    i = 0
    while any(wins[k] for k in keys):
        for k in keys:
            if wins[k]: order.append((k,) + wins[k].pop())
    order = order[: a.n_windows + len(keys)]
    sos_lp = None
    rows, key, labels = [], [], []
    per_page = 4
    npages = math.ceil(len(order) / per_page)
    paths = {p.name: p for p in a.data_root.glob("*.abf")}
    info = rec.set_index("new_name")
    for pg in range(npages):
        chunk = order[pg * per_page:(pg + 1) * per_page]
        fig, axes = plt.subplots(len(chunk) * 2, 1, figsize=(15, 4.3 * len(chunk)), dpi=80,
                                 gridspec_kw=dict(height_ratios=[1] * (len(chunk) * 2)))
        axes = np.atleast_1d(axes)
        for j, (stratum, f, st, why) in enumerate(chunk):
            wid = f"W{pg * per_page + j + 1:03d}"
            ab = pyabf.ABF(str(paths[f])); ab.setSweep(0)
            y = np.asarray(ab.sweepY, float); fs = float(ab.dataRate)
            sos = butter(4, 3000, "low", fs=fs, output="sos")
            i0, i1 = int(st * fs), int((st + WIN) * fs)
            lo, hi = max(0, i0 - int(2 * fs)), min(len(y), i1 + int(2 * fs))
            seg = y[lo:hi]
            lp = sosfiltfilt(sos, seg)
            base = median_filter(lp[::20], size=max(3, int(1.0 * fs / 20)), mode="nearest")
            base = np.interp(np.arange(len(seg)), np.arange(0, len(seg), 20)[:len(base)], base)
            seg_d = seg - base; lp_d = lp - base
            sl = slice(i0 - lo, i1 - lo)
            tt = (np.arange(len(seg)) + lo) / fs
            ylo, yhi = np.percentile(lp_d[sl], 0.2), np.percentile(lp_d[sl], 99.8)
            pad = 0.25 * max(yhi - ylo, 10)
            m = me_by.get(f); o = oe_by.get(f)
            mw = m[(m["peak"] >= st) & (m["peak"] < st + WIN)].sort_values("peak") if m is not None else pd.DataFrame(columns=me.columns)
            ow = o[(o["peak_s"] >= st) & (o["peak_s"] < st + WIN)].sort_values("peak_s") if o is not None else pd.DataFrame(columns=oe.columns)
            ev_amp = np.concatenate([[(lp_d[int(t * fs) - lo] if 0 <= int(t * fs) - lo < len(lp_d) else 0)] for t in list(mw["peak"]) + list(ow["peak_s"])]) if (len(mw) + len(ow)) else np.array([0])
            ylo = min(ylo, ev_amp.min()) - pad; yhi = max(yhi, ev_amp.max()) + pad
            r = info.loc[f]
            for k, (lab, det, tcol) in enumerate([("miniML", mw, "peak"), ("old detector", ow, "peak_s")]):
                ax = axes[j * 2 + k]
                ax.plot(tt[sl], seg_d[sl], color="#c8c8c8", lw=0.5)
                ax.plot(tt[sl], lp_d[sl], color="k", lw=0.8)
                ax.set_xlim(st, st + WIN); ax.set_ylim(ylo, yhi)
                prefix = "M" if k == 0 else "O"
                for n_, (_, e) in enumerate(det.iterrows(), 1):
                    t = float(e[tcol]); idx = int(t * fs) - lo
                    yv = lp_d[idx] if 0 <= idx < len(lp_d) else 0
                    col = "#d6001c" if e["pol"] == "down" else "#1f5fbf"
                    ax.plot([t], [yv], marker="v" if e["pol"] == "down" else "^", color=col, ms=7, mec="k", mew=0.4, ls="none")
                    ax.annotate(f"{n_}", (t, 1.0), xycoords=("data", "axes fraction"), ha="center", va="bottom", fontsize=8, color=col, fontweight="bold")
                    other = oe_by.get(f) if k == 0 else me_by.get(f)
                    ocol = "peak_s" if k == 0 else "peak"
                    matched = bool(other is not None and len(other) and np.min(np.abs(other[ocol].values - t)) <= TOL)
                    key.append(dict(window=wid, id=f"{prefix}{n_}", detector=lab, time_s=t, polarity=e["pol"],
                                    amp_pA=float(e["amplitude"]) if k == 0 else float(e["amp_pA"]),
                                    score_or_snr=float(e["score"]) if k == 0 else float(e["snr"]), other_detector_agrees=matched))
                ax.set_ylabel(f"{lab}\n(pA)", fontsize=9)
                if k == 0:
                    ax.set_title(f"{wid}   {f}   t = {st:.1f}-{st + WIN:.1f} s   holding {r['holding_command_mV']:g} mV   {r['internal']}   [{why}]",
                                 fontsize=10, loc="left", pad=14)
                if k == 1: ax.set_xlabel("time in recording (s)", fontsize=9)
                else: ax.set_xticklabels([])
                for sp in ("top", "right"): ax.spines[sp].set_visible(False)
            rows.append(dict(window=wid, page=pg + 1, file=f, start_s=round(st, 3), end_s=round(st + WIN, 3),
                             holding_mV=r["holding_command_mV"], internal=r["internal"], drug=r["drug_tag"], stratum=stratum,
                             why_chosen=why, directions_shown=r["dirs"], n_miniML=len(mw), n_old=len(ow)))
            labels.append(dict(window=wid, page=pg + 1, file=f, start_s=round(st, 3), miniML_false_ids="", miniML_missed="",
                               old_false_ids="", old_missed="", missed_times_s="", unsure_ids="", note=""))
        fig.tight_layout()
        fig.savefig(a.out / f"page_{pg + 1:02d}.png")
        plt.close(fig)
        print("page", pg + 1, "of", npages, flush=True)
    pd.DataFrame(rows).to_csv(a.out / "audit_windows.csv", index=False)
    pd.DataFrame(key).to_csv(a.out / "audit_events_key.csv", index=False)
    pd.DataFrame(labels).to_csv(a.out / "audit_labels.csv", index=False)
    print("windows", len(rows), "pages", npages, "events in key", len(key))


if __name__ == "__main__":
    main()
