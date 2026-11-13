"""
build_audit.py - hand-audit pages: miniML (top panel) and the second detector (bottom panel) on the SAME 3 s windows.

ONE audit folder, ONE workbook:
  <audit>/page_01.png, page_02.png, ...   4 windows per page, ids W001... in the panel titles
  <audit>/audit_labels.xlsx               sheets: audit_labels (YOU fill the pink columns), audit_windows, audit_events_key, how_to

Everything is read from Ephys_master.xlsx (miniML_events, alt_det_events, files). Running it again APPENDS: new
windows get the next W-numbers and page numbers, labels you already entered are kept untouched, and a window
(file + start) that is already in the audit is never added twice. The audit workbook is backed up first.

Usage (examples)
  python scripts/build_audit.py --n-windows 26 --select-prefix 20261112_ --holdings -80 -50 -40 -30 -20 10 \\
         --min-duration 30 --stratum-by holding --batch "Cs slice 4"
  python scripts/build_audit.py --windows-from some_windows_table.csv --batch "re-render"      (exact windows given)
Options: --master, --data-root, --audit-dir (default: <master folder>/audit), --params-hash (default: the miniML run in `runs`),
         --n-windows, --seed, --min-duration, --select-prefix, --holdings, --stratum-by {class,holding}, --batch
Windows are ordered so that ANY first N pages are already a stratified sample.
Only the direction(s) allowed by the command holding are shown for BOTH detectors.
"""
import argparse, math, shutil, sys, time
from pathlib import Path
import numpy as np, pandas as pd, pyabf
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter
from scipy.signal import butter, sosfiltfilt
from scipy.ndimage import median_filter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ephys_master as em

WIN, TOL, PER_PAGE = 3.0, 0.010, 4
LABEL_COLS = ["miniML_false_ids", "miniML_missed", "old_false_ids", "old_missed", "missed_times_s", "unsure_ids", "note"]
INFO_COLS = ["window", "page", "batch", "holding_mV", "internal", "drug", "why_chosen", "n_miniML", "n_old", "file", "start_s"]
HOW_TO = [
    "Fill the pink columns, one row per window (page_XX.png has the 4 windows of that page, W-numbers in the panel titles).",
    "miniML_false_ids / old_false_ids: numbers of markers that are NOT real events, e.g. 3;7 (leave empty if none).",
    "miniML_missed / old_missed: number of real events that have no marker in that panel. Write 0 if none.",
    "unsure_ids: markers you cannot call, with prefix, e.g. M3;O2.",
    "Unusable window: leave the counts empty and write 'skip' in note.",
    "A row with all four id/count columns empty counts as not audited.",
    "'old' in the column names = the second detector (matched filter), shown in the bottom panel. Top panel = miniML (M-ids).",
    "Keep the file name audit_labels.xlsx and close it in Excel before running build_audit.py or score_audit.py.",
    "batch = which build added the window. New windows are appended to the end; your labels are never touched.",
]


def hold_class(h):
    if h >= 0: return "positive holding (>= 0 mV)"
    if h >= -50: return "intermediate (-50 to -30 mV)"
    return "negative holding (<= -55 mV)"


def drug_tag(d):
    d = str(d if isinstance(d, str) else "").lower()
    if "gabazine" in d or "sr95531" in d: return "gabazine"
    if "dnqx" in d or "ap5" in d: return "DNQX+AP5"
    if "ttx" in d: return "TTX"
    return "no blocker"


# ---------------------------------------------------------------- workbook of the audit
def _write_audit(path: Path, labels, windows, key):
    wb = Workbook(); ws = wb.active; ws.title = "audit_labels"
    head_fill = PatternFill("solid", fgColor="730E04"); pink = PatternFill("solid", fgColor="FDECEA")

    def put(ws, df, widths=None, pinkcols=()):
        ws.append(list(df.columns))
        for r in df.itertuples(index=False):
            ws.append([None if (isinstance(v, float) and np.isnan(v)) or v is pd.NA else (v.item() if hasattr(v, "item") else v) for v in r])
        for c in ws[1]:
            c.font = Font(bold=True, color="FFFFFF"); c.fill = head_fill
        ws.freeze_panes = "C2"
        ws.auto_filter.ref = f"A1:{get_column_letter(len(df.columns))}{len(df) + 1}"
        for j, col in enumerate(df.columns, 1):
            w = max([len(str(col))] + [len(str(v)) for v in df[col].head(200).tolist()])
            ws.column_dimensions[get_column_letter(j)].width = min(max(w + 2, 8), 46)
            if col in pinkcols:
                for i in range(2, len(df) + 2):
                    ws.cell(i, j).fill = pink
    put(ws, labels, pinkcols=LABEL_COLS)
    put(wb.create_sheet("audit_windows"), windows)
    put(wb.create_sheet("audit_events_key"), key)
    h = wb.create_sheet("how_to")
    for t in HOW_TO: h.append([t])
    h.column_dimensions["A"].width = 125
    tmp = path.with_suffix(".tmp.xlsx")
    wb.save(tmp); tmp.replace(path)


def _read_audit(path: Path):
    if not path.exists():
        return None
    x = pd.ExcelFile(path)
    return {s: x.parse(s, dtype=object) for s in ("audit_labels", "audit_windows", "audit_events_key")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", type=Path, default=em.default_master())
    ap.add_argument("--data-root", type=Path, default=None)
    ap.add_argument("--audit-dir", type=Path, default=None)
    ap.add_argument("--params-hash", default=None)
    ap.add_argument("--n-windows", type=int, default=144)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--min-duration", type=float, default=10)
    ap.add_argument("--select-prefix", default=None, help="only recordings whose new_name starts with this (e.g. 20261112_)")
    ap.add_argument("--holdings", type=float, nargs="*", default=None)
    ap.add_argument("--stratum-by", choices=["class", "holding"], default="class")
    ap.add_argument("--batch", default="")
    ap.add_argument("--windows-from", type=Path, default=None, help="csv/xlsx with columns file,start_s[,why_chosen,stratum]: use exactly these windows")
    ap.add_argument("--first-window", type=int, default=None)
    ap.add_argument("--first-page", type=int, default=None)
    a = ap.parse_args()
    root = a.data_root or a.master.parent
    out = a.audit_dir or (a.master.parent / "audit")
    out.mkdir(parents=True, exist_ok=True)
    xl = out / "audit_labels.xlsx"
    rng = np.random.default_rng(a.seed)

    # ---- inputs from the master
    runs = em.read_sheet(a.master, "runs")
    h = a.params_hash or runs[runs.tool == "miniML"].sort_values("last_processed").params_hash.iloc[-1]
    pf = em.read_sheet(a.master, "miniML_per_file"); pf = pf[pf.params_hash == h]
    me = em.read_sheet(a.master, "miniML_events"); me = me[me.params_hash == h].copy()
    oe = em.read_sheet(a.master, "alt_det_events").copy()
    ap_f = set(em.read_sheet(a.master, "alt_det_per_file").new_name)
    print(f"master {a.master.name}: miniML run {h}: {len(pf)} recording-directions, {len(me)} events; second detector {len(ap_f)} recordings, {len(oe)} events")

    rec = pf.drop_duplicates("new_name").copy()
    rec = rec[rec.new_name.isin(ap_f)]
    wf = None
    if a.windows_from:
        wf = pd.read_csv(a.windows_from) if a.windows_from.suffix == ".csv" else pd.read_excel(a.windows_from)
        rec = rec[rec.new_name.isin(wf["file"])]
    else:
        rec = rec[rec.duration_s > a.min_duration]
        if a.select_prefix: rec = rec[rec.new_name.str.startswith(a.select_prefix)]
        if a.holdings: rec = rec[rec.holding_command_mV.isin(a.holdings)]
    rec["drug_tag"] = rec["drug"].map(drug_tag)
    rec["dirs"] = rec["miniML_directions"]
    rec["stratum"] = [f"{r.internal} | {r.drug_tag} | {hold_class(r.holding_command_mV)}" if a.stratum_by == "class"
                      else f"{r.internal} | {r.holding_command_mV:g} mV" for r in rec.itertuples()]
    rec = rec.sort_values("new_name").reset_index(drop=True)
    print(f"{len(rec)} recordings eligible for the audit")
    me["peak"] = me["peak_time_s"]
    me["pol"] = np.where(me["direction"] == "negative", "down", "up")
    oe["pol"] = np.where(oe["polarity"] == "inward", "down", "up")
    allowed = {r.new_name: {("down" if d.strip() == "negative" else "up") for d in str(r.dirs).split(",")} for r in rec.itertuples()}
    oe = oe[[p in allowed.get(f, set()) for f, p in zip(oe["new_name"], oe["pol"])]]
    me = me[[p in allowed.get(f, set()) for f, p in zip(me["new_name"], me["pol"])]]
    me_by = {k: v for k, v in me.groupby("new_name")}
    oe_by = {k: v for k, v in oe.groupby("new_name")}

    def unmatched(f):
        m = me_by.get(f); o = oe_by.get(f)
        mt = m["peak"].values if m is not None else np.array([]); ot = o["peak_s"].values if o is not None else np.array([])
        res = [t for t in mt if len(ot) == 0 or np.min(np.abs(ot - t)) > TOL] + [t for t in ot if len(mt) == 0 or np.min(np.abs(mt - t)) > TOL]
        return np.array(res)

    # ---- which windows
    if wf is not None:
        info = rec.set_index("new_name")
        order = [(getattr(r, "stratum", ""), r.file, float(r.start_s), getattr(r, "why_chosen", "given")) for r in wf.itertuples()]
    else:
        info = rec.set_index("new_name")
        strata = rec.groupby("stratum").size()
        w = np.sqrt(strata); alloc = np.maximum(4, np.round(a.n_windows * w / w.sum())).astype(int)
        wins = {}
        for s, n in alloc.items():
            r_s = rec[rec["stratum"] == s]
            lst, pool = [], []
            for r in r_s.itertuples():
                for t in unmatched(r.new_name): pool.append((r.new_name, t, r.duration_s))
            rng.shuffle(pool)
            for f, t, dur in pool[: n // 2]:
                lst.append((f, float(np.clip(t - WIN * rng.uniform(0.2, 0.8), 0, max(0.1, dur - WIN))), "disagreement"))
            while len(lst) < n:
                r = r_s.iloc[rng.integers(len(r_s))]
                lst.append((r.new_name, float(rng.uniform(0, max(0.1, r.duration_s - WIN))), "random"))
            rng.shuffle(lst)
            wins[s] = lst
        order, keys = [], list(wins)
        while any(wins[k] for k in keys):
            for k in keys:
                if wins[k]: order.append((k,) + wins[k].pop())
        order = order[: a.n_windows + len(keys)]

    # ---- existing audit: continue numbering, skip duplicates
    old = _read_audit(xl)
    if old is not None:
        bk = out / "_backups"; bk.mkdir(exist_ok=True)
        shutil.copy2(xl, bk / f"audit_labels_{time.strftime('%Y%m%d_%H%M%S')}.xlsx")
        for b in sorted(bk.glob("audit_labels_*.xlsx"))[:-5]: b.unlink()
        have = set(zip(old["audit_windows"]["file"], old["audit_windows"]["start_s"].astype(float).round(3)))
        n_w = int(old["audit_windows"]["window"].str[1:].astype(int).max()); n_p = int(old["audit_windows"]["page"].astype(int).max())
    else:
        have, n_w, n_p = set(), 0, 0
    if a.first_window: n_w = a.first_window - 1
    if a.first_page: n_p = a.first_page - 1
    order = [o for o in order if (o[1], round(o[2], 3)) not in have]
    print(f"{len(order)} new windows to draw (existing audit: {n_w} windows, {n_p} pages)")
    if not order:
        return

    paths = {p.name: p for p in root.rglob("*.abf") if p.is_file()}
    rows, key, labels = [], [], []
    npages = math.ceil(len(order) / PER_PAGE)
    for pg in range(npages):
        chunk = order[pg * PER_PAGE:(pg + 1) * PER_PAGE]
        fig, axes = plt.subplots(len(chunk) * 2, 1, figsize=(15, 4.3 * len(chunk)), dpi=80)
        axes = np.atleast_1d(axes)
        for j, (stratum, f, st, why) in enumerate(chunk):
            wid = f"W{n_w + pg * PER_PAGE + j + 1:03d}"; page_no = n_p + pg + 1
            ab = pyabf.ABF(str(paths[f])); ab.setSweep(0)
            y = np.asarray(ab.sweepY, float); fs = float(ab.dataRate)
            sos = butter(4, 3000, "low", fs=fs, output="sos")
            i0, i1 = int(st * fs), int((st + WIN) * fs)
            lo, hi = max(0, i0 - int(2 * fs)), min(len(y), i1 + int(2 * fs))
            seg = y[lo:hi]; lp = sosfiltfilt(sos, seg)
            base = median_filter(lp[::20], size=max(3, int(1.0 * fs / 20)), mode="nearest")
            base = np.interp(np.arange(len(seg)), np.arange(0, len(seg), 20)[:len(base)], base)
            seg_d = seg - base; lp_d = lp - base
            sl = slice(i0 - lo, i1 - lo); tt = (np.arange(len(seg)) + lo) / fs
            ylo, yhi = np.percentile(lp_d[sl], 0.2), np.percentile(lp_d[sl], 99.8)
            pad = 0.25 * max(yhi - ylo, 10)
            m = me_by.get(f); o = oe_by.get(f)
            mw = m[(m["peak"] >= st) & (m["peak"] < st + WIN)].sort_values("peak") if m is not None else pd.DataFrame(columns=list(me.columns))
            ow = o[(o["peak_s"] >= st) & (o["peak_s"] < st + WIN)].sort_values("peak_s") if o is not None else pd.DataFrame(columns=list(oe.columns))
            ts = list(mw["peak"]) + list(ow["peak_s"])
            ev_amp = np.array([lp_d[int(t * fs) - lo] if 0 <= int(t * fs) - lo < len(lp_d) else 0 for t in ts]) if ts else np.array([0])
            ylo = min(ylo, ev_amp.min()) - pad; yhi = max(yhi, ev_amp.max()) + pad
            r = info.loc[f]
            for k, (lab, det, tcol) in enumerate([("miniML", mw, "peak"), ("old detector", ow, "peak_s")]):
                ax = axes[j * 2 + k]
                ax.plot(tt[sl], seg_d[sl], color="#c8c8c8", lw=0.5); ax.plot(tt[sl], lp_d[sl], color="k", lw=0.8)
                ax.set_xlim(st, st + WIN); ax.set_ylim(ylo, yhi)
                prefix = "M" if k == 0 else "O"
                for n_, (_, e) in enumerate(det.iterrows(), 1):
                    t = float(e[tcol]); idx = int(t * fs) - lo
                    yv = lp_d[idx] if 0 <= idx < len(lp_d) else 0
                    col = "#d6001c" if e["pol"] == "down" else "#1f5fbf"
                    ax.plot([t], [yv], marker="v" if e["pol"] == "down" else "^", color=col, ms=7, mec="k", mew=0.4, ls="none")
                    ax.annotate(f"{n_}", (t, 1.0), xycoords=("data", "axes fraction"), ha="center", va="bottom", fontsize=8, color=col, fontweight="bold")
                    other = oe_by.get(f) if k == 0 else me_by.get(f); ocol = "peak_s" if k == 0 else "peak"
                    matched = bool(other is not None and len(other) and np.min(np.abs(other[ocol].values - t)) <= TOL)
                    key.append(dict(window=wid, id=f"{prefix}{n_}", detector="miniML" if k == 0 else "old detector", time_s=t, polarity=e["pol"],
                                    amp_pA=float(e["amplitude"]) if k == 0 else float(e["amp_pA"]),
                                    score_or_snr=float(e["score"]) if k == 0 else float(e["snr"]), other_detector_agrees=matched))
                ax.set_ylabel(f"{lab}\n(pA)", fontsize=9)
                if k == 0:
                    ax.set_title(f"{wid}   {f}   t = {st:.1f}-{st + WIN:.1f} s   holding {r['holding_command_mV']:g} mV   {r['internal']}   [{why}]",
                                 fontsize=10, loc="left", pad=14)
                if k == 1: ax.set_xlabel("time in recording (s)", fontsize=9)
                else: ax.set_xticklabels([])
                for sp in ("top", "right"): ax.spines[sp].set_visible(False)
            rows.append(dict(window=wid, page=page_no, batch=a.batch, file=f, start_s=round(st, 3), end_s=round(st + WIN, 3),
                             holding_mV=r["holding_command_mV"], internal=r["internal"], drug=r["drug_tag"], stratum=stratum,
                             why_chosen=why, directions_shown=r["dirs"], n_miniML=len(mw), n_old=len(ow)))
            labels.append(dict(window=wid, page=page_no, batch=a.batch, holding_mV=r["holding_command_mV"], internal=r["internal"],
                               drug=r["drug_tag"], why_chosen=why, n_miniML=len(mw), n_old=len(ow), file=f, start_s=round(st, 3),
                               **{c: None for c in LABEL_COLS}))
        fig.tight_layout(); fig.savefig(out / f"page_{n_p + pg + 1:02d}.png"); plt.close(fig)
        print("page", n_p + pg + 1, flush=True)

    nw, nk, nl = pd.DataFrame(rows), pd.DataFrame(key), pd.DataFrame(labels)[INFO_COLS + LABEL_COLS]
    if old is not None:
        ow_, ok_, ol_ = old["audit_windows"], old["audit_events_key"], old["audit_labels"]
        for c in nw.columns:
            if c not in ow_.columns: ow_[c] = None
        for c in INFO_COLS:
            if c not in ol_.columns: ol_[c] = None
        nw = pd.concat([ow_, nw], ignore_index=True); nk = pd.concat([ok_, nk], ignore_index=True)
        nl = pd.concat([ol_[[c for c in INFO_COLS + LABEL_COLS if c in ol_.columns]], nl], ignore_index=True)
    _write_audit(xl, nl, nw, nk)
    print(f"audit workbook: {len(nl)} windows, {len(nk)} markers in the key, pages up to {n_p + npages}")


if __name__ == "__main__":
    main()
