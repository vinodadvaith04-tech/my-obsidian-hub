"""
build_audit_tool.py - build ONE self-contained HTML page (audit_tool.html) for clicking through the audit windows.

Reads   <audit>/audit_labels.xlsx   (sheets audit_windows + audit_events_key: the same windows and M/O markers as page_XX.png)
        the .abf recordings          (--data-root, searched recursively by file name)
Writes  <audit>/audit_tool.html      (double-click to open in Chrome/Safari/Firefox; works offline, nothing is uploaded)

The traces are drawn exactly like the PNG pages: 3 kHz low-pass, slow baseline removed. They are stored at half the
sampling rate (10 kHz) as 16-bit numbers inside the HTML, so the file is about 14 MB for 170 windows.

Re-run it whenever build_audit.py has appended new windows; your clicks are kept (they are stored per window id in the
browser and in the JSON you export). Only windows whose file and start time are unchanged keep their clicks.

Usage:  python scripts/build_audit_tool.py --data-root <folder with the .abf files> [--audit-dir <audit folder>] [--out audit_tool.html]
        add  --base-html audit/audit_tool.html  to REUSE the traces already stored in an earlier build (same window id, file,
        start time and number of markers) and read the .abf files only for the new windows. Handy when only some recordings
        are on this computer.
Needs: pyabf scipy numpy pandas openpyxl
"""
import argparse, base64, json, re, sys, time
from pathlib import Path
import numpy as np, pandas as pd, pyabf
from scipy.signal import butter, sosfiltfilt
from scipy.ndimage import median_filter

TOL = 0.010          # s: two markers closer than this are the same event (same value as other_detector_agrees in build_audit.py)
DEC = 2              # keep every 2nd sample (the trace is low-passed at 3 kHz)
VERSION = "1.0"


def window_trace(path, st, win):
    """Same processing as build_audit.py: returns (t0, dt, low-passed detrended trace in the window, helper for marker y)."""
    ab = pyabf.ABF(str(path)); ab.setSweep(0)
    y = np.asarray(ab.sweepY, float); fs = float(ab.dataRate)
    sos = butter(4, 3000, "low", fs=fs, output="sos")
    i0, i1 = int(st * fs), int((st + win) * fs)
    lo, hi = max(0, i0 - int(2 * fs)), min(len(y), i1 + int(2 * fs))
    seg = y[lo:hi]; lp = sosfiltfilt(sos, seg)
    base = median_filter(lp[::20], size=max(3, int(1.0 * fs / 20)), mode="nearest")
    base = np.interp(np.arange(len(seg)), np.arange(0, len(seg), 20)[:len(base)], base)
    lp_d = lp - base
    sl = slice(i0 - lo, min(i1, hi) - lo)
    return fs, i0, lo, lp_d, sl


def main():
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--audit-dir", type=Path, default=here.parent / "audit")
    ap.add_argument("--data-root", type=Path, required=True, help="folder that contains the .abf files (searched recursively)")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--template", type=Path, default=here / "audit_tool_template.html")
    ap.add_argument("--base-html", type=Path, default=None, help="earlier audit_tool.html whose windows are reused when unchanged")
    a = ap.parse_args()
    xl = a.audit_dir / "audit_labels.xlsx"
    win = pd.read_excel(xl, sheet_name="audit_windows")
    key = pd.read_excel(xl, sheet_name="audit_events_key")
    paths = {p.name: p for p in a.data_root.rglob("*.abf") if p.is_file()}
    base = {}
    if a.base_html:
        txt = a.base_html.read_text(encoding="utf-8")
        mm = re.search(r'<script id="data" type="application/json">(.*?)</script>', txt, re.S)
        base = {w["id"]: w for w in json.loads(mm.group(1))["windows"]}
    def reusable(r):
        w = base.get(r.window)
        n_key = int((key["window"] == r.window).sum())
        return (w is not None and w["file"] == r.file and abs(w["start"] - float(r.start_s)) < 1e-6
                and abs(w["end"] - float(r.end_s)) < 1e-6 and len(w["markers"]) == n_key
                and w["batch"] == str(r.batch) and w["page"] == int(r.page))
    miss = sorted({r.file for r in win.itertuples() if not reusable(r)} - set(paths))
    if miss:
        sys.exit(f"{len(miss)} recordings not found under {a.data_root}: {miss[:5]} ...")

    out_w = []
    n_reused = 0
    for r in win.itertuples():
        if reusable(r):
            out_w.append(base[r.window]); n_reused += 1; continue
        st, en = float(r.start_s), float(r.end_s); L = en - st
        fs, i0, lo, lp_d, sl = window_trace(paths[r.file], st, L)
        mk = key[key["window"] == r.window].copy()
        mk["det"] = mk["id"].str[0]
        mk = mk.sort_values(["det", "time_s"], ascending=[False, True])   # M first, then O; time order inside
        ev_amp, markers = [], []
        for m in mk.itertuples():
            idx = int(m.time_s * fs) - lo
            yv = float(lp_d[idx]) if 0 <= idx < len(lp_d) else 0.0
            ev_amp.append(yv)
            markers.append(dict(id=m.id, det=m.det, t=round(float(m.time_s), 5), pol=m.polarity, amp=round(float(m.amp_pA), 2),
                                score=round(float(m.score_or_snr), 3), yv=round(yv, 2)))
        seg = lp_d[sl]
        ylo, yhi = np.percentile(seg, 0.2), np.percentile(seg, 99.8)
        pad = 0.25 * max(yhi - ylo, 10)
        if ev_amp:
            ylo = min(ylo, min(ev_amp)); yhi = max(yhi, max(ev_amp))
        ylo -= pad; yhi += pad
        arr = seg[::DEC]
        scale = float(max(np.abs(arr).max(), 1e-6) / 32000.0)
        q = np.round(arr / scale).astype("<i2")
        out_w.append(dict(
            id=r.window, page=int(r.page), batch=str(r.batch), file=r.file, start=st, end=en,
            holding=float(r.holding_mV), internal=str(r.internal), drug=str(r.drug), why=str(r.why_chosen),
            dirs=str(r.directions_shown), t0=i0 / fs, dt=DEC / fs, scale=scale,
            ylo=round(float(ylo), 2), yhi=round(float(yhi), 2), n=int(len(q)),
            y=base64.b64encode(q.tobytes()).decode("ascii"), markers=markers))
        print(r.window, len(markers), "markers", flush=True)

    data = dict(version=VERSION, tol=TOL, built=time.strftime("%Y-%m-%d %H:%M"), windows=out_w)
    html = a.template.read_text(encoding="utf-8").replace("/*__DATA__*/", json.dumps(data, separators=(",", ":")))
    out = a.out or (a.audit_dir / "audit_tool.html")
    out.write_text(html, encoding="utf-8")
    print(f"wrote {out}  ({out.stat().st_size / 1e6:.1f} MB, {len(out_w)} windows, {n_reused} reused from the earlier build)")


if __name__ == "__main__":
    main()
