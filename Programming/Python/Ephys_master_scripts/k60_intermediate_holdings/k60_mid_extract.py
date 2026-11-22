#!/usr/bin/env python3
"""
k60_0_extract.py  -  stability / noise / drift / transient features for K-gluconate,
no-blocker, no-TTX gap-free recordings at -60 mV and ~0 mV (command holdings).

Run from the Ephys folder:  python3 _work/k60_0/k60_0_extract.py
Reads  Ephys_master.xlsx (read-only) and the raw .abf files.
Writes _work/k60_0/k60_gapfree_features.csv, k60_capcheck.csv, k60_transients.csv
Nothing in the master is changed.
"""
import os, sys, numpy as np, pandas as pd, pyabf
from scipy.signal import butter, sosfiltfilt, find_peaks, peak_widths
from scipy.optimize import curve_fit

ROOT = os.getcwd()
SET = os.environ.get("K60_SET", "kglu")          # kglu (default) | cs60 | mid
OUT = os.path.join(ROOT, "_work", "k60_mid" if SET == "mid" else "k60_0")
os.makedirs(OUT, exist_ok=True)
HOLDS = [-60, 0, 5, 10, 15]

files = pd.read_excel(os.path.join(ROOT, "Ephys_master.xlsx"), sheet_name="files")
SUF = "" if SET == "kglu" else "_" + SET
if SET in ("kglu", "mid"):
    k = files[(files.internal == "KGlu") & (files.drug.isna()) & (files.session != "Sam EPSC (22 Oct)")].copy()
else:
    k = files[(files.internal == "Cs-gluconate") & (files.drug.isna())].copy()
    HOLDS = [-60]
if SET == "mid":
    HOLDS = [-50, -45, -40, -35, -30]            # intermediate command holdings
k["t_s"] = pd.to_datetime(k.date_header.astype(int).astype(str) + " " + k.start_time.astype(str),
                          format="%Y%m%d %H:%M:%S")
g = k[(k.protocol_type == "sPSC gap-free") & k.holding_command_mV.isin(HOLDS)].copy()
cells = g[["date_header", "slice", "cell"]].drop_duplicates()
cc = pd.concat([k[(k.protocol_type == "cap check") & (k.date_header == c.date_header) &
                  (k.slice == c.slice) & (k.cell == c.cell)] for c in cells.itertuples()])

# ------------------------------------------------------------------ seal test (cap check)
def cap_fit(path):
    a = pyabf.ABF(path)
    fs = a.sampleRate
    a.setSweep(0)
    c = a.sweepC
    ch = np.flatnonzero(np.diff(c) != 0)
    vhold_cap = np.nan
    if len(ch) >= 2:                               # step stored in the DAC waveform
        i0, i1 = ch[0] + 1, ch[1] + 1
        dv = float(c[i0] - c[0])
    else:                                          # MTest files: read the step from the voltage monitor (channel 1)
        if a.channelCount < 2:
            return None
        a.setSweep(0, channel=1)
        v = a.sweepY.copy()
        if np.abs(np.median(v)) > 200:             # voltage monitor stored as 'pA' = mV x 20 (e.g. 6 Nov files)
            v = v / 20.0
        e = np.flatnonzero(np.abs(np.diff(v)) > 0.5 * np.ptp(v))
        if len(e) < 2 or np.ptp(v) < 1:
            return None
        i0, i1 = e[0] + 1, e[1] + 1
        dv = float(np.median(v[i0 + 10:i1 - 10]) - np.median(v[:i0 - 5]))
        vhold_cap = float(np.median(v[:i0 - 5]))
        a.setSweep(0, channel=0)
    ys = np.array([(a.setSweep(s, channel=0), a.sweepY.copy())[1] for s in range(a.sweepCount)])
    res = []
    for y in ys:
        ib = y[5:i0 - 5].mean()
        iss = y[i1 - 40:i1 - 5].mean()
        seg = y[i0 + 1:i0 + 1 + 200] - iss      # decay towards steady state
        tt = np.arange(len(seg)) / fs * 1000.0  # ms
        try:
            p, _ = curve_fit(lambda t, A, tau: A * np.exp(-t / tau), tt, seg, p0=[seg[0], 2.0],
                             bounds=([-1e5, 0.05], [1e5, 50]), maxfev=2000)
        except Exception:
            continue
        A, tau = p
        d_iss = iss - ib
        d_0 = d_iss + A                          # instantaneous current at step onset
        if d_iss <= 0 or d_0 <= 0:
            continue
        ra = dv / d_0 * 1000.0                   # MOhm
        rin = dv / d_iss * 1000.0
        rm = rin - ra
        if rm <= 0:
            continue
        cm = tau / ra * 1000.0 * (rin / rm)      # pF  (Cm = tau*(Ra+Rm)/(Ra*Rm))
        res.append((ra, rm, cm, ib, tau))
    if not res:
        return None
    r = np.array(res)
    return dict(Ra_MOhm=np.median(r[:, 0]), Rm_MOhm=np.median(r[:, 1]), Cm_pF=np.median(r[:, 2]),
                Ibase_pA=np.median(r[:, 3]), tau_ms=np.median(r[:, 4]), n_sweeps_ok=len(r), step_mV=dv,
                Vhold_capcheck_mV=vhold_cap)

rows = []
for r in cc.itertuples():
    try:
        d = cap_fit(r.path_on_mac)
    except Exception as e:
        d = None
        print("capcheck fail", r.new_name, e, file=sys.stderr)
    row = dict(new_name=r.new_name, date=int(r.date_header), slice=int(r.slice), cell=int(r.cell),
               time=r.t_s.isoformat())
    if d:
        row.update(d)
    rows.append(row)
pd.DataFrame(rows).to_csv(os.path.join(OUT, "k60_capcheck"+SUF+".csv"), index=False)
print("cap checks:", len(rows))

# ------------------------------------------------------------------ gap-free features
sos3k = None
feat, trans = [], []
for n, r in enumerate(g.itertuples()):
    a = pyabf.ABF(r.path_on_mac)
    fs = a.sampleRate
    a.setSweep(0)
    y = a.sweepY.astype(np.float64)
    N = len(y)
    dur = N / fs
    # --- holding current and drift (1-s block medians)
    nb = int(dur)
    bl = np.median(y[:nb * fs].reshape(nb, fs), axis=1)
    tb = np.arange(nb) + 0.5
    slope, icpt = np.polyfit(tb, bl, 1)           # pA/s
    fitres = bl - (slope * tb + icpt)
    first10, last10 = np.median(y[:10 * fs]), np.median(y[-10 * fs:])
    sm5 = np.convolve(bl, np.ones(5) / 5, mode="valid")
    # --- noise
    d1 = np.diff(y)
    noise_hf = 1.4826 * np.median(np.abs(d1 - np.median(d1))) / np.sqrt(2)
    if sos3k is None or sos3k[0] != fs:
        sos3k = (fs, butter(4, 3000, btype="low", fs=fs, output="sos"))
    ylp = sosfiltfilt(sos3k[1], y)
    blk = 2000 if fs == 20000 else int(0.1 * fs)
    nbk = N // blk
    yb = ylp[:nbk * blk].reshape(nbk, blk)
    res = (yb - np.median(yb, axis=1, keepdims=True)).ravel()
    noise_lp = 1.4826 * np.median(np.abs(res - np.median(res)))
    # raw (unfiltered) detrended SD, robust
    yr = y[:nbk * blk].reshape(nbk, blk)
    resr = (yr - np.median(yr, axis=1, keepdims=True)).ravel()
    noise_raw = 1.4826 * np.median(np.abs(resr - np.median(resr)))
    # --- fast large transients (candidate unclamped spikes / large events)
    thr = 40.0
    pk = []
    for sign in (+1, -1):
        sig = sign * res
        idx, pr = find_peaks(sig, height=thr, prominence=thr, distance=int(0.003 * fs))
        if len(idx):
            w = peak_widths(sig, idx, rel_height=0.5)[0] / fs * 1000.0
            for i, amp, pm, ww in zip(idx, sig[idx], pr["prominences"], w):
                # biphasic check: opposite-sign excursion within +-3 ms
                lo, hi = max(0, i - int(0.003 * fs)), min(len(res), i + int(0.003 * fs))
                opp = (-sign * res[lo:hi]).max()
                pk.append((n, r.new_name, sign, i / fs, amp, pm, ww, opp))
    trans.extend(pk)
    feat.append(dict(
        new_name=r.new_name, session=r.session, date=int(r.date_header), slice=int(r.slice), cell=int(r.cell),
        holding_command_mV=r.holding_command_mV, duration_s=dur, start_time=r.t_s.isoformat(),
        Ihold_median_pA=float(np.median(y)), Ihold_first10s_pA=float(first10), Ihold_last10s_pA=float(last10),
        drift_slope_pA_per_min=slope * 60, drift_total_pA=slope * dur, drift_end_minus_start_pA=float(last10 - first10),
        wander_sd_pA=float(np.std(fitres)), range_5s_median_pA=float(sm5.max() - sm5.min()),
        noise_hf_pA=float(noise_hf), noise_lp3k_pA=float(noise_lp), noise_raw_pA=float(noise_raw),
        y_min=float(y.min()), y_max=float(y.max())))
    # keep 1-s baseline series for plots
    np.save(os.path.join(OUT, "bl" + SUF + "_" + r.new_name.replace(".abf", ".npy")), bl.astype(np.float32))
    print(n + 1, len(g), r.new_name, flush=True)

pd.DataFrame(feat).to_csv(os.path.join(OUT, "k60_gapfree_features"+SUF+".csv"), index=False)
pd.DataFrame(trans, columns=["i", "new_name", "sign", "t_s", "amp_pA", "prom_pA", "fwhm_ms", "opp_pA"]).drop(columns="i") \
    .to_csv(os.path.join(OUT, "k60_transients"+SUF+".csv"), index=False)
print("done")
