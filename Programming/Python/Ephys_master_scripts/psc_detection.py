"""
psc_detection.py  -  spontaneous PSC (sEPSC / sIPSC) detector for gap-free voltage-clamp recordings
=====================================================================================================
Stand-alone event detection (nothing here depends on ephys_analysis.py). Two stages:

  STAGE 1  candidates      Clements & Bekkers (1997) matched filter, run separately for inward (EPSC-like) and
                           outward (IPSC-like) events with a fast and a slow template each, deliberately permissive
                           (criterion >= CAND_CRIT). Baseline = rolling median with the events masked out (2 passes).
  STAGE 2  verification    every candidate is fitted with a rise-and-decay function
                               i(t) = c + A * [exp(-(t-t0)/tau_decay) - exp(-(t-t0)/tau_rise)] / peak   (tau_rise < tau_decay)
                           and accepted only if it passes amplitude, kinetics and fit-quality gates. All reported
                           amplitudes / rise / decay / half-width / charge come from the FITTED curve, not the noisy samples.

Amplitude limit: an event is accepted when its fitted amplitude >= max(ABS_MIN_PA, K_SD[pol] * noise_lp) of that recording
(its own detection limit, written to the output). For comparisons between cells use the columns that count events above one
COMMON cutoff (COMMON_CUTOFF_PA) and the per-recording injected-event detection probabilities (psc_validate.py).

Usage
  python psc_detection.py --root "C:/Users/admin/Documents/1. Ephys" --out "C:/Users/admin/Documents/1. Ephys/_psc_detection"
  python psc_detection.py --file some_GapFree_0003.abf --plot      (one file, quick look png next to the output)
Outputs (in --out):  psc_events.csv, psc_candidates_rejected.csv, psc_recordings.csv, psc_event_averages.npz, psc_detection_log.txt
"""
import argparse
import csv
import json
import math
import re
import sys
import time
import warnings
from collections import OrderedDict
from pathlib import Path

import numpy as np
from scipy.ndimage import median_filter
from scipy.optimize import least_squares
from scipy.signal import butter, fftconvolve, find_peaks, sosfiltfilt

warnings.filterwarnings("ignore")

CFG = dict(
    LP_HZ=3000.0,              # zero-phase low-pass used for fitting / measuring (detection also runs on this trace)
    BASE_BLOCK_S=0.005,        # baseline: median of 5 ms blocks (event samples excluded) ...
    BASE_MEDIAN_S=0.5,         # ... median-filtered over 0.5 s
    CAND_CRIT=3.0,             # stage-1 matched-filter criterion (scale / SE). Permissive on purpose.
    PREGATE=0.6,               # a candidate is only fitted if its rough amplitude >= PREGATE x its recording cutoff
    # template (rise_ms, decay_ms, window_ms); two per polarity (fast + slow)
    TEMPLATES=dict(inward=[(0.6, 3.0, 25.0), (1.5, 7.0, 45.0)], outward=[(1.5, 8.0, 60.0), (4.0, 22.0, 120.0)]),
    PEAK_WIN_S=dict(inward=0.012, outward=0.030),    # peak is searched this long after the matched-filter window start
    PRE_S=dict(inward=0.006, outward=0.025),         # fit window starts this long before the rough peak
    DECAY_GUESS_MS=dict(inward=5.0, outward=14.0),
    RISE_GUESS_MS=dict(inward=0.7, outward=2.5),
    # --- amplitude limit: max(ABS_MIN_PA, K_SD x noise_lp of the recording)
    ABS_MIN_PA=4.0,
    K_SD=dict(inward=4.0, outward=4.5),
    MAX_AMP_PA=250.0,          # bigger = unclamped spike / artefact: excluded, 0.3 s afterwards masked
    MASK_AFTER_S=0.3,
    UNSTABLE_FACTOR=2.0,       # 0.5 s blocks whose noise (MAD) exceeds this x the median block noise (and +1.5 pA) are masked
    # --- kinetic gates on the fitted curve (ms)
    GATE=dict(
        inward=dict(decay=(0.8, 30.0), rise1090=(0.08, 6.0)),
        outward=dict(decay=(2.0, 120.0), rise1090=(0.3, 25.0))),
    MAX_RESID_RATIO=2.2,       # RMS of fit residual / noise_lp over the fit window
    AMP_RATIO=(0.55, 1.8),     # fitted amplitude / rough amplitude
    MIN_PEAK_SEP_S=0.002,      # fitted peaks closer than this = duplicate detection
    TWIN_S=0.010,              # opposite-polarity candidates this close: only the larger is kept
    OVERLAP_FACTOR=3.0,        # an event is 'overlapped' if the next event starts within this x its decay tau
    ISOLATION_S=dict(inward=0.025, outward=0.060),   # isolated event: no other event this long before / after the peak
    COMMON_CUTOFF_PA=10.0,     # columns 'n_ge_common' / 'freq_ge_common_Hz' count only events >= this
    MIN_EVENTS_KINETICS=10,
)

# ------------------------------------------------------------------------------------------------
# signal helpers
# ------------------------------------------------------------------------------------------------


def load_abf(path):
    import pyabf
    abf = pyabf.ABF(str(path))
    abf.setSweep(0, channel=0)
    return np.asarray(abf.sweepY, float), float(abf.dataRate), str(abf.adcUnits[0])


def _template(fs, rise_ms, decay_ms, win_ms):
    t = np.arange(int(win_ms * 1e-3 * fs)) / fs * 1e3
    T = (1 - np.exp(-t / rise_ms)) * np.exp(-t / decay_ms)
    return T / T.max()


def rolling_baseline(y, fs, cfg, event_mask=None):
    """slow baseline: median of 5 ms blocks (masked samples ignored), median-filtered over BASE_MEDIAN_S, interpolated"""
    q = max(1, int(cfg["BASE_BLOCK_S"] * fs))
    nb = len(y) // q
    yy = y[: nb * q].copy()
    if event_mask is not None:
        yy[event_mask[: nb * q]] = np.nan
    blk = np.nanmedian(yy.reshape(nb, q), axis=1)
    good = np.isfinite(blk)
    xc = (np.arange(nb) + 0.5) * q
    if good.sum() < 5:
        return np.full(len(y), np.nanmedian(y))
    blk = np.interp(xc, xc[good], blk[good])
    k = max(3, int(cfg["BASE_MEDIAN_S"] / cfg["BASE_BLOCK_S"]) | 1)
    bm = median_filter(blk, size=k, mode="nearest")
    return np.interp(np.arange(len(y)), xc, bm)


def matched_filter_criterion(y, T):
    """Clements & Bekkers detection criterion for every window start (events must be positive-going in y)"""
    L = len(T)
    n = len(y)
    sT = T.sum()
    sTT = (T ** 2).sum()
    cs = np.concatenate([[0.0], np.cumsum(y)])
    cs2 = np.concatenate([[0.0], np.cumsum(y ** 2)])
    m = n - L + 1
    sy = cs[L:L + m] - cs[:m]
    syy = cs2[L:L + m] - cs2[:m]
    syT = fftconvolve(y, T[::-1], mode="valid")[:m]
    denom = sTT - sT ** 2 / L
    scale = (syT - sT * sy / L) / denom
    c = (sy - scale * sT) / L
    sse = syy + scale ** 2 * sTT + L * c ** 2 - 2 * (scale * syT + c * sy - scale * c * sT)
    se = np.sqrt(np.maximum(sse, 1e-12) / (L - 1))
    crit = scale / se
    return np.concatenate([crit, np.zeros(n - m)]), np.concatenate([scale, np.zeros(n - m)])


def mad_sd(x):
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def unstable_mask(ylp, fs, cfg):
    """mask 0.5 s blocks whose noise is far above the file's median (seal / membrane movement)"""
    n = len(ylp)
    mask = np.zeros(n, bool)
    blk = int(0.5 * fs)
    nb = n // blk
    if nb >= 4:
        seg = ylp[: nb * blk].reshape(nb, blk)
        mad_b = 1.4826 * np.median(np.abs(seg - np.median(seg, axis=1, keepdims=True)), axis=1)
        ref = np.median(mad_b)
        bad = (mad_b > cfg["UNSTABLE_FACTOR"] * ref) & (mad_b > ref + 1.5)
        bad2 = bad.copy()
        bad2[1:] |= bad[:-1]
        bad2[:-1] |= bad[1:]
        for i in np.where(bad2)[0]:
            mask[i * blk:(i + 1) * blk] = True
    return mask


# ------------------------------------------------------------------------------------------------
# the event model
# ------------------------------------------------------------------------------------------------


def _shape(t, t0, tr, td):
    """unnormalised: exp(-s/td) - exp(-s/tr), s = t - t0 (zero before onset); tr = rise time constant < td = decay time constant"""
    s = np.maximum(t - t0, 0.0)
    return np.exp(-s / td) - np.exp(-s / tr)


def _shape_peak(tr, td):
    a = 1.0 / td
    b = 1.0 / tr
    ts = math.log(b / a) / (b - a)
    return ts, math.exp(-a * ts) - math.exp(-b * ts)


def model(t, c, A, t0, tr, td):
    ts, pk = _shape_peak(tr, td)
    return c + A * _shape(t, t0, tr, td) / pk


def curve_kinetics(tr, td):
    """peak latency, 10-90 rise, half-width and area (per unit amplitude, ms and ms) of the model curve"""
    ts, pk = _shape_peak(tr, td)
    dt = min(tr, td) / 40.0
    t = np.arange(0, td * 8 + ts + 5 * tr, dt)
    f = _shape(t, 0.0, tr, td) / pk
    ip = int(np.argmax(f))
    up = f[: ip + 1]

    def cross_up(level):
        j = int(np.searchsorted(up, level))
        if j <= 0:
            return 0.0
        return t[j - 1] + (level - up[j - 1]) / (up[j] - up[j - 1]) * dt

    rise = cross_up(0.9) - cross_up(0.1)
    above = np.where(f >= 0.5)[0]
    hw = (above[-1] - above[0]) * dt
    area = (td - tr) / pk          # ms (x amplitude pA = fC)
    return ts, rise, hw, area


# ------------------------------------------------------------------------------------------------
# stage 1: candidates
# ------------------------------------------------------------------------------------------------


def candidates(y, ylp, fs, cfg, pol, thr=None, templates=None, valid=None):
    """matched-filter candidates for one polarity. Returns list of dict(i0, crit, tmpl). y: baseline-subtracted raw trace."""
    sgn = -1.0 if pol == "inward" else 1.0
    thr = cfg["CAND_CRIT"] if thr is None else thr
    out = {}
    min_sep = max(1, int(0.002 * fs))
    for k, (rise, dec, win) in enumerate(templates or cfg["TEMPLATES"][pol]):
        crit, scale = matched_filter_criterion(sgn * ylp, _template(fs, rise, dec, win))
        crit = np.where(scale > 0, crit, 0.0)
        pk, _ = find_peaks(crit, height=thr, distance=min_sep)
        for p in pk:
            p = int(p)
            if valid is not None and not valid[p]:
                continue
            near = [q for q in out if abs(q - p) < int(0.004 * fs)]
            if near:
                q = near[0]
                if crit[p] > out[q]["crit"]:
                    del out[q]
                    out[p] = dict(i0=p, crit=float(crit[p]), tmpl=k)
            else:
                out[p] = dict(i0=p, crit=float(crit[p]), tmpl=k)
    return [out[k] for k in sorted(out)]


# ------------------------------------------------------------------------------------------------
# stage 2: fit + gates
# ------------------------------------------------------------------------------------------------


def rough_peak(yp, fs, i0, pol, cfg):
    """highest point of the positive-going trace yp in the PEAK_WIN after i0 (start 2 ms earlier); None if on a slope"""
    pw = int(cfg["PEAK_WIN_S"][pol] * fs)
    s0 = max(0, i0 - int(0.002 * fs))
    seg = yp[s0:i0 + pw]
    if len(seg) < 8:
        return None
    ip = s0 + int(np.argmax(seg))
    if ip - s0 < 2 or (s0 + len(seg) - 1) - ip < 2:
        return None
    b0, b1 = max(0, s0 - int(0.004 * fs)), max(1, s0)
    base = float(np.median(yp[b0:b1]))
    return ip, float(yp[ip] - base), base


def fit_event(yp, fs, ip, amp0, base0, noise, pol, cfg, next_ip=None):
    """fit one event in the positive-going low-passed trace yp. Returns dict or None."""
    pre = int(cfg["PRE_S"][pol] * fs)
    td0 = cfg["DECAY_GUESS_MS"][pol]
    tr0 = cfg["RISE_GUESS_MS"][pol]
    post = int(5 * td0 * 1e-3 * fs)
    i_a = max(0, ip - pre)
    i_b = min(len(yp), ip + post)
    truncated = False
    if next_ip is not None and next_ip - int(0.001 * fs) < i_b:
        i_b = max(ip + int(0.003 * fs), next_ip - int(0.001 * fs))
        truncated = True
    if i_b - i_a < int(0.012 * fs) or ip - i_a < 5:
        return None
    t = (np.arange(i_a, i_b) - ip) / fs * 1e3             # ms relative to rough peak
    y = yp[i_a:i_b]
    dec_hi = 60.0 if pol == "inward" else 200.0
    rise_hi = 15.0 if pol == "inward" else 40.0
    t0_lo = -pre / fs * 1e3 + 1.0
    lo = [-max(3 * noise, 2.0), 0.3 * amp0, t0_lo, 0.05, 0.3]
    hi = [max(3 * noise, 2.0), 3.0 * amp0, -0.05, rise_hi, dec_hi]
    p0 = [base0 * 0 + 0.0, amp0, -(2.0 * tr0 + 0.3), tr0, td0]
    p0[2] = min(max(p0[2], t0_lo + 0.1), -0.1)

    # parameters: c, A, t0, tr (rise tau), k = td / tr  (k >= 1.05 so the rise is always the faster process)
    lo[3:5] = [0.04, 1.05]
    hi[3:5] = [rise_hi, 400.0]
    p0[3:5] = [tr0, td0 / tr0]

    def resid(p):
        return model(t, p[0], p[1], p[2], p[3], p[3] * p[4]) - y

    try:
        r = least_squares(resid, p0, bounds=(lo, hi), method="trf", loss="linear", max_nfev=150,
                          x_scale=[noise, max(amp0, noise), 1.0, 1.0, 3.0])
    except Exception:
        return None
    c, A, t0, tr, k = r.x
    td = tr * k
    ts, rise, hw, area = curve_kinetics(tr, td)
    res = r.fun
    return dict(c=float(c), A=float(A), t0_ms=float(t0), tr=float(tr), td=float(td), peak_lat_ms=float(ts), rise1090_ms=float(rise),
                halfwidth_ms=float(hw), charge_fC=float(A * area), resid_rms=float(np.sqrt(np.mean(res ** 2))),
                resid_ratio=float(np.sqrt(np.mean(res ** 2)) / noise), truncated=truncated,
                i_peak=int(round(ip + (t0 + ts) * 1e-3 * fs)), i_onset=int(round(ip + t0 * 1e-3 * fs)),
                at_bound=bool(td > 0.97 * dec_hi or tr > 0.97 * rise_hi or tr < 0.045 or A > 2.95 * amp0))


def amp_limit(noise_lp, pol, cfg, min_amp=None):
    return max(cfg["ABS_MIN_PA"] if min_amp is None else min_amp, cfg["K_SD"][pol] * noise_lp)


# ------------------------------------------------------------------------------------------------
# full recording
# ------------------------------------------------------------------------------------------------


def detect(I, fs, cfg=CFG, min_amp=None, return_aux=False):
    """detect events in one gap-free current trace (pA). Returns (events, rejected, info[, aux])."""
    I = np.asarray(I, float)
    n = len(I)
    sos = butter(4, cfg["LP_HZ"], btype="low", fs=fs, output="sos")

    # ---- pass 1: coarse baseline, strong candidates only -> event mask -> pass 2 baseline without events
    base1 = rolling_baseline(I, fs, cfg)
    y1 = I - base1
    y1lp = sosfiltfilt(sos, y1)
    emask = np.zeros(n, bool)
    for pol in ("inward", "outward"):
        for c in candidates(y1, y1lp, fs, cfg, pol, thr=4.5, templates=cfg["TEMPLATES"][pol][:1]):
            dec = cfg["DECAY_GUESS_MS"][pol]
            emask[max(0, c["i0"] - int(0.003 * fs)): min(n, c["i0"] + int(4 * dec * 1e-3 * fs) + int(0.01 * fs))] = True
    base = rolling_baseline(I, fs, cfg, event_mask=emask)
    y = I - base
    ylp = sosfiltfilt(sos, y)
    mask = unstable_mask(ylp, fs, cfg)
    quiet = ~(mask | emask)
    if quiet.sum() < 5 * fs:
        quiet = ~mask if (~mask).sum() > 5 * fs else np.ones(n, bool)
    noise_lp = mad_sd(ylp[quiet])
    noise_raw = mad_sd(y[quiet])

    events, rejected = [], []
    per_pol = {}
    for pol, sgn in (("inward", -1.0), ("outward", 1.0)):
        yp = sgn * ylp
        limit = amp_limit(noise_lp, pol, cfg, min_amp)
        cands = candidates(y, ylp, fs, cfg, pol, valid=~mask)
        # rough peak + pre-gate
        pre = []
        for c in cands:
            rp = rough_peak(yp, fs, c["i0"], pol, cfg)
            if rp is None:
                continue
            ip, amp, base_local = rp
            if amp < cfg["PREGATE"] * limit:
                continue
            c.update(ip=ip, amp0=amp, base0=base_local)
            pre.append(c)
        # same peak found by two templates -> keep larger criterion
        pre.sort(key=lambda c: c["ip"])
        dd = []
        for c in pre:
            if dd and c["ip"] - dd[-1]["ip"] < int(cfg["MIN_PEAK_SEP_S"] * fs):
                if c["crit"] > dd[-1]["crit"]:
                    dd[-1] = c
            else:
                dd.append(c)
        per_pol[pol] = dict(n_cand=len(cands), n_fitted=len(dd), limit=limit)
        for j, c in enumerate(dd):
            nxt = dd[j + 1]["ip"] if j + 1 < len(dd) else None
            f = fit_event(yp, fs, c["ip"], c["amp0"], c["base0"], noise_lp, pol, cfg, next_ip=nxt)
            rec = dict(polarity=pol, crit=c["crit"], template=c["tmpl"], rough_amp_pA=c["amp0"], limit_pA=limit,
                       next_peak_gap_ms=((nxt - c["ip"]) / fs * 1e3 if nxt is not None else np.nan))
            if f is None:
                rec.update(peak_s=c["ip"] / fs, onset_s=c["i0"] / fs, reason="fit failed / window too short")
                rejected.append(rec)
                continue
            g = cfg["GATE"][pol]
            reason = None
            if f["A"] < limit:
                reason = "fitted amplitude below limit"
            elif f["A"] > cfg["MAX_AMP_PA"]:
                reason = "excluded: > MAX_AMP_PA (artefact)"
            elif not (cfg["AMP_RATIO"][0] <= f["A"] / c["amp0"] <= cfg["AMP_RATIO"][1]):
                reason = "fit amplitude disagrees with data peak"
            elif f["resid_ratio"] > cfg["MAX_RESID_RATIO"]:
                reason = "poor fit (residual too large)"
            elif not (g["rise1090"][0] <= f["rise1090_ms"] <= g["rise1090"][1]):
                reason = "rise time outside physiological range"
            elif not f["truncated"] and not (g["decay"][0] <= f["td"] <= g["decay"][1]):
                reason = "decay outside physiological range"
            elif f["at_bound"] and not f["truncated"]:
                reason = "fit parameter at bound"
            rec.update(f)
            rec.update(peak_s=f["i_peak"] / fs, onset_s=f["i_onset"] / fs, amp_pA=f["A"])
            if reason:
                rec["reason"] = reason
                rejected.append(rec)
            else:
                rec["status"] = "ok"
                events.append(rec)
    # ---- opposite-polarity twins (one biphasic transient seen by both): keep the larger
    ev = sorted(events, key=lambda e: e["peak_s"])
    twin = int(cfg["TWIN_S"] * fs) / fs
    for a, b in zip(ev[:-1], ev[1:]):
        if a["status"] == "ok" and b["status"] == "ok" and a["polarity"] != b["polarity"] and abs(b["peak_s"] - a["peak_s"]) < twin:
            loser = a if a["amp_pA"] < b["amp_pA"] else b
            loser["status"] = "rejected"
            loser["reason"] = "opposite-polarity twin"
    # ---- same polarity duplicates after fitting (two candidates fitting the same event)
    for pol in ("inward", "outward"):
        last = None
        for e in [e for e in ev if e["polarity"] == pol and e["status"] == "ok"]:
            if last is not None and e["peak_s"] - last["peak_s"] < cfg["MIN_PEAK_SEP_S"] * 2:
                loser = e if e["amp_pA"] <= last["amp_pA"] else last
                loser["status"] = "rejected"
                loser["reason"] = "duplicate of neighbouring detection"
                last = e if loser is last else last
            else:
                last = e
    # ---- large-artefact masking, unstable masking
    for e in ev:
        if e["status"] != "ok":
            continue
        i = int(e["peak_s"] * fs)
        if mask[min(n - 1, i)]:
            e["status"] = "rejected"
            e["reason"] = "inside unstable stretch"
    for r in rejected:
        if r.get("reason", "").startswith("excluded: >") and r.get("peak_s") is not None:
            i = int(r["peak_s"] * fs)
            mask[max(0, i - int(0.02 * fs)): min(n, i + int(cfg["MASK_AFTER_S"] * fs))] = True
    for e in ev:
        if e["status"] == "ok" and mask[min(n - 1, int(e["peak_s"] * fs))]:
            e["status"] = "rejected"
            e["reason"] = "inside masked stretch (after artefact)"
    rejected += [e for e in ev if e["status"] == "rejected"]
    ok = sorted([e for e in ev if e["status"] == "ok"], key=lambda e: e["peak_s"])
    # ---- overlap / isolation flags and inter-event intervals (all polarities together)
    for k, e in enumerate(ok):
        nxt = ok[k + 1] if k + 1 < len(ok) else None
        prv = ok[k - 1] if k > 0 else None
        e["iei_prev_ms"] = (e["peak_s"] - prv["peak_s"]) * 1e3 if prv is not None else np.nan
        gap_next = (nxt["onset_s"] - e["onset_s"]) * 1e3 if nxt is not None else np.inf
        e["overlapped"] = bool(gap_next < cfg["OVERLAP_FACTOR"] * e["td"] or e["truncated"])
        iso = cfg["ISOLATION_S"][e["polarity"]]
        e["isolated"] = bool(not e["overlapped"] and (prv is None or e["onset_s"] - prv["peak_s"] > iso) and (nxt is None or nxt["onset_s"] - e["peak_s"] > iso))
        e["noise_lp_pA"] = noise_lp
        e["snr"] = e["amp_pA"] / noise_lp
        e["baseline_pA"] = float(base[min(n - 1, int(e["peak_s"] * fs))])
    info = dict(noise_lp_pA=noise_lp, noise_raw_pA=noise_raw, masked_s=float(mask.sum() / fs), duration_s=n / fs,
                limit_inward_pA=per_pol["inward"]["limit"], limit_outward_pA=per_pol["outward"]["limit"],
                n_cand_inward=per_pol["inward"]["n_cand"], n_cand_outward=per_pol["outward"]["n_cand"],
                n_fit_inward=per_pol["inward"]["n_fitted"], n_fit_outward=per_pol["outward"]["n_fitted"])
    if return_aux:
        return ok, rejected, info, dict(base=base, ylp=ylp, mask=mask, y=y)
    return ok, rejected, info


# ------------------------------------------------------------------------------------------------
# summaries
# ------------------------------------------------------------------------------------------------


def average_event(ylp, fs, evs, pol, pre_ms=5.0, post_ms=None):
    """average waveform of isolated events aligned on the fitted onset (baseline = 3 ms before onset) + biexponential fit"""
    post_ms = post_ms or (45.0 if pol == "inward" else 120.0)
    sgn = -1.0 if pol == "inward" else 1.0
    a, b = int(pre_ms * 1e-3 * fs), int(post_ms * 1e-3 * fs)
    W = []
    for e in evs:
        i = int(round(e["onset_s"] * fs))
        if i - a < 0 or i + b > len(ylp):
            continue
        w = sgn * ylp[i - a:i + b]
        w = w - np.median(w[max(0, a - int(0.003 * fs)):a])
        W.append(w)
    if len(W) < 5:
        return None
    W = np.array(W)
    avg = W.mean(axis=0)
    t = (np.arange(len(avg)) - a) / fs * 1e3
    out = dict(n=len(W), t_ms=t, avg=avg, sem=W.std(axis=0, ddof=1) / np.sqrt(len(W)))
    try:
        A0 = float(avg.max())
        lo = [-0.2 * A0, 0.3 * A0, -2.0, 0.04, 1.05]
        hi = [0.2 * A0, 3 * A0, 2.0, 15.0 if pol == "inward" else 40.0, 400.0]
        tr0, td0 = (0.8, 5.0) if pol == "inward" else (2.5, 14.0)
        r = least_squares(lambda p: model(t, p[0], p[1], p[2], p[3], p[3] * p[4]) - avg, [0, A0, 0.0, tr0, td0 / tr0],
                          bounds=(lo, hi), max_nfev=300)
        c, A, t0, tr, k = r.x
        td = tr * k
        ts, rise, hw, area = curve_kinetics(tr, td)
        out.update(amp=float(A), rise_tau=float(tr), decay_tau=float(td), rise1090=float(rise), halfwidth=float(hw), charge_fC=float(A * area),
                   r2=float(1 - np.sum(r.fun ** 2) / np.sum((avg - avg.mean()) ** 2)))
    except Exception:
        pass
    return out


def event_stats(evs, dur, pre, cfg=CFG):
    o = OrderedDict()
    n = len(evs)
    a = np.array([e["amp_pA"] for e in evs], float)
    o[pre + "n"] = n
    o[pre + "freq_Hz"] = n / dur if dur > 0 else np.nan
    for nm, fn in (("amp_mean_pA", np.mean), ("amp_median_pA", np.median), ("amp_p10_pA", lambda x: np.percentile(x, 10)),
                   ("amp_p90_pA", lambda x: np.percentile(x, 90)), ("amp_max_pA", np.max)):
        o[pre + nm] = float(fn(a)) if n else np.nan
    o[pre + "amp_sd_pA"] = float(a.std(ddof=1)) if n > 1 else np.nan
    o[pre + "amp_cv"] = float(a.std(ddof=1) / a.mean()) if n > 2 else np.nan
    big = a >= cfg["COMMON_CUTOFF_PA"]
    o[pre + "n_ge_common"] = int(big.sum())
    o[pre + "freq_ge_common_Hz"] = big.sum() / dur if dur > 0 else np.nan
    o[pre + "charge_total_pC"] = float(sum(e["charge_fC"] for e in evs) / 1e3)
    o[pre + "charge_rate_fC_per_s"] = float(sum(e["charge_fC"] for e in evs) / dur) if dur > 0 else np.nan
    iso = [e for e in evs if e["isolated"]]
    o[pre + "n_isolated"] = len(iso)
    o[pre + "frac_overlapped"] = float(np.mean([e["overlapped"] for e in evs])) if n else np.nan
    kin = len(iso) >= cfg["MIN_EVENTS_KINETICS"]
    for nm, k in (("rise_10_90_ms", "rise1090_ms"), ("decay_tau_ms", "td"), ("halfwidth_ms", "halfwidth_ms"), ("rise_tau_ms", "tr")):
        v = np.array([e[k] for e in iso], float)
        o[pre + nm + "_median"] = float(np.median(v)) if kin else np.nan
    ie = np.array([e["iei_prev_ms"] for e in evs], float)
    ie = ie[np.isfinite(ie)]
    o[pre + "iei_median_ms"] = float(np.median(ie)) if len(ie) > 2 else np.nan
    o[pre + "iei_cv"] = float(ie.std(ddof=1) / ie.mean()) if len(ie) > 5 else np.nan
    # stability over time: events per 30 s bin (valid bins only handled by caller via dur) -> first-half / second-half freq
    return o


def time_course(evs, dur_total, bin_s=30.0):
    nb = max(1, int(dur_total // bin_s))
    t = np.array([e["peak_s"] for e in evs])
    h, _ = np.histogram(t, bins=np.arange(nb + 1) * bin_s)
    return h / bin_s


# ------------------------------------------------------------------------------------------------
# file bookkeeping
# ------------------------------------------------------------------------------------------------
PROTO = re.compile(r"^(CapCheck|GapFree|HCN|Firing|AHP)", re.I)


def name_info(path, root=None):
    p = Path(path)
    parts = p.stem.split("_")
    key = p.stem
    for i, s in enumerate(parts):
        if i > 0 and PROTO.match(s):
            key = "_".join(parts[:i])
            break
    run = parts[-1]
    date = parts[0][:8]
    rig = ""
    if root is not None:
        try:
            rel = p.relative_to(root)
            rig = rel.parts[0].split()[0] if len(rel.parts) > 1 else ""
        except ValueError:
            pass
    cond = "TTX" if "TTX" in parts else ("Control" if "Control" in parts else "")
    sex = next((s for s in parts if s in ("F", "M")), "")
    div = next((s[3:] for s in parts if s.startswith("DIV")), "")
    return dict(cell_key=key, run=run, date=date, rig=rig, condition=cond, sex=sex, DIV=div, tags=",".join(s for s in parts if s in ("KCl", "bad")))


def load_holdings(path):
    H = {}
    if not path or not Path(path).exists():
        return H
    with open(path, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            H.setdefault((row["date"], row["run"]), []).append(row)
    return H


def lookup_hold(H, info):
    c = H.get((info["date"], info["run"]), [])
    if not c:
        return None
    if len(c) > 1:
        same = [x for x in c if x["rig"].lower() == info["rig"].lower()]
        c = same or c
    return c[0]


def recording_qc(I, fs):
    n = len(I)
    q = int(fs * 0.01)
    blk = I[: (n // q) * q].reshape(-1, q).mean(axis=1)
    drift = float(np.polyfit(np.arange(len(blk)) * 0.01, blk, 1)[0] * 60.0) if len(blk) > 10 else np.nan
    return OrderedDict(duration_s=n / fs, Ihold_median_pA=float(np.median(I)), Ihold_start_pA=float(np.median(I[: int(5 * fs)])),
                       Ihold_end_pA=float(np.median(I[-int(5 * fs):])), Ihold_drift_pA_per_min=drift)


def analyse_file(path, root=None, H=None, cfg=CFG, min_amp=None, keep_aux=False):
    I, fs, units = load_abf(path)
    info = name_info(path, root)
    hold = lookup_hold(H or {}, info)
    row = OrderedDict(cell_key=info["cell_key"], file=Path(path).name, rig=info["rig"], date=info["date"], run=info["run"],
                      condition=info["condition"], sex=info["sex"], DIV=info["DIV"], tags=info["tags"])
    if hold:
        row.update(hold_true_mV=hold["hold_true_mV"], hold_cmd_mV=hold["hold_cmd_mV"], internal=hold["internal"], hold_confidence=hold["confidence"], hold_source=hold["source"])
    else:
        row.update(hold_true_mV="", hold_cmd_mV="", internal="", hold_confidence="", hold_source="")
    if units != "pA":
        row["status"] = "not analysed: units %s" % units
        return row, [], [], None, None
    row.update(recording_qc(I, fs))
    if hold and hold["source"] == "excluded":
        row["status"] = "not analysed: " + hold["comment"]
        return row, [], [], None, None
    ok, rej, di, aux = detect(I, fs, cfg, min_amp=min_amp, return_aux=True)
    dur = di["duration_s"] - di["masked_s"]
    row["status"] = "analysed"
    row.update(analysed_s=dur, masked_s=di["masked_s"], noise_raw_pA=di["noise_raw_pA"], noise_lp_pA=di["noise_lp_pA"],
               limit_inward_pA=di["limit_inward_pA"], limit_outward_pA=di["limit_outward_pA"],
               n_candidates_inward=di["n_cand_inward"], n_candidates_outward=di["n_cand_outward"],
               n_fitted_inward=di["n_fit_inward"], n_fitted_outward=di["n_fit_outward"])
    reasons = {}
    for r in rej:
        reasons[r["reason"]] = reasons.get(r["reason"], 0) + 1
    row["n_artefact_excluded"] = sum(v for k, v in reasons.items() if k.startswith("excluded"))
    row["n_rejected_after_fit"] = sum(reasons.values())
    row["reject_reasons"] = json.dumps(reasons)
    avgs = {}
    for pol in ("inward", "outward"):
        evs = [e for e in ok if e["polarity"] == pol]
        row.update(event_stats(evs, dur, pol + "_"))
        av = average_event(aux["ylp"], fs, [e for e in evs if e["isolated"]], pol)
        if av:
            avgs[pol] = av
            for k in ("amp", "rise1090", "decay_tau", "halfwidth", "charge_fC", "r2"):
                row[pol + "_avg_" + k] = av.get(k, np.nan)
            row[pol + "_avg_n"] = av["n"]
    allv = ok
    tc = time_course(allv, di["duration_s"])
    row["events_per_30s_min"] = float(tc.min()) * 30 if len(tc) else np.nan
    row["events_per_30s_max"] = float(tc.max()) * 30 if len(tc) else np.nan
    row["freq_first_half_Hz"] = sum(1 for e in allv if e["peak_s"] < di["duration_s"] / 2) / (di["duration_s"] / 2)
    row["freq_second_half_Hz"] = sum(1 for e in allv if e["peak_s"] >= di["duration_s"] / 2) / (di["duration_s"] / 2)
    for e in ok + rej:
        e["cell_key"] = info["cell_key"]
        e["file"] = Path(path).name
    return row, ok, rej, avgs, (aux if keep_aux else None)


EVENT_COLS = ["cell_key", "file", "polarity", "peak_s", "onset_s", "amp_pA", "rough_amp_pA", "snr", "rise1090_ms", "tr", "td", "halfwidth_ms",
              "charge_fC", "resid_ratio", "crit", "template", "limit_pA", "overlapped", "isolated", "truncated", "iei_prev_ms", "baseline_pA", "noise_lp_pA"]
REJ_COLS = ["cell_key", "file", "polarity", "peak_s", "reason", "rough_amp_pA", "amp_pA", "limit_pA", "crit", "rise1090_ms", "td", "resid_ratio"]


def _fmt(x):
    if isinstance(x, (float, np.floating)):
        return "" if not np.isfinite(x) else ("%.5g" % x)
    if isinstance(x, (bool, np.bool_)):
        return int(x)
    return x


def write_csv(path, rows, cols=None):
    if not rows:
        return
    if cols is None:
        cols = []
        for r in rows:
            for k in r:
                if k not in cols:
                    cols.append(k)
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([_fmt(r.get(c, "")) for c in cols])


def _work(args):
    path, root, hpath, min_amp = args
    H = load_holdings(hpath)
    t = time.time()
    try:
        row, ok, rej, avgs, _ = analyse_file(path, root, H, CFG, min_amp)
    except Exception as ex:                      # keep the batch going
        info = name_info(path, root)
        return dict(cell_key=info["cell_key"], file=Path(path).name, status="ERROR: %r" % ex), [], [], None, 0.0
    return row, ok, rej, avgs, time.time() - t


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, help="Ephys folder; every *GapFree*.abf below it is analysed")
    ap.add_argument("--file", type=Path, nargs="*", help="analyse only these files")
    ap.add_argument("--list", type=Path, help="text file with one .abf path per line")
    ap.add_argument("--holdings", type=Path, help="holdings.csv (default: <root>/_traces analysis/holdings.csv)")
    ap.add_argument("--out", type=Path, default=Path("psc_detection_out"))
    ap.add_argument("--min-amp", type=float, default=None, help="override ABS_MIN_PA")
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    files = []
    if a.file:
        files += a.file
    if a.list:
        files += [Path(x.strip()) for x in a.list.read_text().splitlines() if x.strip()]
    if a.root:
        files += sorted(a.root.rglob("*GapFree*.abf"))
    files = sorted({f for f in files})
    if not files:
        sys.exit("no gap-free files found")
    hpath = a.holdings or (a.root / "_traces analysis" / "holdings.csv" if a.root else None)
    a.out.mkdir(parents=True, exist_ok=True)
    jobs = [(f, a.root, hpath, a.min_amp) for f in files]
    rows, evrows, rejrows, avgs_all = [], [], [], {}
    t0 = time.time()
    if a.workers > 1 and len(jobs) > 1:
        from multiprocessing import Pool
        with Pool(a.workers) as pool:
            it = pool.imap(_work, jobs)
            res = list(_progress(it, len(jobs), t0))
    else:
        res = list(_progress(map(_work, jobs), len(jobs), t0))
    for (row, ok, rej, avgs, dt), job in zip(res, jobs):
        rows.append(row)
        evrows += ok
        rejrows += rej
        if avgs:
            for pol, av in avgs.items():
                avgs_all[job[0].name + "|" + pol] = av
    write_csv(a.out / "psc_recordings.csv", rows)
    write_csv(a.out / "psc_events.csv", evrows, EVENT_COLS)
    write_csv(a.out / "psc_candidates_rejected.csv", rejrows, REJ_COLS)
    if avgs_all:
        np.savez_compressed(a.out / "psc_event_averages.npz", **{k + "|" + f: v[f] for k, v in avgs_all.items() for f in ("t_ms", "avg", "sem")})
    print("done: %d files, %d events, %d rejected after fit, %.0f s" % (len(rows), len(evrows), len(rejrows), time.time() - t0))


def _progress(it, n, t0):
    for i, r in enumerate(it, 1):
        if i % 10 == 0 or i == n:
            print("  %d/%d  (%.0f s)" % (i, n, time.time() - t0), flush=True)
        yield r


if __name__ == "__main__":
    main()
