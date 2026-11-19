"""
One-off, 18 Nov 2026: manual holdings / notes for the Cs slice 5-6 + main-experiment batch that step 1 could not (or should not) settle alone.
Run once, after step 1:   python scripts/one_off/fill_holdings_20261118.py --master <Ephys_master.xlsx>
It edits only the files sheet (holding_command_mV, holding_flag, analysis_note) of rows from 17/18 Nov. Backup first.
"""
import argparse, sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import ephys_master as em

ap = argparse.ArgumentParser(); ap.add_argument("--master", type=Path, default=em.default_master()); a = ap.parse_args()
f = em.read_sheet(a.master, "files")
f["_d"] = pd.to_numeric(f["date_header"], errors="coerce")


def rows(date, sl, ce, gap=True, minimum=30):
    m = (f["_d"] == date) & (pd.to_numeric(f["slice"], errors="coerce") == sl) & (pd.to_numeric(f["cell"], errors="coerce") == ce)
    if gap:
        m &= (f["protocol_type"] == "sPSC gap-free") & (pd.to_numeric(f["duration_s"], errors="coerce") >= minimum)
    return f[m].sort_values("start_time").index.tolist()


def add_note(i, txt, col="analysis_note"):
    cur = f.at[i, col]
    f.at[i, col] = txt if pd.isna(cur) or not str(cur).strip() else f"{cur}; {txt}"


# ---- Slice 6 (journal slice 2) cell 3: 19 gap-free recordings, journal lists 9 holdings = 18 files -------------------------
ix = rows(20261117, 6, 3)
assert len(ix) == 19, len(ix)
hold = [-60] + [-60, -60] + [-70, -70, -80, -80, -50, -50, -40, -40, -30, -30, 0, 0, 10, 10, -70, -70]
base = ("SET BY HAND 18 Nov (not from the automatic journal match): the journal lists 9 holdings (-60,-70,-80,-50,-40,-30,0,10,-70) but there are 10 "
        "groups of recordings; ")
for k, (i, h) in enumerate(zip(ix, hold)):
    f.at[i, "holding_command_mV"] = float(h); f.at[i, "holding_listed_in_note"] = float(h) if k else None
    if k == 0:
        f.at[i, "holding_flag"] = (base + "this 120-s file belongs to an extra first attempt (its partner stopped after 11.6 s, then a 40 min pause); no holding in the journal, "
                                   "so the default -60 mV (files without a note) was used; baseline current about -90 pA fits -60")
    elif k in (1, 2):
        f.at[i, "holding_flag"] = (base + "this pair (after the pause, two cap checks) was taken as the journal's first -60 (baseline current -70/-82 pA); "
                                   "the next 8 pairs follow the journal order. CHECK the assumption")
    else:
        f.at[i, "holding_flag"] = base + "pairs 2-10 follow the journal order from -70 onwards"
        if h in (-40, -30, -50):
            add_note(i, "baseline holding current is not monotonic with holding in this cell (-50: -39 pA, -40: -78 pA, -30: -64 pA): irregular cell, check holdings and quality")
        # keep holding_listed_in_note as the journal value for k>=3
        f.at[i, "holding_listed_in_note"] = float(h)

# ---- Slice 5 cell 1: the -40 pair has a baseline current like the -60 pair ---------------------------------------------
ix = rows(20261117, 5, 1)
for i in ix:
    if f.at[i, "holding_command_mV"] == -40:
        add_note(i, "CHECK holding: journal says -40 mV, but the baseline holding current (-121/-104 pA) is more negative than at -50 mV (-86/-114 pA) and similar to -60 mV (-131 pA); drift or note order?")

# ---- Slice 6 cell 2: artifact pair at -40 mV -----------------------------------------------------------------------------
ix = [i for i in rows(20261117, 6, 2) if f.at[i, "holding_command_mV"] == -40]
for i in ix:
    add_note(i, "contains large transients/artifact (5th percentile about -860 pA in the first file, SD 51 pA in the second): inspect before analysing")

# ---- main experiment (slice 6): journal remarks ---------------------------------------------------------------------------
for i in rows(20261118, 6, 1):
    add_note(i, "holding current about -1 nA (very leaky); journal: not the best cell, lost after these two recordings")
for i in rows(20261118, 6, 3):
    add_note(i, "journal: continuous firing throughout, a lot of spiking at -40/-50; remove from the final data set if enough N")
for i in rows(20261118, 6, 2):
    if f.at[i, "holding_command_mV"] == 0:
        add_note(i, "journal: many EPSCs visible at 0 mV")
for i in rows(20261118, 6, 2)[-2:]:
    add_note(i, "journal: access resistance rose to 34 MOhm at the end, recording stopped")

f = f.drop(columns=["_d"])
em.write_sheets(a.master, {"files": f})
print("files sheet updated")

# keep the cells sheet in step with the files sheet for the cell set by hand
c = em.read_sheet(a.master, "cells")
ff = em.read_sheet(a.master, "files")
m = (pd.to_numeric(c["date"], errors="coerce") == 20261117) & (pd.to_numeric(c["slice"], errors="coerce") == 6) & (pd.to_numeric(c["cell_new"], errors="coerce") == 3)
g = ff[(pd.to_numeric(ff["date_header"], errors="coerce") == 20261117) & (pd.to_numeric(ff["slice"], errors="coerce") == 6) &
       (pd.to_numeric(ff["cell"], errors="coerce") == 3) & (ff["protocol_type"] == "sPSC gap-free")].sort_values("start_time")
c.loc[m, "holding_commands (gap-free, time order)"] = ", ".join(f"{v:g}" for v in g["holding_command_mV"].dropna())
c.loc[m, "flags"] = ("holdings set by hand 18 Nov (journal lists 9 holdings, 10 groups of recordings exist; see holding_flag); "
                     "baseline current irregular at -50/-40/-30 mV; one extra 120-s recording before a 40 min pause")
em.write_sheets(a.master, {"cells": c}, backup=False)
print("cells sheet updated")
