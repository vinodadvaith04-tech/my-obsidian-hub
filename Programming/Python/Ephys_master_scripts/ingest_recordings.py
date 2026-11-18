"""
ingest_recordings.py - bring recordings from someone else (or a new session) into the master: rename them to the
naming convention, add them to the master workbook, and keep the old names.

Convention (same as every renamed file so far):
    <date>_<S.No>_<T.Code>_<sex>_slice<N>_C<N>_<counter>.abf        e.g. 20261110_Adv00_NT_unk_slice2_C1_057.abf
    date     = date in the ABF header
    slice N  = slices of that DATE numbered in order of first recording; numbering continues after the slices already in the master for that date
    C N      = cells of that slice numbered in order of first recording
    counter  = running file number of that date (continues after the counters already used that date, so names never clash)
    S.No / T.Code default to Adv00 / NT (the values used for Sacha, Gerda and the Cs pilot); sex from the file name if it says so, else unk.
    The original ids (slicing date, slice, cell, original name) stay in the master (old_name, cells.cell_old_folder, analysis_note).

What it does
  1. finds the .abf files under --folder (relative to the Ephys folder) that are not in the master yet
  2. reads each header (date, time, protocol, duration, sweeps) and the recorded voltage monitor (channel 2) = the real holding
  3. dry run: prints the plan + writes rename_plan.csv. With --apply: renames the files in place, adds the rows to sheets
     files / cells / rename_log (backup of the workbook first) and sets use_for_miniML = Y for gap-free sPSC recordings
  Nothing is overwritten; a name that already exists makes the script stop.

Usage
  python scripts/ingest_recordings.py --folder "3. KCl IPSC recordings" --session "Helena KCl IPSC" --internal KCl \\
         --psc-content "KCl IPSC recording (blockers not confirmed)" [--drug "DNQX+AP5"] [--dry-run] [--apply]
Needs: pyabf pandas openpyxl
"""
import argparse, re, sys
from pathlib import Path
import numpy as np, pandas as pd, pyabf

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ephys_master as em

SEX = {"male": "male", "m": "male", "female": "female", "f": "female"}


def parse_old(name):
    stem = name[:-4]
    cnt = re.search(r"_(\d{3,4})$", stem)
    cell = re.search(r"_[cC](\d+)(?=_)", stem + "_")
    sx = re.search(r"_(male|female|M|F)_", stem, re.I)
    sl = re.search(r"[Ss]lice(\d*)", stem.replace("Slicing", "Sxxxing"))
    sd = re.search(r"(\d{4})(?=slice)", stem, re.I) or re.search(r"Slicing(\d{4})", stem)
    return dict(counter_old=int(cnt.group(1)) if cnt else None, cell_old=cell.group(1) if cell else None,
                sex=SEX.get(sx.group(1).lower()) if sx else None, slice_old=sl.group(1) if sl else None,
                slicing=sd.group(1) if sd else None)


def read_header(p):
    ab = pyabf.ABF(str(p), loadData=False)
    proto = ab.protocol or ""
    ptype = "sPSC gap-free" if "GapFree" in proto else ("cap check" if "Cap" in proto else ("spike frequency" if "spike" in proto else "other"))
    try: hh = float(ab.holdingCommand[0])
    except Exception: hh = None
    h = dict(date=int(ab.abfDateTime.strftime("%Y%m%d")), time=ab.abfDateTime.strftime("%H:%M:%S"), protocol=proto, ptype=ptype,
             mode="VC" if ab.adcUnits[0] == "pA" else "CC", dur=round(ab.dataLengthSec, 1), sweeps=ab.sweepCount, hold_header=hh,
             vm=None, vm_note="")
    if ab.channelCount > 1 and ptype == "sPSC gap-free":      # second channel = voltage monitor of the amplifier
        ab2 = pyabf.ABF(str(p)); ab2.setSweep(0, channel=1)
        med = float(np.median(ab2.sweepY[: int(60 * ab2.dataRate)])); unit = ab2.adcUnits[1]; nm = ab2.adcNames[1]
        if unit == "mV":
            h["vm"], h["vm_note"] = med, f"{nm} median {med:.2f} mV"
        elif unit == "pA" and -1500 < med < -300:             # monitor stored with the wrong unit and x20 gain (value / 20 = mV)
            h["vm"], h["vm_note"] = med / 20, f"{nm} labelled pA, median {med:.1f} (/20 = {med / 20:.1f} mV)"
        else:
            h["vm_note"] = f"{nm} {unit} median {med:.2f}: not read as a voltage"
    return h


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", type=Path, default=em.default_master())
    ap.add_argument("--data-root", type=Path, default=None, help="Ephys folder that --folder is relative to (default: the master's folder)")
    ap.add_argument("--folder", required=True)
    ap.add_argument("--session", required=True)
    ap.add_argument("--internal", required=True, help="KCl, KGlu or Cs-gluconate")
    ap.add_argument("--drug", default="")
    ap.add_argument("--psc-content", default="")
    ap.add_argument("--sno", default="Adv00")
    ap.add_argument("--tcode", default="NT")
    ap.add_argument("--holding", type=float, default=None, help="command holding for all (mV); default = read from the voltage monitor channel, else -60")
    ap.add_argument("--note", default="")
    ap.add_argument("--plan-out", type=Path, default=None)
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    root = a.data_root or a.master.parent
    folder = root / a.folder
    files = em.read_sheet(a.master, "files")
    known_new, known_old = set(files.new_name.dropna()), set(files.old_name.dropna())
    known_paths = set(files.path_on_mac.dropna().astype(str))
    new = []
    for p in sorted(folder.rglob("*.abf")):
        rel = p.relative_to(root).as_posix()
        if p.name.startswith("._") or not p.is_file(): continue
        if p.name in known_new or p.name in known_old or rel in known_paths: continue
        new.append(p)
    print(f"{len(new)} new recordings under {a.folder}")
    if not new: return
    info = []
    for p in new:
        h = read_header(p); h.update(parse_old(p.name)); h["path"] = p; info.append(h)
    info.sort(key=lambda h: (h["date"], h["time"]))

    # numbering per date, continuing after what the master already uses that date
    fm = files.copy()
    fm["_date"] = pd.to_numeric(fm["date_header"], errors="coerce")
    rows, plan, slice_idx, cell_idx, next_ctr, cell_seq = [], [], {}, {}, {}, {}
    taken = set(known_new)
    for h in info:
        d = h["date"]
        used = fm[fm["_date"] == d]
        if d not in next_ctr:
            next_ctr[d] = int(pd.to_numeric(used["counter"], errors="coerce").max() or 0) + 1 if len(used) else 1
            cell_seq[d] = int(pd.to_numeric(used["slice"], errors="coerce").max() or 0) if len(used) else 0
        skey = (d, h["slicing"], h["slice_old"])
        if skey not in slice_idx:
            cell_seq[d] += 1; slice_idx[skey] = cell_seq[d]
        si = slice_idx[skey]
        ckey = (skey, h["cell_old"])
        if ckey not in cell_idx:
            cell_idx[ckey] = 1 + sum(1 for k in cell_idx if k[0] == skey)
        ci = cell_idx[ckey]
        ctr = next_ctr[d]; next_ctr[d] += 1
        sex = h["sex"] or "unk"
        nn = f"{d}_{a.sno}_{a.tcode}_{sex}_slice{si}_C{ci}_{ctr:03d}.abf"
        if nn in taken: sys.exit(f"name collision: {nn}")
        taken.add(nn)
        old_rel = h["path"].relative_to(root).as_posix()
        new_rel = (h["path"].parent / nn).relative_to(root).as_posix()
        # holding: from the recorded voltage monitor, else the given value, else -60 (standing rule)
        if a.holding is not None:
            hold, flag = a.holding, f"stated: {a.holding:g} mV for all recordings of this session"
        elif h["vm"] is not None:
            r5 = round(h["vm"] / 5) * 5
            hold = float(r5) if abs(h["vm"] - r5) <= 1.5 else round(h["vm"], 1)
            flag = f"from the recorded voltage monitor (channel 2): {h['vm_note']}; header holding is 0 and not informative"
        elif h["ptype"] == "sPSC gap-free":
            hold, flag = -60.0, "default -60 (no note, no usable voltage monitor)"
        else:
            hold, flag = None, None
        gap = h["ptype"] == "sPSC gap-free"
        rows.append(dict(new_name=nn, old_name=h["path"].name, old_folder=h["path"].parent.relative_to(root).as_posix(),
                         new_folder=h["path"].parent.relative_to(root).as_posix(), path_on_mac=new_rel, status="renamed",
                         date_header=d, start_time=h["time"], session=a.session, **{"S.No": a.sno, "T.Code": a.tcode},
                         sex=sex, sex_source="file name" if h["sex"] else "not stated", slice=si, cell=ci, counter=ctr,
                         protocol=h["protocol"], protocol_type=h["ptype"], mode=h["mode"], duration_s=h["dur"], sweeps=h["sweeps"],
                         psc_content=a.psc_content if gap else None, drug=(a.drug or None) if gap else None, internal=a.internal,
                         internal_source="stated (folder name)", holding_header_mV=h["hold_header"],
                         holding_command_mV=hold, holding_flag=flag, use_for_miniML="Y" if gap else "N",
                         miniML_directions=("negative" if a.internal == "KCl" else None) if gap else None,
                         analysis_note=(f"original ids: slicing {h['slicing'] or '?'}, slice {h['slice_old'] or '?'}, cell c{h['cell_old'] or '?'}. " + a.note).strip()))
        plan.append(dict(old=old_rel, new=new_rel, session=a.session, header_date=d, header_time=h["time"], protocol=h["protocol"],
                         slicing=h["slicing"], slice_old=h["slice_old"], cell_old=h["cell_old"], holding=hold, vm_note=h["vm_note"]))
    pl = pd.DataFrame(plan)
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 80)
    print(pl[["old", "new", "holding"]].assign(old=lambda x: x.old.str.split("/").str[-1], new=lambda x: x.new.str.split("/").str[-1]).to_string())
    out = a.plan_out or (a.master.parent / "_work" / "rename_plan.csv")
    out.parent.mkdir(parents=True, exist_ok=True); pl.to_csv(out, index=False); print("plan written to", out)
    if not a.apply:
        print("dry run (add --apply to rename and write the master)"); return

    for h, r in zip(info, rows):                        # rename in place (never over an existing file)
        dst = h["path"].with_name(r["new_name"])
        if dst.exists(): sys.exit(f"{dst} exists")
    for h, r in zip(info, rows):
        h["path"].rename(h["path"].with_name(r["new_name"]))
    newf = pd.DataFrame(rows)
    # cells rows
    cells = []
    for (skey, cold), ci in cell_idx.items():
        g = newf[(newf["date_header"] == skey[0]) & (newf["slice"] == slice_idx[skey]) & (newf["cell"] == ci)]
        gf = g[g.protocol_type == "sPSC gap-free"]
        cells.append({"date": skey[0], "session": a.session, "slice": slice_idx[skey], "cell_new": ci, "cell_old_folder": cold,
                      "T.Code": a.tcode, "sex": g["sex"].iloc[0], "internal": a.internal, "internal_source": "stated (folder name)",
                      "n_files": len(g), "n_sPSC_gapfree": len(gf),
                      "holding_commands (gap-free, time order)": ", ".join(f"{v:g}" for v in gf["holding_command_mV"].dropna()),
                      "flags": None})
    cf = em.merge_frames(em.read_sheet(a.master, "cells"), pd.DataFrame(cells), em.KEYS["cells"])
    ff = em.merge_frames(files, newf, em.KEYS["files"])
    rl = em.read_sheet(a.master, "rename_log")
    rl = pd.concat([rl, pd.DataFrame({"old_path": pl["old"], "new_path": pl["new"]})], ignore_index=True)
    em.write_sheets(a.master, {"files": ff, "cells": cf, "rename_log": rl})
    print(f"renamed {len(rows)} files, added {len(newf)} rows to files, {len(cells)} to cells, {len(pl)} to rename_log")


if __name__ == "__main__":
    main()
