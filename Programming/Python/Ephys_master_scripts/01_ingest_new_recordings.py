#!/usr/bin/env python3
"""
01_ingest_new_recordings.py  -  STEP 1 of the pipeline: bring new recordings into the master workbook.

What it does, in this order
  1. checks the master is not open in Excel; as soon as it is clear that something will change (and before ANY file or sheet is
     touched) it copies the master to  _backups/pipeline/
  2. looks for .abf files under the Ephys folder that are not in the master's `files` sheet yet
     (a known file that you moved to another folder only gets its new location written)
  3. reads each header (date, time, protocol, duration, sweeps)
  4. checks the file name against the convention
         <date>_Adv<S.No>_<T.Code>_<sex>_slice<N>_C<N>_<counter>.abf      e.g. 20261118_Adv6_1B12_male_slice6_C2_0003.abf
     - follows it  -> the name is kept (the name is compared with the lab journal and differences are flagged, never changed)
     - does not    -> renamed ON DISK, in place (never moved, never over-writing); S.No / T.Code / sex come from the
                      lab journal of that day if it has them, otherwise Adv00 / NT / unk; slice and cell come from the folder
                      names (Slice 5 / Cell 2) or the old file name; the counter continues after the day's highest counter
  5. reads that day's lab journal (Obsidian, Lab Journal/Ephys/...) and writes the holding of every gap-free recording:
     the "cap check / -80 / -80 / cap check / -70 / -70 ..." list is matched to the recordings in order, ONLY when the
     number of listed holdings equals the number of gap-free recordings of that cell. Otherwise the holding stays BLANK
     (and holding_flag says why) so you can type it in Excel. Internal solution and drug are read from lines like
        Internal: K-gluconate        Drug: gabazine
     at the top of the day's journal or in a slice header block; or given with --internal / --drug.
  6. writes old name, new name, old / new location, protocol, holding ... to `files`, updates `cells` and `rename_log`
  7. writes a report (what was found, renamed, left blank and why) to logs/ and prints it

It does NOT run miniML or the second detector: that is  02_run_detection.py.
Rows with a blank holding or internal are simply skipped by step 2 until you fill them in (then run step 2 again).
Re-running step 1 later also tries to fill still-blank holdings from the journal (e.g. after you finished the notes).

Usage (terminal, from the Ephys folder or anywhere):
    python scripts/01_ingest_new_recordings.py --dry-run            # show everything, change nothing
    python scripts/01_ingest_new_recordings.py                      # do it
  Options: --folder "1. Sacha recordings/Caesium pilot"   only look there
           --internal Cs-gluconate  --drug gabazine  --session "Cs internal pilot"   for every new file of this run
           --journal-dir "<Lab Journal folder>"   --master <file>   --data-root <folder>
           --undo-renames logs/renamed_XXXX.csv     put the file names back (then restore the master from _backups/pipeline/)
Needs: pandas openpyxl pyabf
"""
from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ephys_master as em
import journal_parser as jp
import pipeline_common as pc

NAME_RE = re.compile(r"^(?P<date>\d{8})_(?P<sno>Adv\d+)_(?P<tcode>[A-Za-z0-9]+)_(?P<sex>male|female|unk)_slice(?P<slice>\d+)_C(?P<cell>\d+)_(?P<ctr>\d{3,4})\.abf$", re.I)
FILE_COLS = ["new_name", "old_name", "old_folder", "new_folder", "path_on_mac", "status", "date_header", "start_time", "session",
             "S.No", "T.Code", "sex", "sex_source", "slice", "cell", "counter", "protocol", "protocol_type", "mode", "duration_s",
             "sweeps", "psc_content", "drug", "internal", "internal_source", "holding_header_mV", "holding_listed_in_note",
             "ljp_in_note", "holding_command_mV", "holding_flag", "use_for_miniML", "blockers", "Ra", "Rm", "Cm", "Vm",
             "miniML_directions", "analysis_note"]
CELL_COLS = ["date", "session", "slice", "cell_new", "cell_old_folder", "T.Code", "sex", "internal", "internal_source",
             "n_files", "n_sPSC_gapfree", "holding_commands (gap-free, time order)", "flags"]


# ---------------------------------------------------------------------------------------------- headers
def read_header(p: Path, proto_map: dict) -> dict:
    import pyabf
    ab = pyabf.ABF(str(p), loadData=False)
    proto = ab.protocol or ""
    mode = "VC" if (ab.adcUnits[0] == "pA") else "CC"
    if proto in proto_map:
        ptype = proto_map[proto]
    elif "GapFree" in proto:
        ptype = "sPSC gap-free"
    elif "Cap" in proto:
        ptype = "cap check"
    elif "spike" in proto.lower():
        ptype = "spike frequency"
    elif "HCN" in proto:
        ptype = f"HCN ({mode})"
    elif "AHP" in proto:
        ptype = "AHP"
    else:
        ptype = "other"
    try:
        hh = float(ab.holdingCommand[0])
    except Exception:
        hh = None
    h = dict(date=int(ab.abfDateTime.strftime("%Y%m%d")), time=ab.abfDateTime.strftime("%H:%M:%S"), protocol=proto, ptype=ptype,
             mode=mode, dur=round(ab.dataLengthSec, 1), sweeps=ab.sweepCount, hold_header=hh, vm=None)
    if ab.channelCount > 1 and ptype == "sPSC gap-free":      # second channel = voltage monitor (only trusted when it reads in mV)
        try:
            ab2 = pyabf.ABF(str(p)); ab2.setSweep(0, channel=1)
            med = float(np.median(ab2.sweepY[: int(60 * ab2.dataRate)]))
            if ab2.adcUnits[1] == "mV":
                h["vm"] = med
        except Exception:
            pass
    return h


def old_ids(rel: str, name: str):
    """slice / cell numbers from folder names (Slice 5 / Cell 2) or the old file name (..._slice2_c3_...)."""
    sl = ce = None
    for part in rel.split("/")[:-1]:
        m = re.search(r"slice\s*_?(\d+)", part, re.I)
        if m: sl = int(m.group(1))
        m = re.search(r"cell\s*_?(\d+)", part, re.I)
        if m: ce = int(m.group(1))
    stem = name[:-4]
    if sl is None:
        m = re.search(r"slice\s*_?(\d+)", stem, re.I)
        sl = int(m.group(1)) if m else None
    if ce is None:
        m = re.search(r"_c(\d+)(?=_|$)", stem, re.I)
        ce = int(m.group(1)) if m else None
    return sl, ce


def psc_content_for(internal, drug):
    d = (drug or "").lower()
    if not internal:
        return None
    if internal == "KCl":
        return "IPSC (DNQX+AP5, KCl)" if ("dnqx" in d or "ap5" in d) else "mixed (KCl, no blockers)"
    if "gabazine" in d or "sr95531" in d:
        return "EPSC (gabazine)"
    if "dnqx" in d or "ap5" in d:
        return "IPSC (DNQX+AP5)"
    return "mixed E+I" if not d else f"mixed E+I ({drug})"


def fnum(v):
    try:
        x = float(v)
        return None if np.isnan(x) else x
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------------------------- main
def run(a, log):
    cfg = pc.load_config(a.config)
    root, mp = pc.resolve_paths(cfg, a.master, a.data_root)
    jdir = Path(a.journal_dir or cfg["journal_dir"]).expanduser()
    min_gf = float(cfg["min_gapfree_s"])
    stamp = time.strftime("%Y%m%d_%H%M%S")
    log(f"STEP 1 - ingest new recordings   {time.strftime('%Y-%m-%d %H:%M')}")
    log(f"  Ephys folder : {root}")
    log(f"  master       : {mp}")
    log(f"  lab journal  : {jdir} ({'found' if jdir.exists() else 'NOT FOUND: holdings will stay blank'})")
    pc.check_master_free(mp)
    if a.dry_run:
        log("  DRY RUN: nothing will be changed")

    files = em.read_sheet(mp, "files")
    cells_old = em.read_sheet(mp, "cells")
    rlog = em.read_sheet(mp, "rename_log")
    for c in FILE_COLS:
        if c not in files.columns:
            files[c] = np.nan
    proto_map = {}
    if "protocol" in files.columns:
        pm = files.dropna(subset=["protocol", "protocol_type"]).drop_duplicates("protocol")
        proto_map = dict(zip(pm.protocol, pm.protocol_type))

    # ---- 1. find new files -------------------------------------------------------------------
    base = root / a.folder if a.folder else root
    if not base.exists():
        raise pc.UserProblem(f"Folder to scan does not exist: {base}", "Check --folder (relative to the Ephys folder).")
    skip = set(cfg["skip_dirs"])
    disk = []
    for p in sorted(base.rglob("*.abf")):
        rel_parts = p.relative_to(root).parts
        if p.name.startswith("._") or not p.is_file() or (set(rel_parts[:-1]) & skip):
            continue
        disk.append(p)
    known_paths = set(files["path_on_mac"].dropna().astype(str))
    by_new = files.groupby(files["new_name"].astype(str)).groups
    by_old = files.groupby(files["old_name"].astype(str)).groups
    new_files, moved = [], []
    for p in disk:
        rel = p.relative_to(root).as_posix()
        if rel in known_paths:
            continue
        idxs = list(by_new.get(p.name, [])) + [i for i in by_old.get(p.name, []) if i not in by_new.get(p.name, [])]
        if len(idxs) == 1 and not (root / str(files.at[idxs[0], "path_on_mac"])).exists():
            moved.append((idxs[0], rel))
        else:
            new_files.append(p)
    log(f"\n  .abf files on disk (scanned): {len(disk)} | already in the master: {len(disk) - len(new_files) - len(moved)} "
        f"| NEW: {len(new_files)} | known but moved: {len(moved)}")
    for i, rel in moved:
        old = files.at[i, "path_on_mac"]
        files.at[i, "path_on_mac"] = rel
        files.at[i, "new_folder"] = str(Path(rel).parent.as_posix())
        log(f"    location updated: {files.at[i, 'new_name']}: {old} -> {rel}")

    # ---- 2. headers ---------------------------------------------------------------------------
    recs, problems = [], []
    for p in new_files:
        rel = p.relative_to(root).as_posix()
        try:
            h = read_header(p, proto_map)
        except Exception as e:
            problems.append((rel, f"header unreadable: {type(e).__name__}: {e}"))
            continue
        m = NAME_RE.match(p.name)
        sl, ce = (int(m.group("slice")), int(m.group("cell"))) if m else old_ids(rel, p.name)
        recs.append(dict(path=p, rel=rel, name=p.name, h=h, m=m, slice=sl, cell=ce, date=h["date"], time=h["time"], new=True))
    recs.sort(key=lambda r: (r["date"], r["time"]))
    for rel, msg in problems:
        log(f"    PROBLEM {rel}: {msg}  (left out; fix or remove the file, then run again)")

    # ---- 3. journals --------------------------------------------------------------------------
    jcache = {}
    def journals_for(date):
        if date not in jcache:
            out = []
            for pth in jp.find_journals(jdir, str(date)):
                try:
                    out.append(jp.parse_journal(pth))
                except Exception as e:
                    log(f"    journal {pth.name} could not be read: {e}")
            jcache[date] = out
        return jcache[date]

    fm = files.copy()
    fm["_date"] = pd.to_numeric(fm["date_header"], errors="coerce")
    fm["_slice"] = pd.to_numeric(fm["slice"], errors="coerce")
    fm["_cell"] = pd.to_numeric(fm["cell"], errors="coerce")

    # group = (date, slice, cell): existing master rows + new files
    def group_records(date, sl, ce):
        out = []
        ex = fm[(fm["_date"] == date) & (fm["_slice"] == sl) & (fm["_cell"] == ce)]
        for i, r in ex.iterrows():
            out.append(dict(kind="old", idx=i, time=str(r["start_time"]), ptype=r["protocol_type"], dur=fnum(r["duration_s"]) or 0,
                            hold=fnum(r["holding_command_mV"]), mode=r["mode"]))
        for r in recs:
            if r["date"] == date and r["slice"] == sl and r["cell"] == ce:
                out.append(dict(kind="new", rec=r, time=r["time"], ptype=r["h"]["ptype"], dur=r["h"]["dur"], hold=None, mode=r["h"]["mode"]))
        out.sort(key=lambda x: x["time"])
        return out

    def eligible(x):
        return x["ptype"] == "sPSC gap-free" and x["dur"] >= min_gf and x["mode"] == "VC"

    def slices_of_day(date):
        ex = fm[fm["_date"] == date][["_slice", "start_time"]].dropna(subset=["_slice"])
        pairs = [(int(s), str(t)) for s, t in zip(ex["_slice"], ex["start_time"])] + [(r["slice"], r["time"]) for r in recs if r["date"] == date and r["slice"] is not None]
        first = {}
        for s, t in pairs:
            first[s] = min(t, first.get(s, "99:99:99"))
        return [s for s, _ in sorted(first.items(), key=lambda kv: kv[1])]

    day_groups = {}
    for r in recs:
        day_groups.setdefault(r["date"], set()).add((r["slice"], r["cell"]))
    # existing gap-free rows still without a holding: try to fill them too (journal may be complete now)
    if a.refill_blank:
        blank = fm[(fm["use_for_miniML"].astype(str).str.upper() == "Y") & fm["holding_command_mV"].isna() & fm["_date"].notna() & fm["_slice"].notna() & fm["_cell"].notna()]
        for _, r in blank.iterrows():
            day_groups.setdefault(int(r["_date"]), set()).add((int(r["_slice"]), int(r["_cell"])))

    assign = {}          # (id) -> dict(hold, flag, ...)   id = ("new", path) or ("old", idx)
    meta = {}            # (date, slice) -> dict(journal slice info)
    cell_flags = {}      # (date, slice, cell) -> [text]
    cell_report = []
    for date, gset in sorted(day_groups.items(), key=lambda kv: kv[0]):
        jdays = journals_for(date)
        sl_order = slices_of_day(date)
        for sl in sorted({s for s, _ in gset if s is not None}):
            cells_here = sorted({c for s, c in gset if s == sl and c is not None})
            cands = []                                       # (journal day, jslice, how)
            for jd in jdays:
                for js in jd.slices:
                    if js.number == sl:
                        cands.append((jd, js, "slice number"))
            rank = sl_order.index(sl) if sl in sl_order else None
            for jd in jdays:
                if rank is not None and rank < len(jd.slices) and (jd, jd.slices[rank], "slice number") not in cands:
                    if jd.slices[rank].number != sl:
                        cands.append((jd, jd.slices[rank], f"position (journal slice {jd.slices[rank].number} = {rank + 1}. slice of the day)"))
            def score(c):
                jd, js, how = c
                sc = 0
                for ce in cells_here:
                    jc = next((x for x in js.cells if x.number == ce), None)
                    if jc and not jc.suspicious and len(jc.lines) == sum(eligible(x) for x in group_records(date, sl, ce)):
                        sc += 1
                return sc
            best = max(cands, key=lambda c: (score(c), c[2] == "slice number"), default=None) if cands else None
            meta[(date, sl)] = dict(best=best)
            for ce in cells_here:
                grp = group_records(date, sl, ce)
                el = [x for x in grp if eligible(x)]
                key = (date, sl, ce)
                fl = cell_flags.setdefault(key, [])
                if not el:
                    continue
                why = None
                if not jdays:
                    why = f"no lab journal found for {date}"
                elif best is None:
                    why = f"no journal slice matches file slice {sl}"
                else:
                    jd, js, how = best
                    jc = next((x for x in js.cells if x.number == ce), None)
                    if jc is None:
                        why = f"journal ({jd.path.name}) has no cell {ce} in slice {js.number}"
                    elif jc.suspicious:
                        why = f"journal list for cell {ce} does not look like holdings (numbered list / out of range)"
                    elif len(jc.lines) != len(el):
                        why = (f"journal lists {len(jc.lines)} files ({len(jc.groups)} cap-check groups) but {len(el)} gap-free recordings "
                               f">= {min_gf:g} s exist")
                    else:
                        src = f"{jd.path.name}, journal slice {js.number}, cell {ce}; recordings and listed holdings matched one by one in order; slice matched by {how}"
                        for x, hval in zip(el, jc.lines):
                            if x["hold"] is None:
                                assign[("new", x["rec"]["rel"]) if x["kind"] == "new" else ("old", x["idx"])] = dict(hold=hval, flag=f"from lab journal {src}")
                        if jc.notes:
                            fl.append(f"journal note: {jc.notes}")
                        cell_report.append((date, sl, ce, len(el), "OK", src))
                if why:
                    for x in el:
                        k = ("new", x["rec"]["rel"]) if x["kind"] == "new" else ("old", x["idx"])
                        if x["hold"] is None:
                            assign[k] = dict(hold=None, flag=f"NOT SET (fill in by hand): {why}")
                    fl.append(f"holding not set: {why}")
                    cell_report.append((date, sl, ce, len(el), "BLANK", why))

    # ---- 4. names + rows ------------------------------------------------------------------------
    taken = set(files["new_name"].dropna().astype(str)) | {p.name for p in disk}
    next_ctr = {}
    plan, rows = [], []
    for r in recs:
        date, h, m = r["date"], r["h"], r["m"]
        best = (meta.get((date, r["slice"])) or {}).get("best")
        jd, js = (best[0], best[1]) if best else (None, None)
        jsex = jtc = jser = jint = jdrug = None
        if js:
            jsex = jd.slice_header(js, "sex"); jtc = (jd.slice_header(js, "tcode") or "").replace("&", "") or None
            jser = jd.slice_header(js, "serial"); jint = jd.slice_header(js, "internal"); jdrug = jd.slice_header(js, "drug")
        else:
            jl = journals_for(date)
            if jl:
                jint = jl[0].day_header.get("internal"); jdrug = jl[0].day_header.get("drug")
        notes = []
        if m:
            new_name, status = r["name"], "kept (follows the naming convention)"
            sno, tcode, sex, ctr = m.group("sno"), m.group("tcode"), m.group("sex").lower(), int(m.group("ctr"))
            sex_src = "file name"
            if int(m.group("date")) != date:
                notes.append(f"NAME DATE {m.group('date')} differs from the date in the file header {date}")
            if js:
                if jser and re.sub(r"\D", "", sno) != re.sub(r"\D", "", str(jser)):
                    notes.append(f"CHECK NAME: S.No in the name is {sno} but the lab journal says serial number {jser} for slice {r['slice']}")
                if jtc and jtc.lower() != tcode.lower():
                    notes.append(f"CHECK NAME: T.Code in the name is {tcode} but the lab journal says {jtc}")
                if jsex and jsex != sex:
                    notes.append(f"CHECK NAME: sex in the name is {sex} but the lab journal says {jsex}")
        else:
            sno = f"Adv{int(float(jser))}" if jser and str(jser).replace(".", "").isdigit() else cfg["default_S.No"]
            tcode = jtc or cfg["default_T.Code"]
            sex = jsex or cfg["default_sex"]
            sex_src = "lab journal" if jsex else "not stated"
            if date not in next_ctr:
                used = pd.to_numeric(fm.loc[fm["_date"] == date, "counter"], errors="coerce")
                next_ctr[date] = int(used.max()) + 1 if used.notna().any() else 1
            if r["slice"] is None or r["cell"] is None:
                new_name, status, ctr = r["name"], "NEW - slice/cell number not found in folder or file name, NOT renamed", None
                notes.append("slice/cell unknown: put the file in a 'Slice N/Cell N' folder or fill slice and cell by hand, then rename")
            else:
                while True:
                    ctr = next_ctr[date]; next_ctr[date] += 1
                    new_name = f"{date}_{sno}_{tcode}_{sex}_slice{r['slice']}_C{r['cell']}_{ctr:03d}.abf"
                    if new_name not in taken:
                        break
                status = "renamed"
            mo = re.search(r"(\d{8})", r["name"])
            if mo and int(mo.group(1)) != date:
                notes.append(f"old name carried date {mo.group(1)}, the file header says {date}, so the header date was used in the new name")
        taken.add(new_name)
        rel_dir = str(Path(r["rel"]).parent.as_posix())
        new_rel = f"{rel_dir}/{new_name}"
        if new_name != r["name"]:
            plan.append((r["path"], r["path"].with_name(new_name), r["rel"], new_rel))
        # internal / drug / session
        internal = jp.norm_internal(a.internal) if a.internal else (jp.norm_internal(jint) if jint else None)
        isrc = "command line (--internal)" if a.internal else ("lab journal header" if jint else "not given")
        drug = a.drug if a.drug is not None else jdrug
        session = a.session or (cfg["session_by_journal_folder"].get(jd.folder) if jd else None)
        if not session:
            js_all = journals_for(date)
            session = cfg["session_by_journal_folder"].get(js_all[0].folder) if js_all else None
        if not session:
            session = "NEW (" + " / ".join(Path(r["rel"]).parts[:2]) + ")"
        gap = h["ptype"] == "sPSC gap-free"
        use = "Y" if (gap and h["dur"] >= min_gf and h["mode"] == "VC") else "N"
        if gap and use == "N":
            notes.append(f"short gap-free recording ({h['dur']} s < {min_gf:g} s): aborted attempt, not analysed")
        row = dict(new_name=new_name, old_name=r["name"], old_folder=rel_dir, new_folder=rel_dir, path_on_mac=new_rel, status=status,
                   date_header=date, start_time=h["time"], session=session, **{"S.No": sno, "T.Code": tcode}, sex=sex, sex_source=sex_src,
                   slice=r["slice"], cell=r["cell"], counter=ctr, protocol=h["protocol"], protocol_type=h["ptype"], mode=h["mode"],
                   duration_s=h["dur"], sweeps=h["sweeps"], psc_content=psc_content_for(internal, drug) if gap else None,
                   drug=(drug or None) if gap else None, internal=internal, internal_source=isrc,
                   holding_header_mV=h["hold_header"], use_for_miniML=use,
                   miniML_directions=(cfg["directions_by_internal"].get(internal) or None) if (use == "Y" and internal) else None)
        if not internal:
            notes.append("internal solution not given (add 'Internal: ...' to the journal or fill the column); step 2 waits for it")
        a_ = assign.get(("new", r["rel"]))
        if gap:
            if a_ and a_["hold"] is not None:
                row["holding_command_mV"] = a_["hold"]; row["holding_listed_in_note"] = a_["hold"]; row["holding_flag"] = a_["flag"]
                if h["vm"] is not None and abs(h["vm"] - a_["hold"]) > 3:
                    row["holding_flag"] += f"; CHECK: the recorded voltage monitor reads {h['vm']:.1f} mV"
            else:
                row["holding_flag"] = (a_["flag"] if a_ else "NOT SET (fill in by hand): short or not matched") if use == "Y" else None
                if h["vm"] is not None and use == "Y":
                    row["holding_flag"] = (row["holding_flag"] or "") + f"; the voltage monitor reads {h['vm']:.1f} mV (not used)"
        row["analysis_note"] = "; ".join(notes) or None
        r["row"] = row
        rows.append(row)

    newf = pd.DataFrame(rows, columns=FILE_COLS) if rows else pd.DataFrame(columns=FILE_COLS)

    # ---- 5. fill blank holdings of existing rows ---------------------------------------------------
    filled_old = []
    for (kind, idx), v in assign.items():
        if kind == "old":
            if v["hold"] is not None:
                files.at[idx, "holding_command_mV"] = v["hold"]; files.at[idx, "holding_listed_in_note"] = v["hold"]
                files.at[idx, "holding_flag"] = v["flag"] + "; filled by step 1 on re-run"
                filled_old.append(files.at[idx, "new_name"])
            else:
                files.at[idx, "holding_flag"] = v["flag"]

    # ---- 6. summary before anything is changed -------------------------------------------------------
    n_ren, n_keep = len(plan), len(newf) - len(plan)
    log(f"\n  PLAN: {len(newf)} new rows ({n_ren} files to rename, {n_keep} keep their names), {len(filled_old)} older rows get a holding from the journal")
    if len(newf):
        gf = newf[newf.protocol_type == "sPSC gap-free"]
        log(f"        gap-free: {len(gf)} | analysable by step 2 (Y): {(newf.use_for_miniML == 'Y').sum()} | holding set: {gf.holding_command_mV.notna().sum()} | "
            f"holding BLANK: {(gf.use_for_miniML.eq('Y') & gf.holding_command_mV.isna()).sum()}")
    for (date, sl, ce, n, st, why) in cell_report:
        log(f"    {date} slice {sl} cell {ce}: {n} gap-free  holding {st}: {why}")
    flagged = newf[newf.analysis_note.notna()]
    grouped = {}
    for _, r in flagged.iterrows():
        log.file_only(f"    FLAG {r['new_name']}: {r['analysis_note']}")
        for part in str(r["analysis_note"]).split("; "):
            grouped.setdefault(re.sub(r"\b\d{3,8}\b", "#", part), []).append(r["new_name"])
    if grouped:
        log("\n  FLAGS (every file is listed in the log file):")
        for k, v in sorted(grouped.items(), key=lambda kv: -len(kv[1])):
            sample = ", ".join(v[:2]) + (" ..." if len(v) > 2 else "")
            log(f"    [{len(v)} files] {k}   e.g. {sample}")
    plan_csv = root / "logs" / f"rename_plan_{stamp}.csv"
    plan_csv.parent.mkdir(exist_ok=True)
    pd.DataFrame([(o, n) for _, _, o, n in plan], columns=["old_path", "new_path"]).to_csv(plan_csv, index=False)
    if a.dry_run:
        log(f"\n  DRY RUN finished. Rename plan: {plan_csv}\n  Run again without --dry-run to apply.")
        return 0
    if not len(newf) and not moved and not filled_old:
        log("\n  Nothing to do.")
        return 0

    # ---- 7. apply: backup first, then rename on disk (with rollback), then write the master ----------------
    bk = pc.backup_master(mp, cfg, "step1", root)
    log(f"\n  backup of the master (made before anything is changed): {bk}")
    for src, dst, _, _ in plan:
        if dst.exists():
            raise pc.UserProblem(f"Cannot rename {src.name}: {dst.name} already exists in that folder.", "Nothing was renamed. Move or rename one of the two files.")
    done = []
    try:
        for src, dst, o, n in plan:
            src.rename(dst); done.append((src, dst))
        ff = em.merge_frames(files.drop(columns=[c for c in files.columns if c.startswith("_") and c not in FILE_COLS]), newf, em.KEYS["files"])
        # cells: recompute rows of every touched cell (keeps old flags)
        ff["_date"] = pd.to_numeric(ff["date_header"], errors="coerce"); ff["_slice"] = pd.to_numeric(ff["slice"], errors="coerce"); ff["_cell"] = pd.to_numeric(ff["cell"], errors="coerce")
        touched = {(r["date"], r["slice"], r["cell"]) for r in recs if r["slice"] is not None} | {(d, s, c) for d, g in day_groups.items() for s, c in g}
        crow = []
        for (d, s, c) in sorted(touched, key=lambda t: tuple(-1 if v is None else v for v in t)):
            g = ff[(ff["_date"] == d) & (ff["_slice"] == s) & (ff["_cell"] == c)].sort_values("start_time")
            if not len(g):
                continue
            gg = g[g.protocol_type == "sPSC gap-free"]
            sess = g["session"].iloc[0]
            old = cells_old[(pd.to_numeric(cells_old.get("date"), errors="coerce") == d) & (cells_old.get("session") == sess) &
                            (pd.to_numeric(cells_old.get("slice"), errors="coerce") == s) & (pd.to_numeric(cells_old.get("cell_new"), errors="coerce") == c)] if len(cells_old) else cells_old
            oldflags = old["flags"].iloc[0] if len(old) and pd.notna(old["flags"].iloc[0]) else None
            _, folder_cell = old_ids(str(g["path_on_mac"].iloc[0]), str(g["old_name"].iloc[0]))
            newflags = [t for t in cell_flags.get((d, s, c), []) if not (oldflags and t in str(oldflags))]
            fl = "; ".join(([str(oldflags)] if oldflags else []) + newflags) or None
            crow.append(dict(date=d, session=sess, slice=s, cell_new=c,
                             cell_old_folder=(old["cell_old_folder"].iloc[0] if len(old) else folder_cell),
                             **{"T.Code": g["T.Code"].iloc[0]}, sex=g["sex"].iloc[0], internal=g["internal"].dropna().iloc[0] if g["internal"].notna().any() else None,
                             internal_source=g["internal_source"].iloc[0], n_files=len(g), n_sPSC_gapfree=len(gg),
                             **{"holding_commands (gap-free, time order)": ", ".join(f"{v:g}" for v in gg["holding_command_mV"].dropna()) or None},
                             flags=fl))
        cf = em.merge_frames(cells_old, pd.DataFrame(crow, columns=CELL_COLS), em.KEYS["cells"]) if crow else cells_old
        rl = pd.concat([rlog, pd.DataFrame([(o, n) for _, _, o, n in plan], columns=["old_path", "new_path"])], ignore_index=True) if plan else rlog
        ff = ff.drop(columns=["_date", "_slice", "_cell"])
        em.write_sheets(mp, {"files": ff, "cells": cf, "rename_log": rl}, backup=False)
    except BaseException:
        for src, dst in reversed(done):          # undo the renames so disk and master agree again
            try:
                dst.rename(src)
            except Exception:
                pass
        raise
    ren_csv = root / "logs" / f"renamed_{stamp}.csv"
    pd.DataFrame([(o, n) for _, _, o, n in plan], columns=["old_path", "new_path"]).to_csv(ren_csv, index=False)
    log(f"\n  DONE. Renamed {len(plan)} files, added {len(newf)} rows to files, updated {len(crow)} cells, rename_log +{len(plan)}.")
    log(f"  Rename list (to undo: --undo-renames \"{ren_csv}\"): {ren_csv}")
    blanks = newf[(newf.use_for_miniML == "Y") & newf.holding_command_mV.isna()]
    if len(blanks):
        log(f"  TO DO FOR YOU: {len(blanks)} gap-free recordings have no holding; type it into files.holding_command_mV (see holding_flag), then run step 2.")
    noint = newf[(newf.use_for_miniML == "Y") & newf.internal.isna()]
    if len(noint):
        log(f"  TO DO FOR YOU: {len(noint)} recordings have no internal solution; fill files.internal.")
    return 0


def undo(a, log):
    cfg = pc.load_config(a.config)
    root, _ = pc.resolve_paths(cfg, a.master, a.data_root)
    t = pd.read_csv(a.undo_renames)
    n = 0
    for _, r in t.iterrows():
        src, dst = root / r["new_path"], root / r["old_path"]
        if src.exists() and not dst.exists():
            src.rename(dst); n += 1
        else:
            log(f"  skipped {r['new_path']} (missing, or the old name exists)")
    log(f"  {n} files renamed back. Now restore the master from _backups/pipeline/ (the copy 'before_step1').")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", type=Path)
    ap.add_argument("--data-root", type=Path)
    ap.add_argument("--config", type=Path)
    ap.add_argument("--journal-dir", type=Path)
    ap.add_argument("--folder", help="only scan this folder (relative to the Ephys folder)")
    ap.add_argument("--internal"); ap.add_argument("--session")
    ap.add_argument("--drug", default=None)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--no-refill-blank", dest="refill_blank", action="store_false", help="do not try to fill blank holdings of older rows")
    ap.add_argument("--undo-renames", type=Path)
    a = ap.parse_args()
    cfg = pc.load_config(a.config)
    root, _ = pc.resolve_paths(cfg, a.master, a.data_root)
    log = pc.new_log(root, "step1")
    try:
        sys.exit(undo(a, log) if a.undo_renames else run(a, log))
    except SystemExit:
        raise
    except BaseException as e:
        sys.exit(pc.report_error(log, "step 1 (ingest)", e))


if __name__ == "__main__":
    main()
