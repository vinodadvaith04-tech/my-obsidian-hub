"""
ephys_master.py - shared helpers for the ONE master workbook (Ephys_master.xlsx).

Every script reads from and writes to the same workbook, so all data and all analyses live in one
place. One sheet per kind of table; filter the sheets in Excel to choose what to analyse.

Sheets (see the README sheet inside the workbook for the column guide):
  README, runs, files, cells, rename_log, manual_qc,
  miniML_per_file, miniML_events, miniML_per_cell_holding,
  alt_det_per_file, alt_det_events

Rules the helpers follow:
  * results are ADDED (upsert by key); a result is only replaced when the same key is produced again
  * a timestamped backup of the workbook is made before every write (last 5 kept, in _backups/)
  * the workbook is written to a temporary file first and then moved into place
  * you can edit the sheets by hand in Excel (e.g. the holding in `files`); scripts read what is there
"""
from __future__ import annotations

import collections
import os
import shutil
import time
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import Workbook, load_workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

MASTER_NAME = "Ephys_master.xlsx"
SHEET_ORDER = ["README", "runs", "files", "cells", "rename_log", "manual_qc",
               "miniML_per_file", "miniML_events", "miniML_per_cell_holding",
               "alt_det_per_file", "alt_det_events"]
KEYS = {
    "files": ["new_name", "old_name"],                       # new_name is empty for ignored files -> old_name
    "cells": ["date", "session", "slice", "cell_new"],
    "manual_qc": ["day", "cell", "check"],
    "miniML_per_file": ["new_name", "direction", "params_hash"],
    "miniML_events": ["new_name", "direction", "params_hash"],          # many events share one key
    "miniML_per_cell_holding": ["recording_name", "holding_command_mV", "direction", "params_hash"],
    "alt_det_per_file": ["new_name"],
    "alt_det_events": ["new_name"],
}
HEADER_FILL = PatternFill("solid", fgColor="D9D9D9")


# ---------------------------------------------------------------- locating files
def default_master(folder: Path | None = None) -> Path:
    """Ephys_master.xlsx in the Ephys folder (the folder above scripts/ unless given)."""
    base = Path(folder) if folder else Path(__file__).resolve().parent.parent
    return base / MASTER_NAME


def index_abf(root: Path):
    """{file name: path} for every .abf below root, plus names that occur more than once."""
    idx, dups = {}, collections.defaultdict(list)
    for p in Path(root).rglob("*.abf"):
        if p.name.startswith("._") or not p.is_file():   # some "x.abf" are image folders, not recordings
            continue
        if p.name in idx:
            dups[p.name].append(p)
        else:
            idx[p.name] = p
    return idx, dups


def locate(row, root: Path, idx: dict) -> Path | None:
    """Where is this recording on disk? Tries the stored relative path, then new_name, then old_name."""
    rel = row.get("path_on_mac") if hasattr(row, "get") else None
    if isinstance(rel, str) and rel.strip():
        p = Path(root) / rel
        if p.exists():
            return p
    for k in ("new_name", "old_name"):
        n = row.get(k) if hasattr(row, "get") else None
        if isinstance(n, str) and n in idx:
            return idx[n]
    return None


# ---------------------------------------------------------------- reading
def sheet_names(path: Path):
    wb = load_workbook(path, read_only=True)
    try:
        return list(wb.sheetnames)
    finally:
        wb.close()


def read_sheet(path: Path, name: str) -> pd.DataFrame:
    """One sheet as a DataFrame (a missing sheet gives an empty frame)."""
    if name not in sheet_names(path):
        return pd.DataFrame()
    return pd.read_excel(path, sheet_name=name)


# ---------------------------------------------------------------- writing
def _backup(path: Path, keep: int = 5):
    if not path.exists():
        return
    bdir = path.parent / "_backups"
    bdir.mkdir(exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    shutil.copy2(path, bdir / f"{path.stem}_{stamp}.xlsx")
    old = sorted(bdir.glob(f"{path.stem}_*.xlsx"))
    for p in old[:-keep]:
        p.unlink()


def _fmt_sheet(ws, df: pd.DataFrame, wrap_first_col: bool = False):
    for c in ws[1]:
        c.font = Font(bold=True)
        c.fill = HEADER_FILL
    ws.freeze_panes = "A2"
    if len(df.columns) and len(df):
        ws.auto_filter.ref = f"A1:{get_column_letter(len(df.columns))}{len(df) + 1}"
    samp = df.head(300)
    for j, col in enumerate(df.columns, 1):
        vals = [len(str(col))] + [len(f"{v:.4g}") if isinstance(v, float) else len(str(v)) for v in samp[col].tolist()]
        ws.column_dimensions[get_column_letter(j)].width = min(max(max(vals) + 2, 8), 60)


def _clean(df: pd.DataFrame) -> pd.DataFrame:
    """Values openpyxl can write: NaN/NaT -> empty, numpy bool/int -> python."""
    df = df.copy()
    df = df.astype(object).where(df.notna(), None)
    return df


def _put_sheet(wb: Workbook, name: str, df: pd.DataFrame, index: int | None = None):
    if name in wb.sheetnames:
        pos = wb.sheetnames.index(name)
        del wb[name]
    else:
        pos = index if index is not None else len(wb.sheetnames)
    ws = wb.create_sheet(name, pos)
    ws.append(list(df.columns))
    for row in _clean(df).itertuples(index=False, name=None):
        ws.append(list(row))
    _fmt_sheet(ws, df)
    if name == "README":
        ws.column_dimensions["A"].width = 22
        ws.column_dimensions["B"].width = 26
        ws.column_dimensions["C"].width = 130
        for row in ws.iter_rows(min_row=2):
            for c in row:
                c.alignment = Alignment(wrap_text=True, vertical="top")
    return ws


def _order(wb: Workbook):
    names = [s for s in SHEET_ORDER if s in wb.sheetnames] + [s for s in wb.sheetnames if s not in SHEET_ORDER]
    wb._sheets = [wb[n] for n in names]


def write_sheets(path: Path, updates: dict, backup: bool = True):
    """Replace the given sheets (name -> DataFrame) in the workbook; all other sheets stay as they are."""
    path = Path(path)
    if path.exists():
        if backup:
            _backup(path)
        wb = load_workbook(path)
    else:
        wb = Workbook()
        wb.remove(wb.active)
    for name, df in updates.items():
        _put_sheet(wb, name, df)
    _order(wb)
    tmp = path.with_name(path.stem + ".tmp.xlsx")
    wb.save(tmp)
    os.replace(tmp, path)


def _keystr(df: pd.DataFrame, kk: list[str]) -> pd.Series:
    """One text key per row; empty cells are '', 5.0 and 5 give the same text (works on pandas 2 and 3)."""
    def one(v):
        if v is None or (isinstance(v, float) and np.isnan(v)) or v is pd.NA:
            return ""
        if isinstance(v, (float, np.floating)) and float(v).is_integer():
            return str(int(v))
        return str(v)
    return pd.Series(["|".join(one(v) for v in row) for row in df[kk].itertuples(index=False, name=None)], index=df.index)


def merge_frames(old: pd.DataFrame, new: pd.DataFrame, key: list[str], sort: list[str] | None = None) -> pd.DataFrame:
    """Old rows whose key also occurs in `new` are dropped, then the new rows are added (pure function)."""
    if not len(old):
        out = new.reset_index(drop=True)
    else:
        old = old.copy()
        for c in new.columns:
            if c not in old.columns:
                old[c] = np.nan
        new = new.copy()
        for c in old.columns:
            if c not in new.columns:
                new[c] = np.nan
        kk = [c for c in key if c in old.columns and c in new.columns]
        newkeys = set(_keystr(new, kk))
        keep = ~_keystr(old, kk).isin(newkeys)
        out = pd.concat([old[keep], new[old.columns]], ignore_index=True)
    if sort:
        out = out.sort_values([c for c in sort if c in out.columns], kind="stable").reset_index(drop=True)
    return out


def upsert(path: Path, name: str, new: pd.DataFrame, key: list[str] | None = None,
           sort: list[str] | None = None, backup: bool = True) -> pd.DataFrame:
    """Add rows to ONE sheet and save (see merge_frames). Use write_sheets for several sheets in one save."""
    key = key or KEYS[name]
    old = read_sheet(path, name) if Path(path).exists() else pd.DataFrame()
    out = merge_frames(old, new, key, sort)
    write_sheets(path, {name: out}, backup=backup)
    return out


# ---------------------------------------------------------------- shared analysis helpers
def assumed_event_type(direction: str, internal: str) -> str:
    """Working label from polarity: KCl -> downward = IPSC; otherwise downward = EPSC, upward = IPSC.
    This is an ASSUMPTION from the holding design, not a measurement."""
    if str(internal) == "KCl":
        return "IPSC" if direction == "negative" else "?"
    return "EPSC" if direction == "negative" else "IPSC"


def cell_holding_table(per_file: pd.DataFrame, events: pd.DataFrame) -> pd.DataFrame:
    """Per cell x holding x direction x settings: pooled from the per-recording and per-event tables."""
    if not len(per_file):
        return pd.DataFrame()
    keys = ["recording_name", "holding_command_mV", "direction", "params_hash"]
    g = per_file.groupby(keys, dropna=False)
    cell = g.agg(n_files=("new_name", "nunique"), n_events=("n_events", "sum"),
                 valid_time_s=("valid_time_s", "sum"),
                 session=("session", "first"), date_header=("date_header", "first"),
                 sex=("sex", "first"), T_code=("T.Code", "first"), internal=("internal", "first"),
                 psc_content=("psc_content", "first"), holding_flag=("holding_flag", "first")).reset_index()
    cell["frequency_Hz"] = cell["n_events"] / cell["valid_time_s"]
    cell.insert(4, "event_type_assumed", [assumed_event_type(d, i) for d, i in zip(cell["direction"], cell["internal"])])
    if len(events):
        a = (events.groupby(keys, dropna=False)["amplitude_abs"]
             .agg(amplitude_abs_mean="mean", amplitude_abs_median="median").reset_index())
        cell = cell.merge(a, on=keys, how="left")
    else:
        cell["amplitude_abs_mean"] = np.nan
        cell["amplitude_abs_median"] = np.nan
    return cell


# ---------------------------------------------------------------- metadata attached to every result row
META = ["session", "date_header", "T.Code", "sex", "slice", "cell", "internal", "psc_content", "drug",
        "holding_command_mV", "holding_flag", "ljp_in_note", "duration_s", "miniML_directions"]
EV_META = ["recording_name", "holding_command_mV", "holding_flag", "internal", "psc_content", "sex", "T.Code",
           "date_header", "cell"]


def meta_table(files: pd.DataFrame) -> pd.DataFrame:
    """Per-recording metadata (indexed by new_name) taken from the `files` sheet, so result rows are always
    complete and follow any hand edit of the files sheet."""
    fm = files.dropna(subset=["new_name"]).drop_duplicates("new_name").set_index("new_name")
    cols = [c for c in META if c in fm.columns]
    fm = fm[cols].copy()
    stem = fm.index.to_series().str.replace(r"\.abf$", "", regex=True)
    fm.insert(0, "recording_index", stem.str.rsplit("_", n=1).str[1].values)
    fm.insert(0, "recording_name", stem.str.rsplit("_", n=1).str[0].values)
    return fm


def with_meta_per_file(res: pd.DataFrame, fm: pd.DataFrame) -> pd.DataFrame:
    """Order: new_name, recording_name, recording_index, direction, event_type_assumed, meta..., results..."""
    drop = ["recording_name", "recording_index", "event_type_assumed"] + META
    res_cols = [c for c in res.columns if c not in drop + ["new_name", "direction"]]
    out = res[["new_name", "direction"] + res_cols].merge(fm, left_on="new_name", right_index=True, how="left")
    out.insert(4, "event_type_assumed", [assumed_event_type(d, i) for d, i in zip(out.direction, out.internal)])
    lead = ["new_name", "recording_name", "recording_index", "direction", "event_type_assumed"] + \
           [c for c in META if c in out.columns]
    return out[lead + [c for c in out.columns if c not in lead]]


def with_meta_events(ev: pd.DataFrame, fm: pd.DataFrame) -> pd.DataFrame:
    keep = [c for c in ev.columns if c not in EV_META + ["event_type_assumed"]]
    out = ev[keep].merge(fm[[c for c in EV_META if c in fm.columns]], left_on="new_name", right_index=True, how="left")
    out.insert(2, "event_type_assumed", [assumed_event_type(d, i) for d, i in zip(out.direction, out.internal)])
    return out[[c for c in out.columns if c not in EV_META] + [c for c in EV_META if c in out.columns]]
