"""
csv_to_xlsx.py - make a readable Excel copy of every .csv under a folder.

Why: Excel with Dutch/European regional settings expects ';' as the separator, so comma-separated
CSVs open as ONE column. A real .xlsx has no such problem, keeps numbers as numbers, and can be read
and written by scripts (pandas.read_excel / DataFrame.to_excel, openpyxl).

Usage:   python3 csv_to_xlsx.py <folder> [--overwrite]
Rules:
  * the .csv is never changed or deleted
  * the .xlsx gets the same name next to it, one sheet named after the file (max 31 characters)
  * an existing .xlsx with that name is NOT replaced (it may hold your own edits, e.g. audit labels)
    unless --overwrite is given
  * header row: bold, frozen, with filter buttons; column widths fitted; numbers stay numbers
"""
import argparse, sys
from pathlib import Path
import pandas as pd
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

MAX_ROWS = 1_048_575


def convert(csv: Path, overwrite: bool):
    out = csv.with_suffix(".xlsx")
    if out.exists() and not overwrite:
        return "skipped (xlsx already exists)", out
    df = pd.read_csv(csv, low_memory=False)
    if len(df) > MAX_ROWS:
        return f"NOT converted ({len(df)} rows is more than Excel can hold)", out
    sheet = csv.stem[:31]
    with pd.ExcelWriter(out, engine="openpyxl") as xw:
        df.to_excel(xw, sheet_name=sheet, index=False)
        ws = xw.sheets[sheet]
        fill = PatternFill("solid", fgColor="D9D9D9")
        for c in ws[1]:
            c.font = Font(bold=True); c.fill = fill; c.alignment = Alignment(vertical="top", wrap_text=False)
        ws.freeze_panes = "A2"
        if len(df.columns):
            ws.auto_filter.ref = f"A1:{get_column_letter(len(df.columns))}{len(df) + 1}"
        samp = df.head(300)
        for j, col in enumerate(df.columns, 1):
            w = max([len(str(col))] + [len(f"{v:.4g}") if isinstance(v, float) else len(str(v)) for v in samp[col].tolist()])
            ws.column_dimensions[get_column_letter(j)].width = min(max(w + 2, 8), 60)
    return f"ok ({len(df)} rows x {len(df.columns)} columns)", out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", type=Path)
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()
    files = sorted(p for p in a.folder.rglob("*.csv") if not p.name.startswith("._"))
    print(f"{len(files)} csv files under {a.folder}", flush=True)
    for p in files:
        try:
            msg, out = convert(p, a.overwrite)
        except Exception as e:  # keep going, report at the end
            msg, out = f"FAILED: {e}", p
        print(f"{p.relative_to(a.folder)} -> {msg}", flush=True)


if __name__ == "__main__":
    main()
