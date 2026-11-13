"""
parse_qc.py - read the manual patch checks (Manual check N | pip | Cm | Rm | Ra) from the Lab Journal notes
and add them to the master workbook, sheet manual_qc (rows with the same day/cell/check are replaced).

Usage: python scripts/parse_qc.py --journal "<.../Lab Journal/Ephys/3. Cs internal pilot>" [--session "Cs internal pilot"] [--master ...]
The journal files are only read. Notes are expected as 2026-11-DD*.md with headings "Slice N" / "Cell N".
"""
import argparse, os, re, sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ephys_master as em

ap = argparse.ArgumentParser()
ap.add_argument("--journal", type=Path, required=True)
ap.add_argument("--session", default="Cs internal pilot")
ap.add_argument("--master", type=Path, default=em.default_master())
a = ap.parse_args()
num = lambda s: float(s) if re.fullmatch(r"[\d.]+", s or "") else None
rows = []
for f in sorted(a.journal.glob("20*-*-*.md")):
    day = f.name[:10]; cell = None
    for line in open(f, encoding="utf-8"):
        m = re.match(r"#+\s*Cell\s*(\d+)", line)
        if m: cell = int(m.group(1)); continue
        m = re.match(r"\|\s*Manual check (\d+)\s*\|(.*)\|\s*$", line)
        if m and cell:
            c = [x.strip() for x in m.group(2).split("|")]
            rows.append(dict(session=a.session, day=day, cell=cell, check=int(m.group(1)), pip=num(c[0]), Cm=num(c[1]), Rm=num(c[2]), Ra=num(c[3])))
q = pd.DataFrame(rows).dropna(subset=["Cm", "Rm", "Ra"], how="all")
print(f"{len(q)} manual checks found in {a.journal}")
old = em.read_sheet(a.master, "manual_qc")
out = em.merge_frames(old, q, em.KEYS["manual_qc"], sort=["day", "cell", "check"])
em.write_sheets(a.master, {"manual_qc": out})
print(f"manual_qc: {len(old)} -> {len(out)} rows")
