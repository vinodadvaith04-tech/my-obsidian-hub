"""
one_off/rename_adv3_to_adv6.py - 19 Nov 2026: the 37 main-experiment recordings (18 Nov, slice 6, T.Code 1B12) were saved as
Adv3 but the lab journal gives serial number 6 -> renamed on disk to Adv6; this script makes the master follow.
Edits the workbook IN PLACE (only the affected cells): new_name / recording_name / path in files + result sheets, S.No, status,
analysis_note, README text, 37 rows added to rename_log.  Usage: python rename_adv3_to_adv6.py <master.xlsx> <out.xlsx> <renamed_csv>
"""
import sys, re, csv
from openpyxl import load_workbook
src, out, csvp = sys.argv[1:4]
OLD, NEW = "_Adv3_1B12_", "_Adv6_1B12_"
wb = load_workbook(src)
changed = {}
def cols(ws): return {c.value: c.column for c in ws[1] if c.value}
NOTE_RE = re.compile(r"^CHECK NAME: S\.No in the name is Adv3 but the lab journal says serial number 6 for slice 6;?\s*")
for sh, names in {"files": ["new_name", "path_on_mac"], "miniML_per_file": ["new_name", "recording_name"],
                  "miniML_events": ["new_name", "recording_name"], "miniML_per_cell_holding": ["recording_name"],
                  "alt_det_per_file": ["new_name"], "alt_det_events": ["new_name"]}.items():
    ws = wb[sh]; c = cols(ws); n = 0
    for col in names:
        j = c[col]
        for row in ws.iter_rows(min_row=2, min_col=j, max_col=j):
            cell = row[0]
            if isinstance(cell.value, str) and OLD in cell.value:
                cell.value = cell.value.replace(OLD, NEW); n += 1
    if sh == "files":
        for r in range(2, ws.max_row + 1):
            nn = ws.cell(r, c["new_name"]).value
            if isinstance(nn, str) and NEW in nn:
                ws.cell(r, c["S.No"]).value = "Adv6"
                ws.cell(r, c["status"]).value = "renamed on disk 19 Nov (Adv3 -> Adv6, journal serial number 6)"
                note = ws.cell(r, c["analysis_note"]).value or ""
                note = NOTE_RE.sub("", str(note)).strip()
                ws.cell(r, c["analysis_note"]).value = ("name corrected Adv3 -> Adv6 (journal serial 6)" + ("; " + note if note else ""))
    changed[sh] = n
# README
ws = wb["README"]
for row in ws.iter_rows():
    for cell in row:
        v = cell.value
        if not isinstance(v, str): continue
        if "the 37 main-experiment files keep their names" in v:
            cell.value = v.replace("the 37 main-experiment files keep their names (they say Adv3 but the journal gives serial number 6 for slice 6: check)",
                                   "the 37 main-experiment files were renamed Adv3 -> Adv6 on 19 Nov (journal serial number 6 for slice 6)")
        if v.startswith("Old path -> new path of the 995 files"):
            cell.value = ("Old path -> new path of the 1032 files renamed on disk (864 up to 17 Nov, +131 Cs slice 5/6 files renamed 18 Nov, +37 main-experiment files Adv3 -> Adv6 on 19 Nov). "
                          "The 86 slice-4 files and 9 ignored files were never renamed on the Mac.")
# rename_log
ws = wb["rename_log"]
rows = list(csv.reader(open(csvp)))[1:]
for o, n in rows: ws.append([o, n])
wb.save(out)
print("changed cells:", changed, "| rename_log +", len(rows))
