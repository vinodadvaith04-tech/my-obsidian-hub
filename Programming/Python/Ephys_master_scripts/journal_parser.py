"""
journal_parser.py - read the Obsidian lab journal of one recording day.

Understands the layouts used so far (Cs pilot and Main experiment journals):

    **free text at the top**                      (ignored)
    Internal: Cs-gluconate                        (optional day-level lines; also  Drug: gabazine)
    ### Slice 1                                   (Cs pilot style)   or, Main-experiment style, a header block:
    Vgat-cre / DoB: .. / Sex: Male / DoS: .. / Slicer: .. / Plate: 1 / DIV 15 / Animal: 1 / Slice: 6
    Treatment code: 1B1&2 / Serial number: 6
    #### Cell 1
    | Manual check 1 | pipette | Cm | Rm | Ra |     (parsed, not used for the master yet)
    Notes: free text
    cap check
    -80
    -80                                           one line per 2-min gap-free file, in recording order
    cap check
    0 - lots of EPSCs                             (a number followed by a comment is fine; "0mV", "-60 mV" too)

What comes out: JournalDay(date, path, folder, day_header, slices=[JSlice(number, header, cells=[JCell(number, lines, groups, notes, qc)])]).
A holding that cannot be read is simply not there; the caller compares the number of lines with the number of recordings
and leaves the holding blank when they differ.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

KEYS = {
    "sex": "sex", "dob": "dob", "dos": "dos", "slicer": "slicer", "plate": "plate", "div": "div",
    "animal": "animal", "treatment code": "tcode", "t.code": "tcode", "t code": "tcode", "serial number": "serial",
    "s.no": "serial", "serial": "serial", "internal": "internal", "drug": "drug", "drugs": "drug", "blocker": "drug",
    "blockers": "drug", "genotype": "genotype",
}
KEY_RE = re.compile(r"^\s*(sex|dob|dos|slicer|plate|div|animal|treatment code|t\.code|t code|serial number|s\.no|serial|internal|drugs?|blockers?|genotype)\s*[:=]?\s*(.*?)\s*$", re.I)
SLICE_HEAD = re.compile(r"^\s*#{1,6}\s*Slice\s*(\d+)\b", re.I)
SLICE_KEY = re.compile(r"^\s*Slice\s*[:=]\s*(\d+)\b", re.I)
CELL_HEAD = re.compile(r"^\s*#{1,6}\s*Cell\s*(\d+)\b", re.I)
HOLD_RE = re.compile(r"^\s*(-?\d+(?:\.\d+)?)\s*(?:mV)?\b(.*)$", re.I)
QC_RE = re.compile(r"^\s*\|\s*Manual check\s*(\d+)\s*\|(.*)\|\s*$", re.I)
INTERNAL_NAMES = {"kgluconate": "KGlu", "k-gluconate": "KGlu", "kglu": "KGlu", "k gluconate": "KGlu", "k-glu": "KGlu",
                  "csgluconate": "Cs-gluconate", "cs-gluconate": "Cs-gluconate", "cs gluconate": "Cs-gluconate",
                  "cesium": "Cs-gluconate", "caesium": "Cs-gluconate", "cs": "Cs-gluconate",
                  "kcl": "KCl", "k-cl": "KCl"}


def norm_internal(s: str) -> str:
    k = s.strip().lower().rstrip(".")
    return INTERNAL_NAMES.get(k, s.strip())


@dataclass
class JCell:
    number: int
    lines: list = field(default_factory=list)       # one holding (mV) per listed gap-free file, in order
    groups: list = field(default_factory=list)      # holdings split at every "cap check"
    notes: str = ""
    qc: list = field(default_factory=list)
    suspicious: bool = False                         # True = the 'holdings' look like a numbered list or are out of range


@dataclass
class JSlice:
    number: int | None
    header: dict = field(default_factory=dict)
    cells: list = field(default_factory=list)


@dataclass
class JournalDay:
    date: str
    path: Path
    folder: str
    day_header: dict = field(default_factory=dict)
    slices: list = field(default_factory=list)

    def slice_header(self, js: "JSlice", key: str):
        return js.header.get(key) or self.day_header.get(key)


def _num(s):
    s = (s or "").strip().replace(",", ".")
    try:
        return float(s)
    except ValueError:
        return None


def parse_journal(path: Path, journal_root: Path | None = None) -> JournalDay:
    path = Path(path)
    m = re.search(r"(\d{4})-(\d{2})-(\d{2})", path.name)
    date = "".join(m.groups()) if m else ""
    folder = path.parent.name
    day = JournalDay(date=date, path=path, folder=folder)
    pending: dict = {}
    cur_slice: JSlice | None = None
    cur_cell: JCell | None = None
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.rstrip()
        s = line.strip()
        if not s:
            continue
        m = SLICE_HEAD.match(line)
        m2 = SLICE_KEY.match(line)
        if m or m2:
            n = int((m or m2).group(1))
            cur_slice = JSlice(number=n, header=dict(pending))
            pending = {}
            day.slices.append(cur_slice)
            cur_cell = None
            continue
        m = CELL_HEAD.match(line)
        if m:
            if cur_slice is None:                       # cells before any slice heading: one unnamed slice
                cur_slice = JSlice(number=None, header=dict(pending)); pending = {}
                day.slices.append(cur_slice)
            cur_cell = JCell(number=int(m.group(1)))
            cur_slice.cells.append(cur_cell)
            continue
        qm = QC_RE.match(line)
        if qm and cur_cell is not None:
            parts = [p.strip() for p in qm.group(2).split("|")]
            cur_cell.qc.append(dict(check=int(qm.group(1)), pip=_num(parts[0]) if len(parts) > 0 else None,
                                    Cm=_num(parts[1]) if len(parts) > 1 else None, Rm=_num(parts[2]) if len(parts) > 2 else None,
                                    Ra=_num(parts[3]) if len(parts) > 3 else None))
            continue
        if s.startswith("|"):
            continue
        km = KEY_RE.match(line)
        # a "key: value" header line; with no colon only accept the unambiguous ones (e.g. "DIV 15", "Plate 1")
        if km and (":" in line or "=" in line or km.group(1).lower() in ("div", "plate", "animal")):
            key = KEYS.get(km.group(1).lower(), KEYS.get(km.group(1), km.group(1).lower()))
            val = km.group(2).strip()
            if key == "internal":
                val = norm_internal(val)
            if key == "sex":
                val = {"male": "male", "m": "male", "female": "female", "f": "female"}.get(val.lower(), val.lower())
            if key in ("internal", "drug") and cur_slice is None and not pending:
                day.day_header[key] = val              # top-of-day line: applies to every slice
            elif key in ("internal", "drug") and cur_slice is not None and cur_cell is None:
                cur_slice.header[key] = val
            elif cur_slice is not None and cur_cell is None and not cur_slice.cells:
                cur_slice.header[key] = val
            else:
                pending[key] = val                      # a new header block after earlier cells
            continue
        if re.match(r"^\s*cap\s*check", s, re.I):
            if cur_cell is not None:
                cur_cell.groups.append([])
            continue
        if re.match(r"^\s*notes?\s*:", s, re.I):
            if cur_cell is not None:
                cur_cell.notes = (cur_cell.notes + " " + s.split(":", 1)[1].strip()).strip()
            continue
        hm = HOLD_RE.match(line)
        if hm and cur_cell is not None:
            v = float(hm.group(1))
            if -150 <= v <= 100:
                cur_cell.lines.append(v)
                if not cur_cell.groups:
                    cur_cell.groups.append([])
                cur_cell.groups[-1].append(v)
            continue
    for sl in day.slices:                                   # sanity: a list that is not a list of holdings is not used
        for c in sl.cells:
            L = c.lines
            bad = any(not (-120 <= v <= 60) for v in L) or (len(L) >= 3 and L == [float(i) for i in range(1, len(L) + 1)])
            c.suspicious = bool(bad)
    return day


def find_journals(journal_dir: Path, date: str) -> list[Path]:
    """All journal .md files of that day (date = 'YYYYMMDD')."""
    d = f"{date[:4]}-{date[4:6]}-{date[6:]}"
    jd = Path(journal_dir).expanduser()
    if not jd.exists():
        return []
    ephys = jd / "Ephys"
    base = ephys if ephys.exists() else jd
    return sorted(p for p in base.rglob("*.md") if d in p.name)
