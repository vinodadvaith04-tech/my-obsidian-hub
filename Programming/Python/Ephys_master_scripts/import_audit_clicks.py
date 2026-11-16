"""
import_audit_clicks.py - copy the clicks you made in audit_tool.html into audit_labels.xlsx.

1. In audit_tool.html press "Export clicks (JSON)"  ->  audit_clicks_YYYYMMDD_HHMM.json   (save it in the audit folder)
2. python scripts/import_audit_clicks.py                    (takes the newest audit_clicks_*.json in the audit folder)
   python scripts/import_audit_clicks.py path/to/file.json [more.json ...] [--audit-dir <audit folder>] [--dry-run] [--overwrite]

What it does (audit_labels.xlsx must be closed in Excel; a backup goes to <master folder>/_backups/ first):
  * only windows marked Done (or Skip) in the tool are written; windows still open are left empty.
  * fills the same pink columns score_audit.py reads (miniML_false_ids, miniML_missed, miniML_missed_up, old_*, missed_times_s,
    unsure_ids, note), so   python scripts/score_audit.py   works unchanged.
  * never overwrites something you typed by hand: a cell is only replaced if it is empty or still holds what an earlier import wrote
    (logged in the sheet audit_tool_import). Windows with a conflict are reported and skipped unless you pass --overwrite.
  * adds the sheet audit_truth_events: every real event in the finished windows (markers you left as real + events you clicked),
    with which detector found it. The clicked events (with their exact times) are the new training positives for the model.

How the counts are defined (same rules as the tool shows live):
  real event   = a marker you did not flag false/unsure, from either detector, plus every event you clicked.
                 Markers of the two detectors closer than 10 ms are one event.
  missed by D  = a real event with no real marker of detector D (and no unsure marker of D within 10 ms).
  *_missed_up  = how many of those are upward.   missed_times_s = times (s) of the events miniML missed.
Needs: pandas openpyxl
"""
import argparse, json, shutil, sys, time
from pathlib import Path
import pandas as pd
from openpyxl import load_workbook

TOL = 0.010
LABEL_COLS = ["miniML_false_ids", "miniML_missed", "miniML_missed_up", "old_false_ids", "old_missed", "old_missed_up", "missed_times_s", "unsure_ids", "note"]
HOW_TO_ADD = [
    "INTERACTIVE ALTERNATIVE: open audit_tool.html (same windows), click false markers / missed events, press Export clicks, then run scripts/import_audit_clicks.py - it fills these pink columns for you.",
    "Rows filled by the importer are logged in the sheet audit_tool_import; anything you type by hand is never overwritten by the importer.",
]


def norm(v):
    if v is None: return None
    if isinstance(v, float) and v != v: return None
    if isinstance(v, float) and v.is_integer(): v = int(v)
    s = str(v).strip()
    return s if s else None


def events_for_window(markers, state):
    """markers: list of dict(id, det, t, pol, amp) ordered M first (time order), then O. state: false/unsure id lists, missed [{t,pol,y}]."""
    fl, un = set(state.get("false", [])), set(state.get("unsure", []))
    real = lambda d: [m for m in markers if m["det"] == d and m["id"] not in fl and m["id"] not in un]
    ev = [dict(t=m["t"], pol=m["pol"], amp=m["amp"], src={"M": m["id"]}) for m in real("M")]
    for m in sorted(real("O"), key=lambda m: m["t"]):
        best, bd = None, 1e9
        for e in ev:
            if "O" in e["src"] or "M" not in e["src"]: continue
            d = abs(e["t"] - m["t"])
            if d <= TOL and d < bd: bd, best = d, e
        if best is not None: best["src"]["O"] = m["id"]
        else: ev.append(dict(t=m["t"], pol=m["pol"], amp=m["amp"], src={"O": m["id"]}))
    for c in state.get("missed", []):
        best, bd = None, 1e9
        for e in ev:
            d = abs(e["t"] - c["t"])
            if d <= TOL and d < bd: bd, best = d, e
        if best is not None: best["src"]["C"] = 1
        else: ev.append(dict(t=c["t"], pol=c["pol"], amp=c.get("y"), src={"C": 1}, click_amp=True))
    ev.sort(key=lambda e: e["t"])
    res = {"events": ev}
    for d in ("M", "O"):
        um = [m for m in markers if m["det"] == d and m["id"] in un]
        miss = []
        for e in ev:
            if d in e["src"]: continue
            if any(abs(m["t"] - e["t"]) <= TOL for m in um):
                e["unc" + d] = True; continue
            e["miss" + d] = True; miss.append(e)
        res[d] = dict(missed=len(miss), up=sum(1 for e in miss if e["pol"] == "up"), times=[e["t"] for e in miss],
                      n=sum(1 for m in markers if m["det"] == d), false=sum(1 for m in markers if m["det"] == d and m["id"] in fl),
                      unsure=len(um))
    return res


def num_ids(ids, det):
    return sorted(int(i[1:]) for i in ids if i[:1] == det and i[1:].isdigit())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("json", nargs="*", type=Path)
    ap.add_argument("--audit-dir", type=Path, default=Path(__file__).resolve().parent.parent / "audit")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--overwrite", action="store_true", help="also replace cells you typed by hand")
    a = ap.parse_args()
    xl = a.audit_dir / "audit_labels.xlsx"
    files = a.json or sorted(a.audit_dir.glob("audit_clicks_*.json"))[-1:]
    if not files: sys.exit(f"no audit_clicks_*.json found in {a.audit_dir}; export one from audit_tool.html first")
    # merge exports: newest ts per window wins
    wins = {}
    for f in files:
        j = json.loads(Path(f).read_text(encoding="utf-8"))
        if abs(j.get("tol_s", TOL) - TOL) > 1e-9: print(f"WARNING {f}: tolerance {j.get('tol_s')} differs from {TOL}")
        for wid, x in j["windows"].items():
            if wid not in wins or x.get("ts", 0) >= wins[wid].get("ts", 0): wins[wid] = x
        print("read", f, len(j["windows"]), "windows")

    key = pd.read_excel(xl, sheet_name="audit_events_key")
    win = pd.read_excel(xl, sheet_name="audit_windows").set_index("window")
    key["det"] = key["id"].str[0]
    key = key.sort_values(["window", "det", "time_s"], ascending=[True, False, True])
    mk_by = {w: [dict(id=r.id, det=r.det, t=round(float(r.time_s), 5), pol=r.polarity, amp=float(r.amp_pA)) for r in g.itertuples()] for w, g in key.groupby("window")}

    wb = load_workbook(xl)
    ws = wb["audit_labels"]
    hdr = {c.value: c.column for c in ws[1]}
    rowof = {ws.cell(r, hdr["window"]).value: r for r in range(2, ws.max_row + 1)}
    log_name = "audit_tool_import"
    logws = wb[log_name] if log_name in wb.sheetnames else None
    prev = {}
    if logws is not None:
        for r in logws.iter_rows(min_row=2, values_only=True):
            if r and r[0]: prev[r[0]] = json.loads(r[2]) if r[2] else {}

    written, conflicts, skipped_open, mismatch = {}, [], 0, []
    truth_rows, newlog = [], dict(prev)
    for wid, x in sorted(wins.items()):
        if wid not in rowof or wid not in win.index:
            print(f"  {wid}: not in audit_labels.xlsx - ignored"); continue
        if x["file"] != win.loc[wid, "file"] or abs(float(x["start_s"]) - float(win.loc[wid, "start_s"])) > 0.002:
            print(f"  {wid}: file/start differ from the workbook - ignored"); continue
        if x["status"] not in ("done", "skip"):
            skipped_open += 1; continue
        mk = mk_by.get(wid, [])
        if x["status"] == "skip":
            new = {c: None for c in LABEL_COLS}; new["note"] = x.get("note") or "skip"
        else:
            r = events_for_window(mk, x)
            exp = x.get("counts", {})
            for det, nm in (("M", "miniML"), ("O", "old")):
                ec = exp.get(nm, {})
                if ec and (ec.get("missed") != r[det]["missed"] or ec.get("missed_up") != r[det]["up"] or ec.get("false") != r[det]["false"]):
                    mismatch.append(f"{wid} {nm}: tool says {ec}, recomputed missed={r[det]['missed']} up={r[det]['up']} false={r[det]['false']}")
            join = lambda L: ";".join(str(v) for v in L) or None
            new = dict(
                miniML_false_ids=join(num_ids(x["false"], "M")), miniML_missed=r["M"]["missed"], miniML_missed_up=r["M"]["up"] if r["M"]["missed"] else None,
                old_false_ids=join(num_ids(x["false"], "O")), old_missed=r["O"]["missed"], old_missed_up=r["O"]["up"] if r["O"]["missed"] else None,
                missed_times_s=";".join(f"{t:.3f}" for t in r["M"]["times"]) or None,
                unsure_ids=";".join(sorted(x["unsure"], key=lambda s: (s[0], int(s[1:])))) or None,
                note=x.get("note") or None)
            for e in r["events"]:
                truth_rows.append(dict(window=wid, file=x["file"], time_s=round(e["t"], 5), polarity=e["pol"],
                                       found_by_miniML=("M" in e["src"]), found_by_old=("O" in e["src"]), clicked_by_you=("C" in e["src"]),
                                       miniML_marker=e["src"].get("M"), old_marker=e["src"].get("O"),
                                       amp_pA=None if e.get("amp") is None else round(float(e["amp"]), 2),
                                       amp_source="click, low-passed peak" if e.get("click_amp") else ("miniML fit" if "M" in e["src"] else "second detector"),
                                       holding_mV=win.loc[wid, "holding_mV"], internal=win.loc[wid, "internal"], drug=win.loc[wid, "drug"],
                                       batch=win.loc[wid, "batch"]))
        rr = rowof[wid]
        cur = {c: norm(ws.cell(rr, hdr[c]).value) for c in LABEL_COLS}
        newn = {c: norm(new.get(c)) for c in LABEL_COLS}
        if cur == newn:
            newlog[wid] = [wid, None, json.dumps(newn)]; continue
        had = any(v is not None for v in cur.values())
        mine = prev.get(wid)
        if had and not a.overwrite and cur != mine:
            conflicts.append(wid); continue
        if not a.dry_run:
            for c in LABEL_COLS: ws.cell(rr, hdr[c]).value = new.get(c)
        written[wid] = newn
        newlog[wid] = [wid, None, json.dumps(newn)]

    print(f"\n{len(written)} window(s) written, {skipped_open} still open in the tool, {len(conflicts)} conflict(s)")
    if conflicts: print("  conflicts (cells you filled by hand differ from the tool; rerun with --overwrite to replace them):", ", ".join(conflicts))
    if mismatch:
        print("WARNING: the tool and this script disagree on:"); [print("  ", m) for m in mismatch]
    if a.dry_run: print("dry run: nothing saved"); return

    # truth sheet (rebuilt from all windows in these exports, derived data)
    if truth_rows:
        if "audit_truth_events" in wb.sheetnames: del wb["audit_truth_events"]
        t = wb.create_sheet("audit_truth_events")
        df = pd.DataFrame(truth_rows)
        t.append(list(df.columns))
        for row in df.itertuples(index=False): t.append([None if (isinstance(v, float) and v != v) else (bool(v) if isinstance(v, (bool,)) else v) for v in row])
        for i, c in enumerate(df.columns, 1): t.column_dimensions[t.cell(1, i).column_letter].width = max(12, min(40, len(c) + 4))
    # import log
    if log_name in wb.sheetnames: del wb[log_name]
    lg = wb.create_sheet(log_name); lg.append(["window", "imported_at", "written_values_json"])
    stamp = time.strftime("%Y-%m-%d %H:%M")
    for wid, v in sorted(newlog.items()): lg.append([wid, stamp if wid in written else (prev.get(wid) and "earlier") or stamp, v[2] if isinstance(v, list) else json.dumps(v)])
    # how_to lines
    h = wb["how_to"]; have = {str(r[0].value) for r in h.iter_rows(min_row=1) if r and r[0].value}
    for line in HOW_TO_ADD:
        if line not in have: h.append([line])
    # backup, then atomic save
    bdir = a.audit_dir.parent / "_backups"; bdir.mkdir(exist_ok=True)
    shutil.copy2(xl, bdir / f"audit_labels_{time.strftime('%Y%m%d_%H%M%S')}.xlsx")
    tmp = xl.with_name(xl.stem + ".tmp.xlsx"); wb.save(tmp); tmp.replace(xl)
    print("saved", xl, "(backup in", bdir, ")")
    print("next: python scripts/score_audit.py")


if __name__ == "__main__":
    main()
