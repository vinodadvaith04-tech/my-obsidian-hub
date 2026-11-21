"""
pick_audit.py - ONE audit-window picker for every round (replaces pick_audit_round2/3/4.py).
You describe WHAT you want with --rule lines; the script finds matching 3 s windows, never reuses an audited window,
spreads the picks over cells, and writes the csv that build_audit.py reads. With --build it also adds the windows to
audit_labels.xlsx and rebuilds audit_tool.html in one go.

Run from the Ephys folder (the one with Ephys_master.xlsx):

  python scripts/pick_audit.py --list                          what is available (windows per internal / holding / day)
  python scripts/pick_audit.py --name round5 --n-windows 80    default mix (see below), csv only
  python scripts/pick_audit.py --name round5 --build \
      --rule "internal=Cs-gluconate holding=-40,-30 kind=missed n=10 tag=Cs weak spot" \
      --rule "internal=KGlu kind=fp n=12" \
      --rule "internal=Cs-gluconate,KGlu,KCl kind=random n=20"

A rule = words key=value separated by spaces (several values with commas):
  kind=      random   random windows (the unbiased test set; the scorer treats them as unbiased)
             missed   windows where the second detector marks events that miniML does not (ranked by how many)
             fp       windows where miniML marks low-score events (<0.8) that the second detector does not (false-event hunt)
             active   windows with the most events (both detectors together)
  n=         how many windows this rule should give
  internal=  Cs-gluconate | KGlu | KCl           (commas for several)
  holding=   command holding in mV, e.g. -40 or -50,-40,-30
  day=       recording date(s), e.g. 20261118 or 20261117,20261118      (start of the file name)
  cell=      text that must be in the cell key, e.g. slice6_C2
  min_events=  at least this many detections (both detectors together) in the window
  tag=       label shown in the audit sheet (stratum); default: kind, internal, holding, day
Rules are used in the order given; windows taken by an earlier rule are not taken again.
Without any --rule a default mix is used: 40 % random (Cs, KGlu, KCl), 30 % Cs missed, 15 % KGlu fp, 15 % Cs active.
Rules can also be kept in a text file (one per line, # = comment):  --rules-file audit/my_rules.txt

Options: --data-root (folder with the .abf files, default = the master's folder), --n-windows (only for the default mix), --seed, --min-duration (s, default 30), --max-per-recording (default 2),
  --no-fresh-preference (by default cells that have no audited window yet are preferred), --name, --batch, --out,
  --master, --audit-dir, --build, --dry-run (with --build: only print what would be run)
Needs: pandas numpy openpyxl (+ pyabf scipy matplotlib for --build)
"""
import argparse, shlex, subprocess, sys, time, shutil
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ephys_master as em

WIN, TOL, STEP = 3.0, 0.010, 1.0
WHY = {"random": "random", "missed": "disagreement", "fp": "disagreement", "active": "disagreement"}   # labels the tool/scoring know
KINDS = list(WHY)


def parse_rule(text):
    r = dict(kind="random", n=10, internal=None, holding=None, day=None, cell=None, min_events=0, tag=None)
    toks = shlex.split(text)
    # words without "=" belong to the previous value (e.g. tag=Cs weak spot)
    pairs, last = [], None
    for t in toks:
        if "=" in t and t.split("=", 1)[0] in r:
            k, v = t.split("=", 1); pairs.append([k, v]); last = pairs[-1]
        elif last is not None and last[0] == "tag":
            last[1] += " " + t
        else:
            sys.exit(f"cannot read rule '{text}': '{t}' (use key=value; keys: {', '.join(r)})")
    for k, v in pairs:
        if k == "n" or k == "min_events": r[k] = int(v)
        elif k == "holding": r[k] = [float(x) for x in v.split(",")]
        elif k in ("internal", "day"): r[k] = [x.strip() for x in v.split(",")]
        else: r[k] = v
    if r["kind"] not in KINDS: sys.exit(f"kind must be one of {KINDS} (got {r['kind']})")
    r["tag"] = r["tag"] or " ".join(x for x in [r["kind"], ",".join(r["internal"] or []), ("hold " + ",".join(f"{x:g}" for x in r["holding"])) if r["holding"] else "", ("day " + ",".join(r["day"])) if r["day"] else ""] if x)
    return r


def default_rules(n):
    k = lambda f: max(1, int(round(n * f)))
    return [parse_rule(s) for s in (
        f"kind=random internal=Cs-gluconate n={k(.22)} tag=Random",
        f"kind=random internal=KGlu n={k(.12)} tag=Random",
        f"kind=random internal=KCl n={k(.06)} tag=Random",
        f"kind=missed internal=Cs-gluconate n={k(.30)} tag=Detectors disagree",
        f"kind=fp internal=KGlu n={k(.15)} tag=K-gluconate FP hunt",
        f"kind=active internal=Cs-gluconate n={k(.15)} min_events=3 tag=Busy Cs windows")]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", type=Path, default=em.default_master())
    ap.add_argument("--audit-dir", type=Path, default=None)
    ap.add_argument("--data-root", type=Path, default=None, help="folder with the .abf files (default: the master's folder)")
    ap.add_argument("--rule", action="append", default=[])
    ap.add_argument("--rules-file", type=Path, default=None)
    ap.add_argument("--n-windows", type=int, default=80)
    ap.add_argument("--seed", type=int, default=5)
    ap.add_argument("--min-duration", type=float, default=30.0)
    ap.add_argument("--max-per-recording", type=int, default=2)
    ap.add_argument("--no-fresh-preference", action="store_true")
    ap.add_argument("--name", default=None, help="e.g. round5 -> audit/audit_round5_windows.csv and batch label")
    ap.add_argument("--batch", default=None, help="label written into the audit sheet (default: --name)")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--build", action="store_true", help="after picking: build_audit.py + build_audit_tool.py")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    audit_dir = a.audit_dir or (a.master.parent / "audit")

    runs = em.read_sheet(a.master, "runs")
    h = runs[runs.tool == "miniML"].sort_values("last_processed").params_hash.iloc[-1]
    pf = em.read_sheet(a.master, "miniML_per_file"); pf = pf[pf.params_hash == h].drop_duplicates("new_name")
    me = em.read_sheet(a.master, "miniML_events"); me = me[me.params_hash == h]
    oe = em.read_sheet(a.master, "alt_det_events")
    pf = pf[pf.new_name.isin(set(em.read_sheet(a.master, "alt_det_per_file").new_name)) & (pf.duration_s >= a.min_duration)].copy()
    pf["day"] = pf.new_name.str[:8]
    pf["allowed"] = [{("down" if d.strip() == "negative" else "up") for d in str(x).split(",")} for x in pf.miniML_directions]
    me = me.assign(pol=np.where(me.direction == "negative", "down", "up"))
    oe = oe.assign(pol=np.where(oe.polarity == "inward", "down", "up"))
    me_by = {k: v for k, v in me.groupby("new_name")}; oe_by = {k: v for k, v in oe.groupby("new_name")}

    done = {}
    xl = audit_dir / "audit_labels.xlsx"
    if xl.exists():
        w = pd.read_excel(xl, sheet_name="audit_windows")
        for f, s in zip(w.file, w.start_s): done.setdefault(f, []).append(float(s))

    def tables(r):
        m = me_by.get(r.new_name); o = oe_by.get(r.new_name)
        m = m[m.pol.isin(r.allowed)] if m is not None else pd.DataFrame(columns=me.columns)
        o = o[o.pol.isin(r.allowed)] if o is not None else pd.DataFrame(columns=oe.columns)
        mt, ot = m.peak_time_s.values, o.peak_s.values
        o_only = np.array([t for t in ot if len(mt) == 0 or np.min(np.abs(mt - t)) > TOL])
        m_only = m[np.array([len(ot) == 0 or np.min(np.abs(ot - t)) > TOL for t in mt], dtype=bool)]
        return mt, ot, o_only, m_only[m_only.score < 0.8].peak_time_s.values

    def cnt(x, s): return int(((x >= s) & (x < s + WIN)).sum()) if len(x) else 0

    cands = []
    for r in pf.itertuples():
        mt, ot, o_only, m_low = tables(r)
        for s in np.arange(0, r.duration_s - WIN, STEP):
            if any(abs(s - t) < WIN for t in done.get(r.new_name, [])): continue
            cands.append(dict(file=r.new_name, start_s=float(s), internal=r.internal, holding=r.holding_command_mV, day=r.day,
                              cell=r.new_name.rsplit("_", 1)[0], n_m=cnt(mt, s), n_o=cnt(ot, s),
                              n_missed=cnt(o_only, s), n_lowfp=cnt(m_low, s)))
    c = pd.DataFrame(cands)
    c["act"] = c.n_m + c.n_o
    audited_cells = {f.rsplit("_", 1)[0] for f in done}
    c["fresh"] = ~c.cell.isin(audited_cells)
    print(f"{len(c)} candidate windows (not yet audited, >= {a.min_duration:g} s recordings) in {c.file.nunique()} recordings, "
          f"{c.cell.nunique()} cells ({int(c.drop_duplicates('cell').fresh.sum())} cells with no audited window yet)")

    if a.list:
        pd.set_option("display.width", 200)
        print("\nwindows (1 s steps, so neighbours overlap) by internal x holding (mV):")
        print(c.groupby(["internal", "holding"]).agg(windows=("file", "size"), cells=("cell", "nunique"),
                                                      with_missed=("n_missed", lambda s: int((s >= 2).sum())),
                                                      with_lowfp=("n_lowfp", lambda s: int((s > 0).sum()))).to_string())
        print("\nwindows by day:"); print(c.groupby("day").agg(windows=("file", "size"), cells=("cell", "nunique")).to_string())
        return

    rules = [parse_rule(r) for r in a.rule]
    if a.rules_file:
        rules += [parse_rule(l) for l in a.rules_file.read_text().splitlines() if l.strip() and not l.strip().startswith("#")]
    if not rules: rules = default_rules(a.n_windows)

    taken, rows = {}, []
    free = lambda f, s: len(taken.get(f, [])) < a.max_per_recording and all(abs(s - t) >= WIN for t in taken.get(f, []))

    def subset(r):
        d = c
        if r["internal"]: d = d[d.internal.isin(r["internal"])]
        if r["holding"]: d = d[d.holding.isin(r["holding"])]
        if r["day"]: d = d[d.day.isin(r["day"])]
        if r["cell"]: d = d[d.cell.str.contains(r["cell"], regex=False)]
        if r["min_events"]: d = d[d.act >= r["min_events"]]
        if r["kind"] == "missed": d = d[d.n_missed >= 1]
        if r["kind"] == "fp": d = d[d.n_lowfp >= 1]
        if not a.no_fresh_preference and r["kind"] in ("random", "missed", "fp") and d[d.fresh].cell.nunique() >= 4: d = d[d.fresh]
        return d

    summary = []
    for r in rules:
        df = subset(r).copy()
        score = {"random": None, "missed": "n_missed", "fp": "n_lowfp", "active": "act"}[r["kind"]]
        df["_s"] = (df[score] if score else 0) + rng.uniform(0, 0.5, len(df))
        by_cell = {k: g.sort_values("_s", ascending=False) for k, g in df.groupby("cell")}
        order = list(by_cell); rng.shuffle(order); idx = {k: 0 for k in order}
        got = 0
        while got < r["n"] and order:           # round robin over cells so no cell takes everything
            for k in list(order):
                g = by_cell[k]
                while idx[k] < len(g) and not free(g.iloc[idx[k]].file, g.iloc[idx[k]].start_s): idx[k] += 1
                if idx[k] >= len(g): order.remove(k); continue
                x = g.iloc[idx[k]]; idx[k] += 1
                taken.setdefault(x.file, []).append(x.start_s)
                rows.append(dict(file=x.file, start_s=round(x.start_s, 3), why_chosen=WHY[r["kind"]],
                                 stratum=f"{r['tag']} | {x.internal} | {x.holding:g} mV | {r['kind']}"))
                got += 1
                if got >= r["n"]: break
        summary.append((r["tag"], r["n"], got, len(df)))
    out = pd.DataFrame(rows)
    if not len(out): sys.exit("no windows matched: check the rules (python scripts/pick_audit.py --list shows what exists)")
    out["_g"] = out.groupby("stratum").cumcount()                      # interleave strata: any first N windows are a mix
    out = out.sort_values(["_g", "stratum"], kind="stable").drop(columns="_g").reset_index(drop=True)
    name = a.name or time.strftime("pick_%Y%m%d_%H%M")
    path = a.out or (audit_dir / f"audit_{name}_windows.csv")
    path.parent.mkdir(exist_ok=True); out.to_csv(path, index=False)
    print("\nrule                                         wanted  got  (matching windows)")
    for tag, n, got, pool in summary: print(f"  {tag[:42]:42s} {n:6d} {got:4d}  ({pool})" + ("   <-- fewer than wanted" if got < n else ""))
    print(f"\n{len(out)} windows in {out.file.nunique()} recordings -> {path}")

    batch = a.batch or a.name or "picked"
    cmds = [[sys.executable, str(Path(__file__).resolve().parent / "build_audit.py"), "--master", str(a.master), "--audit-dir", str(audit_dir),
             "--data-root", str(a.data_root or a.master.parent), "--windows-from", str(path), "--batch", batch],
            [sys.executable, str(Path(__file__).resolve().parent / "build_audit_tool.py"), "--audit-dir", str(audit_dir),
             "--data-root", str(a.data_root or a.master.parent), "--base-html", str(audit_dir / "audit_tool.html"), "--out", str(audit_dir / "audit_tool.html")]]
    if not a.build:
        print("\nNext:  python scripts/build_audit.py --windows-from " + str(path) + f' --batch "{batch}"\n'
              "       python scripts/build_audit_tool.py --data-root . --base-html audit/audit_tool.html --out audit/audit_tool.html\n"
              "(or run this again with --build to do both)")
        return
    tool = audit_dir / "audit_tool.html"
    if a.dry_run:
        print("\nwould run:"); [print("  " + " ".join(shlex.quote(x) for x in cm)) for cm in cmds]; return
    if tool.exists():                       # keep the old tool: it holds the traces of all earlier windows
        bk = audit_dir / "_backups"; bk.mkdir(exist_ok=True)
        shutil.copy2(tool, bk / f"audit_tool_before_{name}_{time.strftime('%Y%m%d_%H%M')}.html")
    for cm in cmds:
        print("\n$ " + " ".join(shlex.quote(x) for x in cm), flush=True)
        if subprocess.call(cm) != 0: sys.exit("a build step failed - audit_labels.xlsx is unchanged (backups are in audit/_backups, the old audit_tool.html too). Fix the error above and run again.")
    print(f"\nDone. Open {tool} in the same browser as before, press J to jump to the first open window, export when finished.")


if __name__ == "__main__":
    main()
