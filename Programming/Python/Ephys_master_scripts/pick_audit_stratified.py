"""
pick_audit_stratified.py - an UNBIASED audit round: random 3 s windows, spread evenly over cells, in a fixed number per
stratum (internal x holding class). Meant for estimating what miniML finds and misses, so that frequencies can be corrected.

  python scripts/pick_audit_stratified.py --name round6 --list-plan          show the plan and what is already audited
  python scripts/pick_audit_stratified.py --name round6 --dry-run            pick + write audit/audit_round6_windows.csv only
  python scripts/pick_audit_stratified.py --name round6 --build              pick, then build_audit.py + build_audit_tool.py (backs up the tool first)

How the windows are chosen (differs from pick_audit.py's kind=random):
  * strata = internal x holding class: negative (<= -55 mV), intermediate (-50 to -30 mV), positive (>= 0 mV)
  * the number per stratum is a plan (--plan "Cs-gluconate:neg=12,int=25,pos=20;KCl:neg=6;KGlu:neg=40,int=32,pos=20"); the default plan is
    the one in DEFAULT_PLAN below (more windows where earlier rounds left few true events: K-Glu, Cs intermediate/positive)
  * inside a stratum every cell gets one window before any cell gets a second (cells with the FEWEST random windows audited so far first), max --max-per-cell (3) per cell;
    inside a cell a recording is drawn in proportion to how many free windows it has; at most --max-per-recording windows per recording
  * never reuses or overlaps an audited / already queued window (audit_labels.xlsx, sheet audit_windows)
  * every window is "random" (why_chosen=random), i.e. NO selection on what the detectors found. The csv is ordered so that the first N windows
    are already a stratified mix: you can stop early.
Needs: pandas numpy openpyxl (+ pyabf scipy matplotlib for --build)
"""
import argparse, shlex, subprocess, sys, time, shutil
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ephys_master as em

WIN, STEP = 3.0, 1.0
DEFAULT_PLAN = "Cs-gluconate:neg=12,int=25,pos=20;KCl:neg=6;KGlu:neg=40,int=32,pos=20"
HC = {"neg": "negative (<= -55 mV)", "int": "intermediate (-50 to -30 mV)", "pos": "positive (>= 0 mV)"}


def hclass(h): return "pos" if h >= 0 else ("neg" if h <= -55 else "int")


def parse_plan(s):
    plan = {}
    for part in s.split(";"):
        if not part.strip(): continue
        name, items = part.split(":")
        for it in items.split(","):
            k, v = it.split("="); plan[(name.strip(), k.strip())] = int(v)
    return plan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", type=Path, default=em.default_master())
    ap.add_argument("--audit-dir", type=Path, default=None)
    ap.add_argument("--data-root", type=Path, default=None)
    ap.add_argument("--name", required=True)
    ap.add_argument("--plan", default=DEFAULT_PLAN)
    ap.add_argument("--seed", type=int, default=6)
    ap.add_argument("--min-duration", type=float, default=30.0)
    ap.add_argument("--max-per-recording", type=int, default=2)
    ap.add_argument("--max-per-cell", type=int, default=3)
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--list-plan", action="store_true")
    ap.add_argument("--out", type=Path, default=None)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    audit_dir = a.audit_dir or (a.master.parent / "audit")
    plan = parse_plan(a.plan)

    runs = em.read_sheet(a.master, "runs")
    h = runs[runs.tool == "miniML"].sort_values("last_processed").params_hash.iloc[-1]
    pf = em.read_sheet(a.master, "miniML_per_file"); pf = pf[pf.params_hash == h].drop_duplicates("new_name")
    pf = pf[pf.new_name.isin(set(em.read_sheet(a.master, "alt_det_per_file").new_name)) & (pf.duration_s >= a.min_duration)].copy()
    pf["cell"] = pf.new_name.str.rsplit("_", n=1).str[0]
    pf["hc"] = [hclass(x) for x in pf.holding_command_mV]

    done = {}
    w = pd.read_excel(audit_dir / "audit_labels.xlsx", sheet_name="audit_windows")
    for f, s in zip(w.file, w.start_s): done.setdefault(f, []).append(float(s))
    w["cell"] = w.file.str.rsplit("_", n=1).str[0]
    w["hc"] = [hclass(x) for x in w.holding_mV]
    n_rand = w[w.why_chosen == "random"].groupby("cell").size().to_dict()          # random windows audited per cell so far

    # free candidate windows (1 s steps; not overlapping an audited window)
    cands = []
    for r in pf.itertuples():
        for s in np.arange(0, r.duration_s - WIN, STEP):
            if any(abs(s - t) < WIN for t in done.get(r.new_name, [])): continue
            cands.append((r.new_name, float(s), r.internal, r.hc, r.holding_command_mV, r.cell))
    c = pd.DataFrame(cands, columns=["file", "start_s", "internal", "hc", "holding", "cell"])

    if a.list_plan:
        t = pf.groupby(["internal", "hc"]).agg(recordings=("new_name", "size"), cells=("cell", "nunique"), minutes=("duration_s", lambda s: round(s.sum() / 60)))
        t["random_windows_audited"] = w[w.why_chosen == "random"].groupby(["internal", "hc"]).size()
        t["planned_new"] = [plan.get(k, 0) for k in t.index]
        print(t.fillna(0).astype({"random_windows_audited": int}).to_string()); print("planned total:", sum(plan.values())); return

    taken, rows, summary = {}, [], []
    free = lambda f, s: len(taken.get(f, [])) < a.max_per_recording and all(abs(s - t) >= WIN for t in taken.get(f, []))
    for (internal, hcn), n in plan.items():
        d = c[(c.internal == internal) & (c.hc == hcn)]
        cells = list(d.cell.unique()); got_cell = {k: 0 for k in cells}; got = 0
        while got < n and cells:
            # one window per cell first (cells with the fewest random windows audited so far first), then a second round, ...
            key = lambda k: (got_cell[k], n_rand.get(k, 0), rng.random())
            k = min(cells, key=key)
            dd = d[(d.cell == k)]
            dd = dd[[free(f, s) for f, s in zip(dd.file, dd.start_s)]]
            if got_cell[k] >= a.max_per_cell or not len(dd): cells.remove(k); continue
            x = dd.iloc[rng.integers(len(dd))]       # drawn over free windows => recordings in proportion to their free time
            taken.setdefault(x.file, []).append(x.start_s); got_cell[k] += 1; got += 1
            rows.append(dict(file=x.file, start_s=round(x.start_s, 3), why_chosen="random", stratum=f"Random {a.name} | {internal} | {x.holding:g} mV | random",
                             _stratum=f"{internal}|{hcn}", cell=k))
        summary.append((internal, hcn, n, got, int(d.cell.nunique()), len(d)))
    out = pd.DataFrame(rows)
    out["_g"] = out.groupby("_stratum").cumcount()
    out = out.sort_values(["_g", "_stratum"], kind="stable").reset_index(drop=True)      # interleave strata: any first N windows are a mix
    path = a.out or (audit_dir / f"audit_{a.name}_windows.csv")
    out[["file", "start_s", "why_chosen", "stratum"]].to_csv(path, index=False)
    print("stratum                         planned  got  cells with free windows  free windows")
    for i, hcn, n, g, nc, nw in summary: print(f"  {i:12s} {hcn:4s}  {n:10d} {g:4d} {nc:12d} {nw:20d}" + ("   <-- fewer than planned" if g < n else ""))
    print(f"\n{len(out)} windows in {out.file.nunique()} recordings, {out.cell.nunique()} cells -> {path}")
    out.assign(hold_class=out._stratum.str.split("|").str[1]).drop(columns=["_g"]).to_csv(Path(str(path).replace("_windows.csv", "_plan_check.csv")), index=False)

    if not a.build:
        print("\nNext: run again with --build (adds the windows to audit_labels.xlsx and rebuilds audit_tool.html)"); return
    tool = audit_dir / "audit_tool.html"
    cmds = [[sys.executable, str(Path(__file__).resolve().parent / "build_audit.py"), "--master", str(a.master), "--audit-dir", str(audit_dir),
             "--data-root", str(a.data_root or a.master.parent), "--windows-from", str(path), "--batch", a.name],
            [sys.executable, str(Path(__file__).resolve().parent / "build_audit_tool.py"), "--audit-dir", str(audit_dir),
             "--data-root", str(a.data_root or a.master.parent), "--base-html", str(tool), "--out", str(tool)]]
    if a.dry_run:
        print("\nwould run:"); [print("  " + " ".join(shlex.quote(x) for x in cm)) for cm in cmds]; return
    bk = audit_dir / "_backups"; bk.mkdir(exist_ok=True)
    for f in (tool, audit_dir / "audit_labels.xlsx"):
        if f.exists(): shutil.copy2(f, bk / f"{f.stem}_before_{a.name}_{time.strftime('%Y%m%d_%H%M')}{f.suffix}")
    for cm in cmds:
        print("\n$ " + " ".join(shlex.quote(x) for x in cm), flush=True)
        if subprocess.call(cm) != 0: sys.exit("a build step failed - backups are in audit/_backups")
    print(f"\nDone. Open {tool} in the same browser as before, press J to jump to the first open window, export when finished.")


if __name__ == "__main__":
    main()
