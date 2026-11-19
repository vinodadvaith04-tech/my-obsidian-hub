#!/usr/bin/env python3
"""
02_run_detection.py  -  STEP 2 of the pipeline: run miniML and the second detector on the recordings that are new, and
put the results back into the master workbook.

What it does, in this order
  1. checks everything is in place (master closed in Excel, python packages, TensorFlow, miniML folder + model) and tells
     you ALL problems at once instead of crashing halfway
  2. reads the master (sheet files) and works out what is still to do:
         rows with use_for_miniML = Y, found on disk, with a holding AND an internal solution
         minus what is already in miniML_per_file (same model/settings) / alt_det_per_file
     rows that are Y but have no holding or no internal solution WAIT (listed in logs/step2_waiting_*.csv) - fill them in
     Excel, run step 2 again, and they are picked up
  3. copies the master to _backups/pipeline/ (before the detectors touch it)
  4. runs miniML  (miniml_batch_master.py - kept as its own program) and then the second detector (run_alt_detector.py),
     each as a separate process: if one fails, the other still runs and you get a clear message for each
  5. checks that every requested recording now has a result in the master and prints what is missing and why

miniML itself is NOT part of this file: step 2 only calls  miniml_batch_master.py  with the model from pipeline_config.json
(miniml_model / miniml_core). A transfer-learned model later = change that path (or the miniML script); nothing else changes.
Detection can be stopped (create a file named STOP in _work/miniML_per_file/ or _work/alt_det_per_file/) or limited with
--max-hours; run the same command again to continue, finished recordings are cached.

Usage:
    python scripts/02_run_detection.py --dry-run          # show what would run and what is waiting
    python scripts/02_run_detection.py                    # run everything that is new
  Options: --only TEXT (name contains)   --limit N (test with N recordings)   --max-hours H   --skip-miniml   --skip-alt
           --redo (run again even if the master already has a result)   --master <file>   --data-root <folder>
Needs: pandas openpyxl pyabf  + for miniML: tensorflow (and miniML-main/core, models/GC_lstm_model.h5)
"""
from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ephys_master as em
import pipeline_common as pc

MINIML_S_PER_RUN = 40      # rough, 2 min recording, 1 direction (cloud 2 cores); a modern Mac is usually faster
ALT_S_PER_FILE = 20


def find_miniml(cfg, root):
    core = Path(cfg["miniml_core"]).expanduser() if cfg.get("miniml_core") else None
    if not core or not (core / "miniML.py").exists():
        cands = []
        for pat in ("~/Downloads/miniML-main*/core", "~/Downloads/*/miniML-main*/core", "~/miniML-main*/core", f"{root.parent}/miniML-main*/core", f"{root}/miniML-main*/core"):
            cands += glob.glob(os.path.expanduser(pat))
        core = next((Path(c) for c in sorted(cands) if (Path(c) / "miniML.py").exists()), None)
    model = Path(cfg["miniml_model"]).expanduser() if cfg.get("miniml_model") else None
    if (not model or not model.exists()) and core:
        m = core.parent / "models" / "GC_lstm_model.h5"
        model = m if m.exists() else model
    return core, model


def run_checked(cmd, env, log, tag):
    """Run a command, stream its output to the screen and the log, return (exit code, lines)."""
    log(f"\n  >>> {tag}: {' '.join(str(c) for c in cmd[:3])} ...")
    lines = []
    p = subprocess.Popen([str(c) for c in cmd], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env, bufsize=1)
    for line in p.stdout:
        line = line.rstrip("\n")
        lines.append(line)
        if line.startswith(("[", "Master", "Detection", "Settings", "Done", "recordings", "master updated", "STOP", "--max", "WARNING", "ERROR", "Traceback")) or "Error" in line:
            log("     " + line)
        else:
            log.file_only("     " + line)
    p.wait()
    return p.returncode, lines


def run(a, log):
    cfg = pc.load_config(a.config)
    root, mp = pc.resolve_paths(cfg, a.master, a.data_root)
    py = cfg.get("python") or sys.executable
    stamp = time.strftime("%Y%m%d_%H%M%S")
    log(f"STEP 2 - run detection   {time.strftime('%Y-%m-%d %H:%M')}")
    log(f"  Ephys folder : {root}\n  master       : {mp}")

    # ---------------- 1. preflight ----------------
    problems = []
    try:
        pc.check_master_free(mp)
    except pc.UserProblem as e:
        problems.append(f"{e}  -> {e.hint}")
    for mod in ("numpy", "pandas", "openpyxl", "pyabf"):
        try:
            __import__(mod)
        except ImportError:
            problems.append(f"python package '{mod}' missing  -> pip install {mod}")
    core = model = None
    env = dict(os.environ)
    if not a.skip_miniml:
        core, model = find_miniml(cfg, root)
        if core is None:
            problems.append("miniML folder not found  -> set \"miniml_core\" in scripts/pipeline_config.json to the miniML-main/core folder")
        if model is None or not Path(model).exists():
            problems.append("miniML model file not found  -> set \"miniml_model\" in scripts/pipeline_config.json (GC_lstm_model.h5)")
        t = subprocess.run([py, "-c", "import tensorflow" ], capture_output=True, text=True, env=env)
        if t.returncode != 0:
            problems.append("TensorFlow cannot be imported by " + py + "  -> install it (pip install tensorflow), or use --skip-miniml to run only the second detector. "
                            + (t.stderr.strip().splitlines()[-1] if t.stderr.strip() else ""))
        else:
            tk = subprocess.run([py, "-c", "import tf_keras"], capture_output=True, text=True, env=env)
            if tk.returncode == 0 and "TF_USE_LEGACY_KERAS" not in env:
                env["TF_USE_LEGACY_KERAS"] = "1"        # miniML's model files need the legacy Keras
        env.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
    if not a.skip_alt and not (HERE / "psc_detection.py").exists():
        problems.append("psc_detection.py not found next to this script")
    if not (HERE / "miniml_batch_master.py").exists() or not (HERE / "run_alt_detector.py").exists():
        problems.append("miniml_batch_master.py or run_alt_detector.py missing from the scripts folder")
    if problems:
        log("\n  CANNOT START - please fix:")
        for p in problems:
            log("   - " + p)
        return 2
    if core:
        log(f"  miniML       : {core}\n  model        : {model}")

    # ---------------- 2. what is to do ----------------
    files = em.read_sheet(mp, "files")
    y = files[files["use_for_miniML"].astype(str).str.upper() == "Y"].copy()
    idx, _ = em.index_abf(root)
    y["_path"] = [em.locate(r, root, idx) for _, r in y.iterrows()]
    y["_why"] = ""
    y.loc[y["_path"].isna(), "_why"] += "file not found on disk; "
    y.loc[y["holding_command_mV"].isna(), "_why"] += "no holding_command_mV; "
    y.loc[y["internal"].isna() | (y["internal"].astype(str).str.strip() == ""), "_why"] += "no internal solution; "
    seen = set()
    for sh in ("miniML_per_file", "alt_det_per_file"):
        t = em.read_sheet(mp, sh)
        if len(t) and "new_name" in t.columns:
            seen |= set(t["new_name"])
    waiting = y[(y["_why"] != "") & ~y["new_name"].isin(seen)]     # rows that already have results are not "waiting"
    ready = y[y["_why"] == ""].copy()
    if a.only:
        ready = ready[ready["new_name"].str.contains(a.only, regex=False)]

    sys.path.insert(0, str(HERE))
    import miniml_batch_master as mm
    ready["_dirs"] = ready.apply(lambda r: (str(r["miniML_directions"]).strip() if str(r.get("miniML_directions")).strip() not in ("", "nan", "None") else mm.default_directions(r)), axis=1)
    h = mm.params_hash(mm.PARAMS, Path(model)) if (model and Path(model).exists()) else None
    done_ml, done_alt = set(), set()
    if h and not a.redo:
        pf = em.read_sheet(mp, "miniML_per_file")
        if len(pf):
            sel = pf[pf.params_hash == h]
            done_ml = set(zip(sel.new_name, sel.direction))
    if not a.redo:
        af = em.read_sheet(mp, "alt_det_per_file")
        done_alt = set(af["new_name"]) if len(af) and "new_name" in af.columns else set()
    todo_ml = [(r["new_name"], d.strip()) for _, r in ready.iterrows() for d in r["_dirs"].split(",") if (r["new_name"], d.strip()) not in done_ml]
    todo_alt = [r["new_name"] for _, r in ready.iterrows() if r["new_name"] not in done_alt]
    if a.limit:
        keep = list(dict.fromkeys([n for n, _ in todo_ml] + todo_alt))[: a.limit]
        todo_ml = [t for t in todo_ml if t[0] in keep]; todo_alt = [n for n in todo_alt if n in keep]
    names_ml = list(dict.fromkeys(n for n, _ in todo_ml))
    log(f"\n  recordings marked Y in the master: {len(y)} | ready (holding + internal + file found): {len(ready)} | WAITING: {len(waiting)}")
    log(f"  miniML runs to do: {len(todo_ml)} ({len(names_ml)} recordings; settings hash {h}) | second detector: {len(todo_alt)} recordings")
    est = (0 if a.skip_miniml else len(todo_ml) * MINIML_S_PER_RUN) + (0 if a.skip_alt else len(todo_alt) * ALT_S_PER_FILE)
    log(f"  rough time: {est / 60:.0f} min")
    if len(waiting):
        wf = root / "logs" / f"step2_waiting_{stamp}.csv"
        wf.parent.mkdir(exist_ok=True)
        waiting[["new_name", "session", "slice", "cell", "_why"]].rename(columns={"_why": "waiting_for"}).to_csv(wf, index=False)
        summ = waiting.groupby("_why").size()
        for k, v in summ.items():
            log(f"    waiting ({v}): {k.strip('; ')}")
        log(f"    full list: {wf}  - fill these columns in the files sheet, then run step 2 again")
    if a.dry_run:
        log("\n  DRY RUN finished (nothing was run).")
        return 0
    if not todo_ml and not todo_alt:
        log("\n  Nothing new to analyse.")
        return 0

    # ---------------- 3. backup ----------------
    bk = pc.backup_master(mp, cfg, "step2", root)
    log(f"  backup       : {bk}")
    t0 = time.time()
    status = {}
    maxh = ["--max-hours", str(a.max_hours)] if a.max_hours else []
    redo = ["--redo"] if a.redo else []

    # ---------------- 4a. miniML ----------------
    if not a.skip_miniml and todo_ml:
        cmd = [py, HERE / "miniml_batch_master.py", "--master", mp, "--data-root", root, "--miniml-core", core, "--model", model,
               "--names", ",".join(names_ml), "--notes", f"step 2 run {stamp}"] + cfg.get("miniml_extra_args", []) + maxh + redo
        code, lines = run_checked(cmd, env, log, "miniML")
        errs = [l for l in lines if "-> ERROR" in l]
        status["miniML"] = (code, errs, lines)
    # ---------------- 4b. second detector ----------------
    if not a.skip_alt and todo_alt:
        cmd = [py, HERE / "run_alt_detector.py", "--master", mp, "--data-root", root, "--names", ",".join(todo_alt)] + maxh + redo
        code, lines = run_checked(cmd, env, log, "second detector")
        errs = [l for l in lines if "-> ERROR" in l or l.startswith("MISSING")]
        status["alt"] = (code, errs, lines)

    # ---------------- 5. check ----------------
    log("\n  ---------------- CHECK ----------------")
    rc = 0
    if "miniML" in status:
        pf = em.read_sheet(mp, "miniML_per_file")
        got = set(zip(pf[pf.params_hash == h].new_name, pf[pf.params_hash == h].direction)) if len(pf) else set()
        miss = [t for t in todo_ml if t not in got]
        code, errs, _l = status["miniML"]
        log(f"  miniML: {len(todo_ml) - len(miss)}/{len(todo_ml)} runs are in the master" + (f" | process exit code {code}" if code else ""))
        for e in errs[:15]:
            log("    " + e)
        if miss:
            rc = 1
            log(f"    NOT in the master yet: {len(miss)} (e.g. {miss[0][0]} [{miss[0][1]}]).")
            log("    Reasons: stopped by --max-hours / STOP file (just run step 2 again), a per-file error above, or the process crashed"
                + (f" (exit code {code}; see the log file {log.path})" if code else "") + ".")
    if "alt" in status:
        af = em.read_sheet(mp, "alt_det_per_file")
        got = set(af.new_name) if len(af) else set()
        miss = [n for n in todo_alt if n not in got]
        code, errs, _l = status["alt"]
        log(f"  second detector: {len(todo_alt) - len(miss)}/{len(todo_alt)} recordings are in the master" + (f" | process exit code {code}" if code else ""))
        for e in errs[:15]:
            log("    " + e)
        if miss:
            rc = 1
            log(f"    NOT in the master yet: {len(miss)} (e.g. {miss[0]}). Run step 2 again to continue; per-file errors are listed above.")
    for k, (code, _e, ls) in status.items():
        if code not in (0, None):
            tail = [l for l in ls if l.strip()][-8:]
            log(f"  {k} stopped with exit code {code}. Last lines of its output:")
            for l in tail:
                log("      " + l[:200])
            h2 = pc.hint_for("\n".join(tail))
            log("  What to try: " + (h2 or "see the log file; typical causes are the TensorFlow/Keras version, out of memory, or the master open in Excel."))
    log(f"\n  FINISHED in {(time.time() - t0) / 60:.1f} min." + ("" if rc == 0 else "  Some recordings are not done (see above)."))
    return rc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--master", type=Path)
    ap.add_argument("--data-root", type=Path)
    ap.add_argument("--config", type=Path)
    ap.add_argument("--only"); ap.add_argument("--limit", type=int); ap.add_argument("--max-hours", type=float)
    ap.add_argument("--skip-miniml", action="store_true"); ap.add_argument("--skip-alt", action="store_true")
    ap.add_argument("--redo", action="store_true"); ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    cfg = pc.load_config(a.config)
    root, _ = pc.resolve_paths(cfg, a.master, a.data_root)
    log = pc.new_log(root, "step2")
    try:
        sys.exit(run(a, log))
    except SystemExit:
        raise
    except BaseException as e:
        sys.exit(pc.report_error(log, "step 2 (detection)", e))


if __name__ == "__main__":
    main()
