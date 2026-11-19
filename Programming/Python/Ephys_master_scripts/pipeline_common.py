"""
pipeline_common.py - small shared helpers for the two pipeline programs
    01_ingest_new_recordings.py   (new .abf -> names, files / cells / rename_log, holdings from the lab journal)
    02_run_detection.py           (miniML + second detector on what is new, results back into the master)

Contents: config loading, a log file that also prints to the screen, the master backup, the "is Excel still open?"
check, and one place that turns an exception into a message you can act on.
Needs: pandas openpyxl (pyabf for step 1)
"""
from __future__ import annotations

import json
import os
import shutil
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_CONFIG = {
    # Folder that holds Ephys_master.xlsx and the recordings (default: the folder above scripts/)
    "data_root": None,
    "master_name": "Ephys_master.xlsx",
    # Obsidian lab journal (searched recursively for files named like 2026-11-18.md or '2026-10-07 Sacha recording.md')
    "journal_dir": "/Users/advaithvinod/Downloads/Untitled/Lab Journal",
    # Where the safety copy of the master goes before anything is changed. None = <data_root>/_backups/pipeline
    "backup_dir": None,
    "keep_backups": 15,
    # A gap-free recording shorter than this (s) is treated as an aborted attempt: not analysed, not counted against the journal
    "min_gapfree_s": 30,
    # Defaults for the file name when neither the name, the journal nor the command line says otherwise
    "default_S.No": "Adv00",
    "default_T.Code": "NT",
    "default_sex": "unk",
    # session label by the journal sub-folder a day's journal lives in
    "session_by_journal_folder": {
        "_Main Experiment": "Main experiment",
        "3. Cs internal pilot": "Cs internal pilot",
        "2. TTX pilot": "TTX pilot",
        "IPSC recordings": "Sam IPSC",
        "Longer recordings": "Long holding",
    },
    # miniML directions written to files.miniML_directions by internal (empty = the holding rule in miniml_batch_master.py decides)
    "directions_by_internal": {"Cs-gluconate": "negative,positive", "KCl": "negative"},
    # Step 2: where miniML lives. Empty = look in the usual places (Downloads/miniML-main*, next to the Ephys folder)
    "miniml_core": "",
    "miniml_model": "",
    "miniml_extra_args": [],
    "python": "",             # python used for the detectors; empty = the one that runs step 2
    "skip_dirs": ["_backups", "_work", "_archive", "_archive_before_master", "scripts", "audit", "logs", "Cs_pilot_package"],
}


def load_config(path: Path | None = None) -> dict:
    cfg = json.loads(json.dumps(DEFAULT_CONFIG))
    p = Path(path) if path else HERE / "pipeline_config.json"
    if p.exists():
        try:
            user = json.loads(p.read_text())
        except json.JSONDecodeError as e:
            raise ConfigError(f"{p} is not valid JSON ({e}). Fix the file or delete it to use the defaults.")
        for k, v in user.items():
            if k.startswith("_"):
                continue
            if isinstance(v, dict) and isinstance(cfg.get(k), dict):
                cfg[k].update(v)
            else:
                cfg[k] = v
    cfg["_config_path"] = str(p)
    return cfg


def resolve_paths(cfg: dict, master=None, data_root=None) -> tuple[Path, Path]:
    root = Path(data_root or cfg.get("data_root") or HERE.parent).expanduser().resolve()
    mp = Path(master).expanduser() if master else root / cfg["master_name"]
    return root, mp


class ConfigError(Exception):
    pass


class UserProblem(Exception):
    """A problem the user can fix (closed Excel, missing folder ...). `hint` says how."""
    def __init__(self, msg, hint=""):
        super().__init__(msg)
        self.hint = hint


# ----------------------------------------------------------------------------- log
class Tee:
    """print() to the screen and to a log file at the same time."""
    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path = path
        self.f = open(path, "a", encoding="utf-8")

    def __call__(self, *a, **k):
        s = " ".join(str(x) for x in a)
        print(s, **k, flush=True)
        self.f.write(s + ("\n" if not k.get("end") else k["end"]))
        self.f.flush()

    def file_only(self, s):
        self.f.write(s + "\n")
        self.f.flush()

    def close(self):
        self.f.close()


def new_log(root: Path, step: str) -> Tee:
    return Tee(root / "logs" / f"{step}_{time.strftime('%Y%m%d_%H%M%S')}.log")


# ----------------------------------------------------------------------------- master safety
def excel_lock_file(mp: Path) -> Path:
    return mp.with_name("~$" + mp.name)


def check_master_free(mp: Path):
    if not mp.exists():
        raise UserProblem(f"Master workbook not found: {mp}",
                          "Check --master / data_root in pipeline_config.json, or the folder you run the script from.")
    if excel_lock_file(mp).exists():
        raise UserProblem(f"{mp.name} looks open in Excel (found {excel_lock_file(mp).name}).",
                          "Close the workbook in Excel (and save nothing), then run the script again.")
    if not os.access(mp, os.W_OK):
        raise UserProblem(f"{mp} is not writable.", "Check the file permissions or whether it is on a read-only disk.")


def backup_master(mp: Path, cfg: dict, tag: str, root: Path) -> Path:
    bdir = Path(cfg["backup_dir"]).expanduser() if cfg.get("backup_dir") else root / "_backups" / "pipeline"
    bdir.mkdir(parents=True, exist_ok=True)
    dst = bdir / f"{mp.stem}_{time.strftime('%Y%m%d_%H%M%S')}_before_{tag}.xlsx"
    shutil.copy2(mp, dst)
    if dst.stat().st_size != mp.stat().st_size:
        raise UserProblem(f"The backup {dst.name} has a different size than the master; stopping before anything is changed.",
                          "Check free disk space.")
    old = sorted(bdir.glob(f"{mp.stem}_*_before_*.xlsx"))
    for p in old[: -int(cfg.get("keep_backups", 30))]:
        p.unlink()
    return dst


# ----------------------------------------------------------------------------- errors
HINTS = [
    ("ModuleNotFoundError: No module named 'pyabf'", "Install it:  pip install pyabf"),
    ("No module named 'openpyxl'", "Install it:  pip install openpyxl pandas"),
    ("No module named 'tensorflow'", "TensorFlow is not installed in this python. Install it (pip install tensorflow) or run step 2 where miniML works."),
    ("No module named 'miniML'", "miniml_core in pipeline_config.json does not point to miniML-main/core."),
    ("BadZipFile", "The master file looks damaged. Restore the newest copy from _backups/pipeline/ and tell Claude."),
    ("PermissionError", "A file is locked or read-only. Close Excel / Clampfit and check the permissions."),
    ("Keras 3", "Set the environment variable TF_USE_LEGACY_KERAS=1 (step 2 does this automatically when tf_keras is installed)."),
    ("legacy", "Set the environment variable TF_USE_LEGACY_KERAS=1 (pip install tf_keras)."),
    ("OOM", "Out of memory: close other programs, or run fewer recordings at once with --limit 10."),
]


def hint_for(text: str) -> str:
    for key, h in HINTS:
        if key.lower() in text.lower():
            return h
    return ""


def report_error(log: Tee | None, stage: str, exc: BaseException, extra: str = "", exit_code: int = 2):
    tb = traceback.format_exc()
    hint = getattr(exc, "hint", "") or hint_for(tb + str(exc))
    lines = ["", "=" * 78, f"STOPPED in stage: {stage}", f"  {type(exc).__name__}: {exc}"]
    if hint:
        lines.append(f"  What to try: {hint}")
    if extra:
        lines.append(f"  {extra}")
    if log is not None:
        lines.append(f"  Full details (copy this file to Claude if you are stuck): {log.path}")
    lines.append("  The master is only replaced at the very end of a successful run; if this run got as far as the backup step, a safety copy is in _backups/pipeline/.")
    lines.append("=" * 78)
    msg = "\n".join(lines)
    if log is not None:
        log.file_only(tb)
        log(msg)
    else:
        print(msg + "\n" + tb)
    return exit_code
