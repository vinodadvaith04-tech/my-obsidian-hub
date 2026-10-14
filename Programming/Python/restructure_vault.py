#!/usr/bin/env python3
"""
Restructure the Obsidian vault (agreed with Claude on 2026-10-14).

Run it with Obsidian CLOSED (and Obsidian Sync paused):

    python3 restructure_vault.py --dry-run     # preview only, changes nothing
    python3 restructure_vault.py               # do it

Before changing anything it commits a snapshot of the vault in git, so the whole
restructure can be undone. It moves/renames notes, removes the empty Maths and
Papers folders, cleans up patch.md, adds a "Vault guide" note and updates the
folder icons. It does not edit the text of any other note.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

TRAILERS = (
    "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>\n"
    "Claude-Session: https://claude.ai/code/session_01SUZ5EBrMDLF8KwD85nRBEi"
)

SKIP_DIRS = {".git", ".obsidian"}

GUIDE = """# Vault guide

How this vault is organised (set up 2026-10-14). Claude reads this before touching the vault.

## Folders
- **Temp**: landing spot for new and random notes. Sort them out later.
- **Miscellaneous**: notes worth keeping that don't belong to any other folder.
- **PhD?**: everything PhD related (ideas, positions, applications).
- **Internship**: supervisor meetings (Sam Meetings, Corette Meetings), literature and method notes.
- **Lab Journal**: shown to supervisors, so keep it tidy. `Recordings/` has one note per recording day, plus `Staining/`, `Pilot studies/` and `Treatments/`.
- **BSc Psychology**, **MSc Cognitive Neuroscience**, **Biology**: course notes, kept as separate folders. MSc and BSc keep their Semester folders.
- **Programming**: language notes. Scripts Claude writes for you are saved here too (`Programming/<language>/`).
- **To do**: task lists.
- **z_Attachments**: pasted images (set in Obsidian's attachment settings, do not rename).

## Naming
- Dates go first and are written `YYYY-MM-DD`, for example `2026-10-13 Sacha recording`.
- Recordings: `YYYY-MM-DD Gerda recording` or `YYYY-MM-DD Sacha recording`, one note per day, all in `Lab Journal/Recordings`.
- Meetings: `YYYY-MM-DD Sam meeting`.
- Sentence case, no typos in names, and avoid special characters like `?` in note names (the `PhD?` folder is a deliberate exception).
"""


def say(msg=""):
    print(msg, flush=True)


def walk_files(root):
    for dp, dn, fn in os.walk(root):
        dn[:] = [d for d in dn if d not in SKIP_DIRS]
        for f in fn:
            yield Path(dp) / f


def md_files(root):
    return [p for p in walk_files(root) if p.suffix == ".md"]


def link_state(root):
    """Return (number of [[links]]/embeds, set of unresolved targets)."""
    stems, names = set(), set()
    for p in walk_files(root):
        names.add(p.name.lower())
        if p.suffix == ".md":
            stems.add(p.stem.lower())
    total, unresolved = 0, set()
    for p in md_files(root):
        text = p.read_text(encoding="utf-8", errors="replace")
        for m in re.finditer(r"!?\[\[([^\]|#]+)", text):
            total += 1
            x = m.group(1).strip().lower()
            base = Path(x).name
            if base not in names and Path(base).stem not in stems and x not in stems:
                unresolved.add(x)
    return total, unresolved


def build_moves(root):
    moves = [
        # root notes that belong to a course
        ("Integrative transcriptomics and electrophysiological profiling of hiPSC-derived neurons identifies novel druggable pathways in Koolen-de Vries Syndrome.md",
         "MSc Cognitive Neuroscience/Semester 1/Neurogenetics project proposal.md"),
        ("Project.md", "MSc Cognitive Neuroscience/Semester 1/Animal models project.md"),
        ("Note-qbn add.md", "MSc Cognitive Neuroscience/Semester 2/Note-qbn add.md"),
        ("Temp/Social neurocognition presentation.md",
         "MSc Cognitive Neuroscience/Semester 2/Social neurocognition presentation.md"),
        ("Lab rotations info.md", "MSc Cognitive Neuroscience/Lab rotations info.md"),
        ("Proposal.md", "MSc Cognitive Neuroscience/Proposal.md"),
        # random notes go to Temp
        ("NAM talk.md", "Temp/NAM talk.md"),
        ("Neuroscience mind map.md", "Temp/Neuroscience mind map.md"),
        ("Sensation and perception of touch.md", "Temp/Sensation and perception of touch.md"),
        # PhD
        ("xPhD Ideas.md", "PhD?/PhD ideas.md"),
        ("xPhD positions?.md", "PhD?/PhD positions.md"),
        # keep, but no natural home
        ("Temp/CV update 09-26.md", "Miscellaneous/CV update 09-26.md"),
        ("Temp/Email format.md", "Miscellaneous/Email format.md"),
        ("Temp/Gate Prep.md", "Miscellaneous/Gate Prep.md"),
        # typo fixes
        ("BSc Psychology/Semester 5/Counselling an Psychotherapy.md",
         "BSc Psychology/Semester 5/Counselling and Psychotherapy.md"),
        ("MSc Cognitive Neuroscience/Semester 2/Whats Next.md",
         "MSc Cognitive Neuroscience/Semester 2/What's Next.md"),
    ]
    # Lab journal: every recording day into one folder, rig in the title
    for rig in ("Gerda", "Sacha"):
        d = root / "Lab Journal" / rig / "Trial"
        if d.is_dir():
            for f in sorted(d.glob("*.md")):
                m = re.fullmatch(r"(\d{4})(\d{2})(\d{2})", f.stem)
                if m:
                    y, mo, da = m.groups()
                    moves.append((f"Lab Journal/{rig}/Trial/{f.name}",
                                  f"Lab Journal/Recordings/{y}-{mo}-{da} {rig} recording.md"))
    # Meetings with Sam: YYMMDD -> YYYY-MM-DD
    sam = root / "Internship" / "Sam Meetings"
    if sam.is_dir():
        for f in sorted(sam.glob("*.md")):
            m = re.fullmatch(r"(\d{2})(\d{2})(\d{2})", f.stem)
            if m:
                yy, mo, da = m.groups()
                moves.append((f"Internship/Sam Meetings/{f.name}",
                              f"Internship/Sam Meetings/20{yy}-{mo}-{da} Sam meeting.md"))
    return moves


def only_empty_stubs(folder):
    for p in folder.rglob("*"):
        if p.is_file() and p.name != ".DS_Store":
            if not (p.suffix == ".md" and p.stat().st_size == 0):
                return False
    return True


def git(root, *args, check=True):
    cmd = ["git", "-C", str(root)]
    have_name = subprocess.run(cmd + ["config", "user.name"], capture_output=True, text=True).stdout.strip()
    if not have_name:
        cmd += ["-c", "user.name=Vault restructure", "-c", "user.email=vault@localhost"]
    return subprocess.run(cmd + list(args), check=check, capture_output=True, text=True)


def commit_all(root, message):
    git(root, "add", "-A")
    if git(root, "diff", "--cached", "--quiet", check=False).returncode == 0:
        return None
    git(root, "commit", "-m", f"{message}\n\n{TRAILERS}")
    return git(root, "rev-parse", "HEAD").stdout.strip()


def clean_patch(root, dry):
    p = root / "Internship" / "patch.md"
    if not p.exists():
        return "patch.md not found, skipped"
    text = p.read_text(encoding="utf-8")
    if text.startswith("# Patch clamp interview prep"):
        return "patch.md already cleaned"
    lines = text.split("\n")
    for i, line in enumerate(lines[:40]):
        if re.match(r"\s*[\d.]+\s*KB\s*[•·]", line):
            body = "\n".join(lines[i + 1:]).lstrip("\n")
            if not dry:
                p.write_text("# Patch clamp interview prep\n\n" + body, encoding="utf-8")
            return f"patch.md: removed {i + 1} lines of pasted interface text at the top"
    return "patch.md: marker line not found, left untouched"


def update_icons(root, dry):
    p = root / ".obsidian" / "plugins" / "obsidian-icon-folder" / "data.json"
    if not p.exists():
        return "icon plugin config not found, skipped"
    data = json.loads(p.read_text(encoding="utf-8"))
    changed = []
    if "Maths" in data:
        del data["Maths"]
        changed.append("removed Maths icon")
    for folder, icon in (("PhD?", "LiGraduationCap"), ("Miscellaneous", "LiArchive")):
        if folder not in data:
            data[folder] = icon
            changed.append(f"added icon for {folder}")
    if changed and not dry:
        p.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    return "icons: " + (", ".join(changed) if changed else "nothing to change")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true", help="show what would happen, change nothing")
    ap.add_argument("--vault", help="vault folder (default: two levels above this script)")
    ap.add_argument("--no-git", action="store_true", help="skip the git snapshot (not recommended)")
    ap.add_argument("--skip-obsidian-check", action="store_true")
    a = ap.parse_args()

    root = Path(a.vault).resolve() if a.vault else Path(__file__).resolve().parents[2]
    dry = a.dry_run
    say(f"Vault: {root}" + ("   [DRY RUN, nothing will change]" if dry else ""))

    # ---------- pre-flight (nothing is changed here) ----------
    problems = []
    if not (root / ".obsidian").is_dir():
        problems.append("This does not look like an Obsidian vault (no .obsidian folder).")
    if not a.skip_obsidian_check and sys.platform == "darwin":
        if subprocess.run(["pgrep", "-x", "Obsidian"], capture_output=True).returncode == 0:
            problems.append("Obsidian is still running. Quit it (Cmd+Q), pause Sync, then run this again.")
    if not a.no_git:
        if not shutil.which("git") or not (root / ".git").exists():
            problems.append("No git repository found, so no safety snapshot is possible. (Use --no-git to override.)")

    moves = build_moves(root)
    to_do, done = [], []
    for src, dst in moves:
        s, d = root / src, root / dst
        if s.exists() and d.exists():
            problems.append(f"Both exist, refusing to overwrite: {src}  ->  {dst}")
        elif s.exists():
            to_do.append((src, dst))
        elif d.exists():
            done.append((src, dst))
        else:
            problems.append(f"Missing: {src}")
    for folder in ("Maths", "Papers"):
        f = root / folder
        if f.exists() and not only_empty_stubs(f):
            problems.append(f"{folder}/ contains something other than empty notes, not removing it.")

    if problems:
        say("\nStopping before any change:")
        for p in problems:
            say(f"  - {p}")
        sys.exit(1)

    notes_before = len(md_files(root))
    guide_existed = (root / "Vault guide.md").exists()
    links_before, unresolved_before = link_state(root)
    say(f"\nFound {notes_before} notes, {links_before} links/embeds "
        f"({len(unresolved_before)} unresolved target(s) already: {sorted(unresolved_before)}).")
    say(f"{len(to_do)} moves/renames to do, {len(done)} already done.")

    # ---------- snapshot ----------
    snapshot = None
    if not dry and not a.no_git:
        snapshot = commit_all(root, "Pre-restructure snapshot")
        if snapshot is None:
            snapshot = git(root, "rev-parse", "HEAD").stdout.strip()
        say(f"\nSnapshot saved in git: {snapshot[:10]}")

    # ---------- do the work ----------
    for folder in ("Temp", "Miscellaneous", "PhD?", "Lab Journal/Recordings"):
        if not (root / folder).exists():
            say(f"+ new folder  {folder}/")
            if not dry:
                (root / folder).mkdir(parents=True)
    for src, dst in to_do:
        say(f"  move  {src}\n    ->  {dst}")
        if not dry:
            (root / dst).parent.mkdir(parents=True, exist_ok=True)
            os.rename(root / src, root / dst)

    removed_empty_notes = 0
    for folder in ("Maths", "Papers"):
        f = root / folder
        if f.exists():
            removed_empty_notes += sum(1 for p in f.rglob("*.md"))
            say(f"- remove empty folder  {folder}/")
            if not dry:
                shutil.rmtree(f)
    for folder in ("Lab Journal/Gerda/Trial", "Lab Journal/Gerda",
                   "Lab Journal/Sacha/Trial", "Lab Journal/Sacha"):
        f = root / folder
        if f.exists():
            leftovers = [p for p in f.iterdir() if p.name != ".DS_Store"]
            if leftovers and not dry:
                say(f"! left {folder}/ in place, it still contains: {[p.name for p in leftovers]}")
            else:
                say(f"- remove emptied folder  {folder}/")
                if not dry:
                    ds = f / ".DS_Store"
                    if ds.exists():
                        ds.unlink()
                    f.rmdir()

    say("\n" + clean_patch(root, dry))
    guide = root / "Vault guide.md"
    if not guide.exists():
        say("+ new note  Vault guide.md")
        if not dry:
            guide.write_text(GUIDE, encoding="utf-8")
    say(update_icons(root, dry))

    if dry:
        say("\nDry run finished. Nothing was changed.")
        return

    # ---------- verify ----------
    notes_after = len(md_files(root))
    links_after, unresolved_after = link_state(root)
    loose = sorted(p.name for p in root.glob("*.md") if p.name != "Vault guide.md")
    say("\n=== Check ===")
    added_guide = 0 if guide_existed else 1
    say(f"Notes: {notes_before} before, {notes_after} after "
        f"({removed_empty_notes} empty stubs removed with Maths/Papers, {added_guide} guide note added).")
    say(f"Links/embeds: {links_before} before, {links_after} after.")
    say(f"Unresolved targets unchanged: {unresolved_after == unresolved_before}  {sorted(unresolved_after)}")
    say(f"Loose notes left in the vault root: {loose if loose else 'none'}")
    ok = (notes_after == notes_before - removed_empty_notes + added_guide
          and links_after == links_before and unresolved_after == unresolved_before and not loose)
    say("All checks passed." if ok else "!! Something differs from what was expected. Do not reopen Obsidian; tell Claude.")

    final = None
    if not a.no_git:
        final = commit_all(root, "Restructure vault: Temp, Miscellaneous and PhD? folders, "
                                 "Lab Journal/Recordings, ISO dates in names")
        say(f"\nRestructure committed in git: {final[:10] if final else 'no changes'}")
        say(f"To undo everything:  git -C \"{root}\" reset --hard {snapshot}")
    say("\nDone. You can reopen Obsidian now.")


if __name__ == "__main__":
    main()
