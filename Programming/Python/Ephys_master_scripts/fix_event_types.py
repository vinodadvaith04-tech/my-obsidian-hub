"""One-off (13 Nov 2026): event_type_assumed follows the rule 'KCl without a blocker = event (unresolved)'.
Recomputes event_type_assumed in miniML_per_file, miniML_events, miniML_per_cell_holding from the files sheet; adds the drug column
to miniML_events / miniML_per_cell_holding. Backup in _backups/ first. Safe to run again."""
import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ephys_master as em
mp = em.default_master()
files = em.read_sheet(mp, "files")
drug = files.dropna(subset=["new_name"]).drop_duplicates("new_name").set_index("new_name")["drug"]
out = {}
for sh in ("miniML_per_file", "miniML_events"):
    d = em.read_sheet(mp, sh)
    d["drug"] = d["new_name"].map(drug) if "drug" not in d.columns or d["drug"].isna().all() else d["drug"]
    if sh == "miniML_events":   # keep the layout: meta columns at the end
        pass
    d["event_type_assumed"] = [em.assumed_event_type(a, b, c) for a, b, c in zip(d["direction"], d["internal"], d["drug"])]
    out[sh] = d
pf, ev = out["miniML_per_file"], out["miniML_events"]
ch = em.read_sheet(mp, "miniML_per_cell_holding")
dm = pf.drop_duplicates("recording_name").set_index("recording_name")["drug"]
ch["drug"] = ch["recording_name"].map(dm)
ch["event_type_assumed"] = [em.assumed_event_type(a, b, c) for a, b, c in zip(ch["direction"], ch["internal"], ch["drug"])]
out["miniML_per_cell_holding"] = ch
for k, v in out.items():
    print(k, v.event_type_assumed.value_counts().to_dict())
em.write_sheets(mp, out)
print("master updated")
