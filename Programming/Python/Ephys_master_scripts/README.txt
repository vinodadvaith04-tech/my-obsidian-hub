Scripts for Ephys_master.xlsx  (all read and write the ONE master workbook; close it in Excel before running)
Run from the Ephys folder (the one that contains Ephys_master.xlsx):  python scripts/<name>.py

=== THE TWO PROGRAMS YOU NORMALLY RUN, IN THIS ORDER =====================================================================
1) 01_ingest_new_recordings.py      new .abf files -> master
     python scripts/01_ingest_new_recordings.py --dry-run     (look first: nothing changes)
     python scripts/01_ingest_new_recordings.py [--internal KGlu] [--drug gabazine] [--folder "<subfolder>"]
   - safety copy of the master in _backups/pipeline/ before anything is changed
   - finds .abf files that are not in the master yet (a known file you moved only gets its new location written)
   - name already follows  <date>_Adv<S.No>_<T.Code>_<sex>_slice<N>_C<N>_<counter>.abf  -> kept (compared with the journal, differences flagged)
     otherwise -> renamed on disk, in place (S.No/T.Code/sex from the day's lab journal, else Adv00/NT/unk; slice/cell from the folder names)
   - writes files (old name, new name, old/new folder, path, protocol, duration, holding ...), cells and rename_log
   - holding per gap-free recording from the lab journal ("cap check / -80 / -80 ..." matched in order, only when the number of
     listed holdings = number of recordings of that cell); otherwise BLANK + holding_flag says why -> you type it in Excel
   - internal / drug from "Internal: ..." / "Drug: ..." lines in the journal (day header or slice header) or --internal / --drug
   - report in logs/step1_*.log, list of renames in logs/renamed_*.csv   (undo: --undo-renames logs/renamed_XXXX.csv)
2) 02_run_detection.py              miniML + second detector on what is new -> master
     python scripts/02_run_detection.py --dry-run             (shows what would run and what is WAITING for a holding/internal)
     python scripts/02_run_detection.py [--limit 3] [--max-hours 2] [--skip-miniml] [--skip-alt]
   - preflight check lists all problems at once (Excel open, packages, TensorFlow, miniML folder, model)
   - runs miniml_batch_master.py and run_alt_detector.py as separate processes; one failing does not stop the other
   - recordings without holding or internal wait; fill them in, run again
   - ends with a check that every requested recording is in the master; log in logs/step2_*.log
   Settings: scripts/pipeline_config.json (journal folder, miniML folder/model, minimum gap-free duration, defaults)

=== KEEP SEPARATE (so changing miniML later does not touch the rest) ===================================================
miniml_batch_master.py   miniML only: runs the model on rows marked Y, writes the miniML_* sheets + runs. Detection settings and the
                         direction-by-holding rule live here; the model path comes from pipeline_config.json (miniml_model).
run_alt_detector.py      second detector only (cross-check; not for sharing); code in psc_detection.py
ephys_master.py          shared helpers for reading/writing the master (not run directly)
pipeline_common.py, journal_parser.py   helpers of the two programs (config, backup, error messages / reading the lab journal)

=== OLDER / OTHER SCRIPTS ======================================================================================================
add_new_recordings.py    older: adds NEW rows to files without renaming or reading the journal (superseded by 01_)
ingest_recordings.py     older: one folder of someone else's recordings -> renamed + master rows (superseded by 01_)
build_master.py          one-off: how Ephys_master.xlsx was first built from the older tables
build_audit.py           adds audit windows/pages to audit/ (appends, keeps your labels)
build_audit_tool.py      builds audit/audit_tool.html from audit_labels.xlsx + the .abf files (needs pyabf). Re-run after build_audit.py adds windows
audit_tool_template.html the page's source (build_audit_tool.py fills in the traces)
import_audit_clicks.py   copies the JSON exported from audit_tool.html into audit_labels.xlsx (+ sheet audit_truth_events)
pick_audit_windows.py    picks audit windows from recordings already run through both detectors; feed the csv to build_audit.py --windows-from
score_audit.py           scores audit/audit_labels.xlsx
cs_report.py             rebuilds Cs_pilot_package/Cs_pilot_tables.xlsx and the figures from the master
fix_event_types.py       re-labels the assumed EPSC/IPSC event type after drug/internal were changed
parse_qc.py              reads manual QC checks from the lab journal into sheet manual_qc
csv_to_xlsx.py           makes readable .xlsx copies of any .csv
one_off/                 scripts used once for a specific batch (e.g. fill_holdings_20261118.py)

Backups: the two programs copy the master to _backups/pipeline/ (last 15 kept) before they change anything; older scripts use _backups/ (last 5).
