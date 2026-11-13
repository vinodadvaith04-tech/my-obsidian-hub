Scripts for Ephys_master.xlsx  (all read and write the ONE master workbook; close it in Excel before running)
Run from the Ephys folder:  python scripts/<name>.py

ephys_master.py          shared helpers (not run directly)
build_master.py          one-off: how Ephys_master.xlsx was first built from the older tables
add_new_recordings.py    finds new .abf files and adds rows to sheet "files" (you then fill in holding, internal, use_for_miniML = Y)
miniml_batch_master.py   runs miniML on every row marked Y that is not in the master yet; adds results to the master
run_alt_detector.py      same for the second detector (cross-check; not for sharing)
build_audit.py           adds audit windows/pages to audit/ (appends, keeps your labels)
score_audit.py           scores audit/audit_labels.xlsx
cs_report.py             rebuilds Cs_pilot_package/Cs_pilot_tables.xlsx and the figures from the master
parse_qc.py              reads manual QC checks from the lab journal into sheet manual_qc
csv_to_xlsx.py           makes readable .xlsx copies of any .csv
psc_detection.py         the second detector's code
Every write makes a backup in _backups/ first (last 5 kept).
