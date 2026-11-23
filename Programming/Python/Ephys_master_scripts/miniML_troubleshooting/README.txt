miniML troubleshooting (22 Nov 2026)
Report (Doc): https://claude.ai/artifact/8uTguftpH95JNZTuGMxJdW
Inputs: Ephys pilot & practice/audit/audit_labels.xlsx (audit_windows, audit_truth_events, audit_events_key, audit_labels).
Run order: load.py (module) -> stats.py (recall/precision per type) -> analyse2.py (overlap, holding, classifier-output bands)
           -> figs.py (Figures 1-5; paths inside point at the cloud session, edit for your machine).
param_sweep_exp.py : miniML re-run per setting on 54 audited windows (window, threshold, peak_w, rel_prom, low-pass, inversion); results in param_sweep_results.txt.
resample_test.py   : decimate x2/x4 test. classifier_output_at_misses.py : network output at audited events. event_kinetics.py : rise/decay/SNR per audited event.
The sweep and probability scripts need TensorFlow, the miniML repo and the raw .abf files; the others need only the audit workbook.
