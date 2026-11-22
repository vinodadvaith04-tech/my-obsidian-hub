Lab standard values: EPSC/IPSC frequency and E/I ratio (22 Nov 2026)
Past lab (Peerboom et al. 2023 controls) vs current lab (Helena, Sam, Advaith in Sam's experiment) vs Advaith Cs-gluconate.
Run order: paperdots.py -> paperdots2.py (digitise paper Fig 3C,G / 6E,I from the PDF) -> build_lab_data.py (per-cell table from Ephys_master.xlsx + Cs analysis)
 -> stats_lab.py -> build_tables.py -> fig_lab.py -> build_workbook.py
Inputs: Ephys_master.xlsx (miniML_per_cell_holding, sam_EPSC_control_cells, alt_det_per_file), Cs analysis csv (cs_cell_holding.csv, cs_ei_per_cell_ref-80_0.csv) from the Cs holdings analysis, the paper PDF.
Paths in the scripts are the cloud-session paths (/mnt/user-data/uploads, /home/claude/work/cs); edit them to run on the Mac.
Outputs: analysis/Lab_standards/ (workbook + 4 figures).
