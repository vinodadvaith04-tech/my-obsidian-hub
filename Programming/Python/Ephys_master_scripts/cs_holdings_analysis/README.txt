Cs-gluconate holdings analysis (EPSC/IPSC isolation and E/I ratio) - run order
1. cs_holdings_extract.py   (run on the Mac from the Ephys folder: seal test, noise, drift, transients -> _work/cs_holdings/)
2. analyse_cs.py            (merge features + miniML both directions + second detector -> cs_master_table.csv)
3. stats_cs.py, stats_cs2.py (cell x holding table, Tables 1-4), stats_cs3.py (isolation, E/I ratio, matching holding)
4. fig_cs.py                (Figures 1-3, box stats), build_cs_tables.py (Tables 1-7 as markdown)
analyse_cs.py reads pickles exported from Ephys_master.xlsx (files, miniML_per_file, alt_det_per_file) from /tmp; adjust paths if re-running.
