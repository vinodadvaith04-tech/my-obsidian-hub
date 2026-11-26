Shared figure style for the four ephys reports (K-Glu -60/0, K-Glu intermediate, Cs-gluconate, Lab standards).
figstyle.py   shared helpers: one box-plot look (box = IQR, thick line = median, diamond = mean, whiskers = min-max,
              one dot per cell/recording, n under each box), linear axes only, full-scale + zoom panels (up/down-k counts
              replace the old clipping triangles).
kglu60.py, kglumid.py, cs.py, cs_opts.py, lab.py  rebuild the figures of each report from the existing analysis CSVs.
Data paths in these scripts point at the cloud working folders (/home/claude/work/...); to run on the Mac, change the
folder constants at the top of each script (K, C, L, O) to the matching folders under analysis/.
Statistics and tables were not changed; only the figures.

NOISE STANDARD (decided 26 Nov 2026, use for all reports and recordings):
  Noise = peak-to-peak = 6 x robust SD (MAD-SD) of the 3 kHz low-passed trace after a 100 ms block-median detrend,
  taken over the whole 2-min recording. Column noise_lp3k_pA in the analysis CSVs is the MAD-SD; noise_p2p = 6 * that.
  Axon Guide: peak-to-peak ~ 6 x RMS. The old 4 pA MAD-SD threshold is 24 pA peak-to-peak.
kglu60.py: Figure 2 leaves out 6 Nov cell 2 (EXCLUDE_FIG2; 7.5 Hz, not physiological for an EPSC) and has a 0-1 Hz zoom.
