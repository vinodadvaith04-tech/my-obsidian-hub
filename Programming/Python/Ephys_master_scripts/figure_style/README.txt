Shared figure style for the four ephys reports (K-Glu -60/0, K-Glu intermediate, Cs-gluconate, Lab standards).
figstyle.py   shared helpers: one box-plot look (box = IQR, thick line = median, diamond = mean, whiskers = min-max,
              one dot per cell/recording, n under each box), linear axes only, full-scale + zoom panels (up/down-k counts
              replace the old clipping triangles).
kglu60.py, kglumid.py, cs.py, cs_opts.py, lab.py  rebuild the figures of each report from the existing analysis CSVs.
Data paths in these scripts point at the cloud working folders (/home/claude/work/...); to run on the Mac, change the
folder constants at the top of each script (K, C, L, O) to the matching folders under analysis/.
Statistics and tables were not changed; only the figures.
