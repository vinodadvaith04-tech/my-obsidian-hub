import sys; sys.path.insert(0,'/home/claude/work/style')
import numpy as np, pandas as pd
from figstyle import *
setup(); K='/home/claude/work/k60/'; O='/home/claude/work/style/out/'
GE='#7a8a99'; GI='#b9987a'
# ===== Figure 1
ch=pd.read_csv(K+'mid/mid_cell_holding.csv'); ref=pd.read_csv(K+'k60_event_freq_per_cell.csv')
H=[-50,-45,-40,-35,-30]
E=[('-60',ref[ref.grp=='-60 mV'].f.values,GE)]+[(str(h),ch[ch.holding_command_mV==h].f_in.dropna().values,BL) for h in H]
I=[(str(h),ch[ch.holding_command_mV==h].f_out.dropna().values,RD) for h in H]+[('0',ref[ref.grp=='0 mV'].f.values,GI),('+5 to +15',ref[ref.grp=='+5 to +15 mV'].f.values,GI)]
fig,axs=plt.subplots(2,2,figsize=(11.6,8.6),gridspec_kw=dict(width_ratios=[6,7]))
yl='miniML event frequency per cell (Hz)'; xl='Command holding (mV)'
box_panel(axs[0,0],E,full_ylim(E),yl,'A  Inward events (EPSC), full scale',xlabel=xl)
box_panel(axs[0,1],I,full_ylim(I),yl,'B  Outward events (IPSC), full scale',xlabel=xl)
box_panel(axs[1,0],E,(0,2.0),yl,'C  Inward events, zoom to 0-2 Hz',xlabel=xl,zoom=True)
box_panel(axs[1,1],I,(0,.75),yl,'D  Outward events, zoom to 0-0.75 Hz',xlabel=xl,zoom=True)
axs[0,0].annotate('6 Nov C2\n7.5 Hz',xy=(0.05,7.5),xytext=(0.55,6.9),fontsize=8,arrowprops=dict(arrowstyle='-',color='#777',lw=.8),va='center')
legend(axs[0,1],'upper left',(0.0,1.0),extra=[Line2D([0],[0],color='none',label='↑k = k points above the zoom axis')])
foot(fig,'Greyed boxes: -60 mV (A, C) and 0, +5 to +15 mV (B, D) from the first report, for reference. Outward events at -60 mV and inward events at 0 mV were not analysed by miniML.\nLower row: same data as the upper row on a narrower scale.',.004)
save(fig,O+'Kmid_fig1.png',rect=(0,.05,1,1))
# ===== Figure 2
g=pd.read_csv(K+'mid/mid_master_table.csv'); a=pd.read_csv(K+'k60_analysis_set3.csv'); a=a[a.grp=='-60 mV']
def cell(d): return d.groupby('cellkey').agg(noise=('noise_lp3k_pA','mean'),drift=('drift_total_pA','mean'))
D={'-60':cell(a)}
for h in H: D[str(h)]=cell(g[g.holding_command_mV==h])
order=['-60']+[str(h) for h in H]
def grp(m,col): return [(k,D[k][m].values,GE if k=='-60' else col) for k in order]
fig,axs=plt.subplots(1,3,figsize=(13.2,5.0))
gn=grp('noise',BL); gd=grp('drift',RD)
box_panel(axs[0],gn,full_ylim(gn),'Noise per cell (pA)','A  Noise (MAD-SD, 3 kHz low-pass)',xlabel=xl)
box_panel(axs[1],gd,full_ylim(gd,zero=False),'Drift per cell (pA over the recording)','B  Drift over 2 min, full scale',xlabel=xl,zero=True)
box_panel(axs[2],gd,(-60,60),'Drift per cell (pA over the recording)','C  Drift, zoom to ±60 pA',xlabel=xl,zero=True,zoom=True)
legend(axs[1],'lower left',(0.0,0.0),extra=[Line2D([0],[0],color='none',label='↑k / ↓k = k points beyond the zoom axis')])
foot(fig,'Grey = -60 mV cells from the first report. Each dot is one cell (mean of its recordings at that holding). Drift is signed: negative = current falling over the 2 min.'); save(fig,O+'Kmid_fig2.png')
print('ok')
