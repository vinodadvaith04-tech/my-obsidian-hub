import sys; sys.path.insert(0,'/home/claude/work/style')
import numpy as np, pandas as pd
from figstyle import *
from matplotlib.gridspec import GridSpec
setup(); K='/home/claude/work/k60/'; O='/home/claude/work/style/out/'
ONLY2='--only2' in sys.argv      # rebuild just Figure 2
# cells left out of Figure 2 (per user, 26 Nov 2026): 6 Nov slice 3 cell 2 = 7.5 Hz, 879 + 924 events in its two -60 mV recordings, not physiological for an EPSC rate
EXCLUDE_FIG2=['20261106_s3_c2']
a=pd.read_csv(K+'k60_analysis_set3.csv'); z=np.load(K+'k60_example_traces.npz')
# NOISE STANDARD (decided 26 Nov 2026): peak-to-peak = 6 x robust SD (MAD-SD of the 3 kHz low-passed trace after a 100 ms block-median detrend, whole recording)
a['noise_p2p']=6*a.noise_lp3k_pA
GRP=['-60 mV','0 mV','+5 to +15 mV']; COL={'-60 mV':BL,'0 mV':RD,'+5 to +15 mV':OC}; TL=['-60','0','+5 to +15']
if not ONLY2:
    # ================= Figure 1
    fig=plt.figure(figsize=(9.2,13.8)); gs=GridSpec(6,6,figure=fig,hspace=.95,wspace=1.6,height_ratios=[1,1,1,1,1.35,1.35],left=.09,right=.985,top=.962,bottom=.085)
    ex=[('A  -60 mV, typical','20260929_Adv00_NT_unk_slice1_C5_039.abf',BL),('B  -60 mV, unstable cell','20261007_Adv00_NT_unk_slice1_C5_048.abf',BL),
        ('C  0 mV, drifting and noisy','20261002_Adv00_NT_unk_slice1_C1_004.abf',RD),('D  0 mV, spike-like transients','20260925_Adv00_NT_unk_slice1_C5_067.abf',RD)]
    for i,(lab,nm,col) in enumerate(ex):
        ax=fig.add_subplot(gs[i,:4]); y=z[nm+'_dec']; t=np.arange(len(y))/2000; r=a[a.new_name==nm].iloc[0]
        ax.plot(t,y,lw=.4,color=col); ax.set_title(lab+'  ('+nm[:8]+', slice '+nm.split('_')[4]+' '+nm.split('_')[5]+')'+f'\nholding current {r.Ihold_median_pA:.0f} pA, drift {r.drift_total_pA:+.0f} pA over 2 min, noise {r.noise_p2p:.1f} pA peak-to-peak',loc='left',fontsize=8.5)
        ax.set_xlim(0,120); ax.set_ylabel('pA')
        if i==3: ax.set_xlabel('Time (s)')
        ax2=fig.add_subplot(gs[i,4:]); y=z[nm+'_zoom']; ax2.plot(np.arange(len(y))/20,y,lw=.4,color=col); ax2.set_title('1 s, 20 kHz',loc='left',fontsize=8.5)
        if i==3: ax2.set_xlabel('ms')
    mets=[('Ihold_median_pA','Holding current (pA)','E  Holding current, full scale',(-500,700),'F  Holding current, zoom'),
          ('drift_total_pA','Drift (pA over 2 min)','Drift over 2 min, full scale',(-120,120),'Drift, zoom'),
          ('noise_p2p','Noise (pA, peak-to-peak)','Noise, full scale',None,None)]
    for j,(c,yl,t1,zl,t2) in enumerate(mets):
        gr=[(TL[k],a[a.grp==g][c].values,COL[g]) for k,g in enumerate(GRP)]
        ax=fig.add_subplot(gs[4,j*2:(j+1)*2]); box_panel(ax,gr,full_ylim(gr,zero=(c=='noise_p2p')),yl,t1,xlabel='Command holding (mV)',zero=(c!='noise_p2p'),xtick_fs=8)
        if zl:
            ax=fig.add_subplot(gs[5,j*2:(j+1)*2]); box_panel(ax,gr,zl,yl,t2,xlabel='Command holding (mV)',zero=True,zoom=True,xtick_fs=8)
    ax=fig.add_subplot(gs[5,4:]); ax.axis('off'); legend(ax,'center',(0.4,0.5),extra=[Line2D([0],[0],color='none',label='\u2191k / \u2193k = k points beyond the axis')])
    foot(fig,'A-D: whole recording as 2 kHz averages (left) and 1 s of raw 20 kHz data from t = 20 s (right).\nE-F: one dot per recording (n = 67, 10 and 5 recordings); the lower row shows the same data on a narrower scale.',.004)
    fig.savefig(O+'K60_fig1.png',dpi=200); plt.close(fig)
# ================= Figure 2
c=pd.read_csv(K+'k60_event_freq_per_cell.csv'); n_all=len(c); c=c[~c.cellkey.isin(EXCLUDE_FIG2)].copy(); print('Figure 2: cells left out:',n_all-len(c))
gr=[(TL[k],c[c.grp==g].f.values,COL[g]) for k,g in enumerate(GRP)]
fig,axs=plt.subplots(1,2,figsize=(9.4,4.8))
box_panel(axs[0],gr,full_ylim(gr),'miniML event frequency per cell (Hz)','A  Full scale',xlabel='Command holding (mV)')
box_panel(axs[1],gr,(0,1.0),'miniML event frequency per cell (Hz)','B  Zoom to 0-1 Hz (same data)',xlabel='Command holding (mV)',zoom=True)
legend(axs[1],'upper right',extra=[Line2D([0],[0],color='none',label='↑k = k points above axis')])
foot(fig,'Each dot is one cell (mean of its recordings). Zoom panel: same data on a narrower scale. 6 Nov cell 2 (7.5 Hz) is left out of this figure.'); save(fig,O+'K60_fig2.png')
c.to_csv(K+'k60_event_freq_per_cell_fig2.csv',index=False)
print(c.groupby('grp').f.agg(['count','mean','median',lambda v:v.quantile(.25),lambda v:v.quantile(.75),'min','max']).round(3))
print('done')
