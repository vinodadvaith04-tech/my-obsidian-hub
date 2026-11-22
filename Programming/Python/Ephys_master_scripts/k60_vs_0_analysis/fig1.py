import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
a=pd.read_csv('k60_analysis_set3.csv'); cs=pd.read_csv('k60_gapfree_features_cs60.csv')
z=np.load('k60_example_traces.npz')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(7.4,10.4)); gs=GridSpec(5,6,figure=fig,hspace=0.85,wspace=1.2,height_ratios=[1,1,1,1,1.15])
ex=[('A  -60 mV, typical','20260929_Adv00_NT_unk_slice1_C5_039.abf','#1f4e79'),
    ('B  -60 mV, unstable cell','20261007_Adv00_NT_unk_slice1_C5_048.abf','#1f4e79'),
    ('C  0 mV, drifting and noisy','20261002_Adv00_NT_unk_slice1_C1_004.abf','#b5451b'),
    ('D  0 mV, spike-like transients','20260925_Adv00_NT_unk_slice1_C5_067.abf','#b5451b'),
    ]
for i,(lab,nm,col) in enumerate(ex):
    ax=fig.add_subplot(gs[i,:4]); y=z[nm+'_dec']; t=np.arange(len(y))/2000
    r=a[a.new_name==nm].iloc[0]; ax.plot(t,y,lw=.4,color=col); ax.set_title(lab+'  ('+nm[:8]+', slice '+nm.split('_')[4]+' '+nm.split('_')[5]+')'+f'\nholding current {r.Ihold_median_pA:.0f} pA, drift {r.drift_total_pA:+.0f} pA over 2 min, noise {r.noise_lp3k_pA:.1f} pA',loc='left',fontsize=7.5)


    ax.set_xlim(0,120); ax.set_ylabel('pA'); 
    if i==3: ax.set_xlabel('Time (s)')
    ax2=fig.add_subplot(gs[i,4:]); y=z[nm+'_zoom']; ax2.plot(np.arange(len(y))/20,y,lw=.4,color=col); ax2.set_title('1 s, 20 kHz',loc='left',fontsize=7)
    ax2.set_xlabel('ms') if i==3 else None
# 4th row strip plots
grp=['-60 mV','0 mV','+5 to +15 mV']; cols={'-60 mV':'#1f4e79','0 mV':'#b5451b','+5 to +15 mV':'#c9892b'}
mets=[('Ihold_median_pA','Holding current (pA)'),('drift_total_pA','Drift over 2 min (pA)'),('noise_lp3k_pA','Noise (pA, MAD-SD)')]
rng=np.random.default_rng(1)
for j,(c,lab) in enumerate(mets):
    ax=fig.add_subplot(gs[4,j*2:(j+1)*2])
    for k,gname in enumerate(grp):
        d=a[a.grp==gname][c].values; x=k+rng.uniform(-.15,.15,len(d))
        ax.scatter(x,d,s=7,color=cols[gname],alpha=.7,lw=0)
        ax.hlines(np.median(d),k-.3,k+.3,color='k',lw=1.2)
    ax.set_xticks(range(3)); ax.set_xticklabels(['-60','0','+5…15'],fontsize=7); ax.set_xlabel('Command holding (mV)',fontsize=7); ax.set_title(('E  ' if j==0 else '')+lab,fontsize=8,loc='left')
    ax.axhline(0,color='#999',lw=.5) if c!='noise_lp3k_pA' else None
    if c=='Ihold_median_pA': ax.set_ylim(-1050,1300)
    if c=='drift_total_pA': ax.set_ylim(-800,450)
fig.text(0.01,0.005,'A-D: whole recording as 2 kHz averages (left) and 1 s of raw 20 kHz data from t = 20 s (right). E: one dot per recording, bar = median; n = 67, 10 and 5 recordings.',fontsize=6.5,color='#444')
plt.savefig('fig1.png',dpi=170,bbox_inches='tight'); print('ok')
