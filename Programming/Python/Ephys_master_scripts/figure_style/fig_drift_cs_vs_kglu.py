import sys; sys.path.insert(0,'/home/claude/work/style')
import numpy as np, pandas as pd
from figstyle import *
setup()
cell=pd.read_csv('drift_cs_vs_kglu_per_cell.csv'); R=pd.read_csv('drift_cs_vs_kglu_by_holding.csv')
H=[-60,-50,-40,-30,0]; off=.2
fig,axs=plt.subplots(1,2,figsize=(15,5.8),gridspec_kw=dict(width_ratios=[1,1]))
for pa,(ax,zoom) in enumerate(zip(axs,[None,150])):
    rng=np.random.default_rng(5)
    for i,h in enumerate(H):
        for j,(lab,col) in enumerate([('K-Glu',BL),('Cs',RD)]):
            v=cell[(cell.h==h)&(cell.internal==lab)].ad.values; x=i+(-off if j==0 else off)
            box(ax,x,v,col,rng,w=.16)
            ax.text(x,-0.035,f'n={len(v)}',ha='center',va='top',fontsize=8,transform=ax.get_xaxis_transform(),color='#333')
            if zoom and (v>zoom).sum(): ax.text(x,.985,f'↑{(v>zoom).sum()}',ha='center',va='top',fontsize=8,color=col,transform=ax.get_xaxis_transform(),fontweight='bold')
        if pa==1:
            r=R[R.h==h].iloc[0]; star='*' if r.p<.05 else ''
            ax.text(i,1.0,f'Cs/K-Glu {r.ratio:.2f}\np = {r.p:.3f}{star}',ha='center',va='bottom',fontsize=8,transform=ax.get_xaxis_transform(),color='#333')
    ax.set_xticks(range(len(H))); ax.set_xticklabels([f'{h}' for h in H]); ax.tick_params(axis='x',pad=14)
    ax.set_xlim(-.6,len(H)-.4); ax.set_ylim(0,zoom if zoom else nice_ceil(cell[cell.h.isin(H)].ad.max()*1.04))
    ax.set_xlabel('Command holding (mV)',labelpad=2); ax.set_ylabel('Absolute drift per cell (pA over 2 min)')
    ax.set_title('A  Full scale' if not zoom else 'B  Zoom to 0-150 pA (same data)',loc='left',pad=30 if zoom else 6)
legend(axs[1],loc='upper left',anchor=(1.01,1.0),extra=[Line2D([0],[0],marker='s',color='w',mfc=BL,ms=8,label='K-gluconate'),Line2D([0],[0],marker='s',color='w',mfc=RD,ms=8,label='Cs-gluconate'),Line2D([0],[0],color='none',label='↑k = k cells above the axis')])
foot(fig,'Drift = absolute change of the fitted linear trend of 1-s block medians over the 2-min recording. One dot per cell (mean of its recordings at that holding).\nTop of each group in B: ratio of Cs to K-Glu medians and Mann-Whitney p on cells (* p < 0.05, not corrected for 5 comparisons).\nAll holdings pooled (log(|drift|+1), holding as a factor, random cell intercept): Cs / K-Glu = 1.64 (95% CI 0.98-2.73, p = 0.06); without leaky or high-Ra Cs recordings 1.43 (0.87-2.35, p = 0.15).',y=.005)
save(fig,'/home/claude/work/noise/Drift_Cs_vs_KGlu.png',rect=(0,.1,1,1))
