import sys; sys.path.insert(0,'/home/claude/work/style')
import numpy as np, pandas as pd
from figstyle import *
setup()
cell=pd.read_csv('noise_cs_vs_kglu_per_cell.csv'); R=pd.read_csv('noise_cs_vs_kglu_by_holding.csv')
H=[-60,-50,-40,-30,0]; rng=np.random.default_rng(5); off=.2
fig,ax=plt.subplots(figsize=(11.2,5.8))
for i,h in enumerate(H):
    for j,(lab,col) in enumerate([('K-Glu',BL),('Cs',RD)]):
        v=cell[(cell.h==h)&(cell.internal==lab)].p2p.values
        box(ax,i+(-off if j==0 else off),v,col,rng,w=.16)
        ax.text(i+(-off if j==0 else off),-0.035,f'n={len(v)}',ha='center',va='top',fontsize=8,transform=ax.get_xaxis_transform(),color='#333')
    r=R[R.h==h].iloc[0]; star='*' if r.p<.05 else ''
    ax.text(i,1.0,f'Cs/K-Glu {r.ratio:.2f}\np = {r.p:.2f}{star}',ha='center',va='bottom',fontsize=8,transform=ax.get_xaxis_transform(),color='#333')
ax.set_xticks(range(len(H))); ax.set_xticklabels([f'{h}' for h in H]); ax.tick_params(axis='x',pad=14)
ax.set_xlim(-.6,len(H)-.4); ax.set_ylim(0,nice_ceil(cell[cell.h.isin(H)].p2p.max()*1.04))
ax.set_xlabel('Command holding (mV)',labelpad=2); ax.set_ylabel('Noise per cell (pA, peak-to-peak)')
ax.set_title('Baseline noise, Cs-gluconate vs K-gluconate internal, at matched holding',loc='left',pad=30)
legend(ax,loc='upper left',anchor=(1.01,1.0),extra=[Line2D([0],[0],marker='s',color='w',mfc=BL,ms=8,label='K-gluconate'),Line2D([0],[0],marker='s',color='w',mfc=RD,ms=8,label='Cs-gluconate')])
foot(fig,'Noise = 6 x MAD-SD of the 3 kHz low-passed trace, whole 2-min recording. One dot per cell (mean of its recordings at that holding).\nTop of each group: ratio of Cs to K-Glu medians and Mann-Whitney p on cells (* p < 0.05, not corrected for 5 comparisons).\nAll holdings pooled (log noise, holding as a factor, random cell intercept): Cs / K-Glu = 0.95 (95% CI 0.81-1.12, p = 0.54).',y=.005)
save(fig,'/home/claude/work/noise/Noise_Cs_vs_KGlu.png',rect=(0,.1,1,1))
