"""Cs-gluconate report additions (26 Nov 2026): per-cell lines across holdings (Figures 4, 5) and the -80 vs -70 / 0 vs +10 comparisons (Figure 6).
Unit = cell (mean of its recordings at a holding). Five recordings with holding current beyond +-1 nA (cell lost) are left out, as in Figures 1-2."""
import sys; sys.path.insert(0,'/home/claude/work/style')
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
from figstyle import *
setup(); C='/home/claude/work/cs/'; O='/home/claude/work/cs_new/'
a=pd.read_csv(C+'cs_analysis_set.csv')
LOST=(a[['Ihold_median_pA','Ihold_first10s_pA','Ihold_last10s_pA']].abs().max(axis=1)>1000); a=a[~LOST].copy()
a['p2p']=6*a.noise_lp3k_pA; a['aI']=a.Ihold_median_pA.abs()
H=[-80,-70,-60,-50,-40,-30,0,10]; xl='Command holding (mV)'
cm=lambda c: a.groupby(['cellkey','holding_command_mV'])[c].mean().unstack()
def stars(p): return '***' if p<.001 else '**' if p<.01 else '*' if p<.05 else 'ns'
def wp(x,y):
    d=pd.concat([x,y],axis=1).dropna()
    if len(d)<5 or (d.iloc[:,0]-d.iloc[:,1]).abs().sum()==0: return np.nan,len(d)
    return wilcoxon(d.iloc[:,0],d.iloc[:,1]).pvalue,len(d)
def lines(ax,df,col,ylim,title,ylab,zoom=False,ref=None,hline=None,stars_on=True):
    """one grey line per cell through the holdings it has; bold median line; stars between neighbouring holdings"""
    lo,hi=ylim; xs={h:i for i,h in enumerate(H)}
    for ck,r in df.iterrows():
        r=r.reindex(H).dropna()
        ax.plot([xs[h] for h in r.index],r.values,color='#9aa3ad',lw=.8,alpha=.8,zorder=1,marker='o',ms=2.6,mec='none')
    med=df.reindex(columns=H).median()
    ax.plot([xs[h] for h in med.index],med.values,color=col,lw=2.8,marker='D',ms=6,mfc='white',mec=col,mew=1.6,zorder=5)
    for h in H:
        n=int(df[h].notna().sum()); ax.text(xs[h],-0.035,f'n={n}',ha='center',va='top',fontsize=8,transform=ax.get_xaxis_transform(),color='#333')
        if zoom:
            up=int((df[h]>hi).sum())
            if up: ax.text(xs[h],.985,f'↑{up}',ha='center',va='top',fontsize=8,color=col,transform=ax.get_xaxis_transform(),fontweight='bold')
    if stars_on:   # paired Wilcoxon between neighbouring holdings (cells with both), or against the reference holding
        for i in range(len(H)-1):
            if ref is None: p,n=wp(df[H[i]],df[H[i+1]]); x=i+.5
            else: continue
            if np.isnan(p): continue
            ax.text(x,1.012,stars(p),ha='center',va='bottom',fontsize=9,color='#222' if p<.05 else '#888',transform=ax.get_xaxis_transform())
    if ref is not None and stars_on:
        for h in H:
            if h==ref: continue
            p,n=wp(df[ref],df[h])
            if np.isnan(p): continue
            ax.text(xs[h],1.012,stars(p),ha='center',va='bottom',fontsize=9,color='#222' if p<.05 else '#888',transform=ax.get_xaxis_transform())
    ax.set_xticks(range(len(H))); ax.set_xticklabels([str(h) for h in H]); ax.tick_params(axis='x',pad=14)
    ax.set_xlim(-.5,len(H)-.5); ax.set_ylim(lo,hi); ax.set_xlabel(xl,labelpad=2); ax.set_ylabel(ylab); ax.set_title(title,loc='left',pad=16)
    if hline is not None: ax.axhline(hline,color='#999',lw=.9,ls='--',zorder=0)
LEG=[Line2D([0],[0],color='#9aa3ad',lw=1,marker='o',ms=3,label='One cell (line through its holdings)'),Line2D([0],[0],color='#555',lw=2.8,marker='D',mfc='white',ms=6,label='Median of cells'),
     Line2D([0],[0],color='none',label='↑k = k cells above the zoom axis')]
fi,fo=cm('f_in'),cm('f_out')
# ---------- Figure 4: raw frequency, every cell
fig,axs=plt.subplots(2,2,figsize=(13,9.4))
lines(axs[0,0],fi,BL,(0,3.3),'A  Inward events (EPSC), full scale','miniML frequency per cell (Hz)')
lines(axs[0,1],fo,RD,(0,11.8),'B  Outward events (IPSC), full scale','miniML frequency per cell (Hz)')
lines(axs[1,0],fi,BL,(0,.6),'C  Inward events, zoom to 0-0.6 Hz','miniML frequency per cell (Hz)',zoom=True,stars_on=False)
lines(axs[1,1],fo,RD,(0,2.0),'D  Outward events, zoom to 0-2 Hz','miniML frequency per cell (Hz)',zoom=True,stars_on=False)
axs[0,1].legend(handles=LEG,frameon=False,loc='upper left',fontsize=7.5)
foot(fig,'Each grey line is one cell; a holding where a cell was not recorded is skipped and the neighbouring points are joined. Bold line = median of the cells at each holding. -20 and +15 mV (1 cell each) are not drawn.\nStars (A, B) = paired Wilcoxon signed-rank between neighbouring holdings, cells with both only (* p < 0.05, ** p < 0.01, *** p < 0.001, ns = not significant; no correction for multiple comparisons).\nLower row: same data as the upper row on a narrower scale. 5 recordings with a holding current beyond 1 nA (cell lost) are left out.',.004)
save(fig,O+'Cs_Figure4.png',rect=(0,.075,1,1))
# ---------- Figure 5: each cell relative to its own reference holding
rin=fi.div(fi[-80],axis=0)*100; rout=fo.div(fo[0],axis=0)*100
rin=rin[fi[-80]>0]; rout=rout[fo[0]>0]
fig,axs=plt.subplots(2,2,figsize=(14,10.5))
lines(axs[0,0],rin,BL,(0,max(250,np.nanmax(rin.values)*1.03)),'A  Inward frequency / its own value at -80 mV','Share of the -80 mV frequency (%)',ref=-80,hline=100)
lines(axs[0,1],rout,RD,(0,max(120,np.nanpercentile(rout.values,99.5)*1.05)),'B  Outward frequency / its own value at 0 mV','Share of the 0 mV frequency (%)',ref=0,hline=100)
lines(axs[1,0],rin,BL,(0,200),'C  Inward share, zoom to 0-200%','Share of the -80 mV frequency (%)',zoom=True,ref=-80,hline=100,stars_on=False)
lines(axs[1,1],rout,RD,(0,150),'D  Outward share, zoom to 0-150%','Share of the 0 mV frequency (%)',zoom=True,ref=0,hline=100,stars_on=False)
axs[0,0].legend(handles=LEG[:2],frameon=False,loc='upper right',fontsize=7.5)
foot(fig,'Each grey line is one cell, divided by its own reference holding (inward at -80 mV, outward at 0 mV) so cells of very different frequency sit on one scale; skipped holdings are joined. Bold = median.\nStars (A, B) = paired Wilcoxon signed-rank of each holding against the reference holding, cells with both only (* p < 0.05, ** p < 0.01, *** p < 0.001, ns; no correction for multiple comparisons).\nLower row: same data on a narrower scale. 5 recordings with a holding current beyond 1 nA (cell lost) are left out.',.004)
save(fig,O+'Cs_Figure5.png',rect=(0,.075,1,1))
# ---------- Figure 6: -80 vs -70, 0 vs +10 (paired cells)
rows=[(-80,-70,[('f_in','Inward frequency (EPSC)\n(Hz)',BL),('f_out','Outward frequency (IPSC left over)\n(Hz)',RD),('aI','Holding current, absolute\n(pA)',GR),('p2p','Noise, peak-to-peak\n(pA)',GY)],'EPSC isolation'),
      (0,10,[('f_out','Outward frequency (IPSC)\n(Hz)',RD),('f_in','Inward frequency (EPSC left over)\n(Hz)',BL),('aI','Holding current, absolute\n(pA)',GR),('p2p','Noise, peak-to-peak\n(pA)',GY)],'IPSC isolation')]
fig,axs=plt.subplots(2,4,figsize=(15,8.6)); out=[]
for r,(h1,h2,pan,name) in enumerate(rows):
    for c,(col,lab,colr) in enumerate(pan):
        ax=axs[r,c]; df=cm(col)[[h1,h2]].dropna(); p,n=wp(df[h1],df[h2]); rng=np.random.default_rng(3)
        for i,(h,sub) in enumerate(zip([h1,h2],[df[h1],df[h2]])): box(ax,i,sub.values,colr,rng,w=.22,dots=False)
        for _,rw in df.iterrows(): ax.plot([0,1],[rw[h1],rw[h2]],color='#9aa3ad',lw=.9,marker='o',ms=3,mec='none',zorder=4)
        top=df.values.max(); ax.set_ylim(0,top*1.22 if top>0 else 1); ytxt=top*1.1
        ax.plot([0,1],[ytxt,ytxt],color='#333',lw=1); ax.text(.5,ytxt*1.01,stars(p),ha='center',va='bottom',fontsize=11)
        ax.set_xticks([0,1]); ax.set_xticklabels([f'{h1}',f'{h2}' if h2!=10 else '+10']); ax.set_xlim(-.6,1.6); ax.set_ylabel(lab,fontsize=8.5)
        ax.text(.5,-0.08,f'n={n} cells, p = {p:.3g}',ha='center',va='top',fontsize=8,transform=ax.transAxes,color='#333'); ax.set_xlabel(xl,labelpad=14)
        ax.set_title(('ABCD'[c] if r==0 else 'EFGH'[c])+f'  {name}'+(': target' if c==0 else (': contamination' if col in('f_in','f_out') else '')) if c<2 else ('ABCD'[c] if r==0 else 'EFGH'[c])+'  '+lab.split('\n')[0],loc='left',fontsize=9)
        out.append((h1,h2,col,n,df[h1].median(),df[h2].median(),p))
foot(fig,'Each grey line is one cell recorded at both holdings (paired); boxes = IQR, bar = median, diamond = mean, whiskers = min, max. Stars from a paired Wilcoxon signed-rank test:\n* p < 0.05, ** p < 0.01, *** p < 0.001, ns = not significant (no correction for multiple comparisons). Top row: -80 vs -70 mV (EPSC isolation). Bottom row: 0 vs +10 mV (IPSC isolation). 5 lost-cell recordings are left out.',.004)
save(fig,O+'Cs_Figure6.png',rect=(0,.07,1,1))
pd.DataFrame(out,columns=['h1','h2','measure','n','median_h1','median_h2','p']).to_csv(O+'cs_fig6_stats.csv',index=False)
print(pd.DataFrame(out,columns=['h1','h2','measure','n','median_h1','median_h2','p']).round(4).to_string())
