"""Shared figure style for all four ephys reports (K-Glu -60/0, K-Glu intermediate, Cs-Glu, Lab standards).
Rules: linear axes only; one box-plot look (box = IQR, thick line = median, diamond = mean, whiskers = min-max,
one dot per cell/recording, n under each box); outliers are never clipped in the full-scale panel; a zoom panel
(same data, narrower y-range) sits next to/under it; dots beyond a zoom axis are counted as 'up k' text, not triangles."""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
BL='#1f4e79'; RD='#b5451b'; GR='#4a7c59'; PU='#5b4a8a'; OC='#c9892b'; GY='#8a8a8a'; GREYBOX='#7b8794'
FOOT=7.5
def setup():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.titlesize':9.5,'axes.labelsize':9,'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':7.5})
def nice_ceil(x):
    if x<=0: return 1
    e=10**np.floor(np.log10(x)); m=x/e
    for s in (1,1.2,1.5,2,2.5,3,4,5,6,8,10):
        if m<=s+1e-9: return s*e
def full_ylim(groups,zero=True,pad=.06):
    v=np.concatenate([np.asarray(g[1],float)[~np.isnan(np.asarray(g[1],float))] for g in groups])
    lo,hi=v.min(),v.max(); span=hi-lo if hi>lo else abs(hi) or 1
    lo=0 if (zero and lo>=0) else lo-pad*span
    return lo,hi+pad*span
def box(ax,k,v,col,rng,w=.3,alpha=.2,zorder=1,dots=True,ls='-',marker='o'):
    v=np.asarray(v,float); v=v[~np.isnan(v)]
    if len(v)==0: return None
    q1,med,q3=np.percentile(v,[25,50,75]); mn,mx,mu=v.min(),v.max(),v.mean()
    ax.add_patch(plt.Rectangle((k-w,q1),2*w,q3-q1,fc=col,alpha=alpha,ec=col,lw=1.2,ls=ls,zorder=zorder))
    ax.plot([k-w,k+w],[med,med],color=col,lw=2.4,zorder=zorder+2)
    ax.plot([k,k],[mn,q1],color=col,lw=1.2,zorder=zorder); ax.plot([k,k],[q3,mx],color=col,lw=1.2,zorder=zorder)
    for y in (mn,mx): ax.plot([k-w/2,k+w/2],[y,y],color=col,lw=1.2,zorder=zorder)
    ax.scatter([k],[mu],marker='D',s=40,fc='white',ec=col,lw=1.6,zorder=zorder+4)
    if dots: ax.scatter(k+rng.uniform(-.15,.15,len(v)),v,s=20,color=col,alpha=.85,lw=.4,ec='white',zorder=zorder+3,marker=marker)
    return (len(v),mu,med,q1,q3,mn,mx)
def box_panel(ax,groups,ylim,ylabel,title,xlabel=None,hline=None,zero=False,band=None,stats=None,tag=None,
              zoom=False,xtick_fs=9,xpad=14,seed=11):
    """groups: list of (xticklabel, values, colour[, alpha])  -> one box per group."""
    rng=np.random.default_rng(seed); lo,hi=ylim
    for k,g in enumerate(groups):
        lab,v,col=g[:3]; al=g[3] if len(g)>3 else .2; ls=g[4] if len(g)>4 else '-'; mk=g[5] if len(g)>5 else 'o'
        v=np.asarray(v,float); v=v[~np.isnan(v)]
        r=box(ax,k,v,col,rng,alpha=al,ls=ls,marker=mk)
        if r is None: continue
        if stats is not None: stats.append((tag,str(lab).replace('\n',' '),*r))
        ax.text(k,-0.035,f'n={len(v)}',ha='center',va='top',fontsize=8,transform=ax.get_xaxis_transform(),color='#333')
        if zoom:
            up=int((v>hi).sum()); dn=int((v<lo).sum())
            if up: ax.text(k,.985,f'↑{up}',ha='center',va='top',fontsize=8,color=col,transform=ax.get_xaxis_transform(),fontweight='bold')
            if dn: ax.text(k,.015,f'↓{dn}',ha='center',va='bottom',fontsize=8,color=col,transform=ax.get_xaxis_transform(),fontweight='bold')
    ax.set_xticks(range(len(groups))); ax.set_xticklabels([g[0] for g in groups],fontsize=xtick_fs); ax.tick_params(axis='x',pad=xpad)
    ax.set_xlim(-.6,len(groups)-.4); ax.set_ylim(lo,hi)
    if xlabel: ax.set_xlabel(xlabel,labelpad=2)
    ax.set_ylabel(ylabel); ax.set_title(title,loc='left')
    if zero: ax.axhline(0,color='#999',lw=.7,zorder=0)
    if hline is not None: ax.axhline(hline,color='#999',lw=.9,ls='--',zorder=0)
    if band is not None: ax.axhspan(band[0],band[1],color=band[2] if len(band)>2 else '#888',alpha=.08,zorder=0,lw=0)
def legend(ax,loc='upper left',anchor=None,extra=None):
    h=[Line2D([0],[0],color='#555',lw=2.4,label='Median'),Line2D([0],[0],marker='D',color='w',mfc='white',mec='#555',mew=1.6,ms=6,label='Mean'),
       Line2D([0],[0],marker='s',color='w',mfc='#ccc',mec='#555',ms=8,label='Box = IQR'),Line2D([0],[0],color='#555',lw=1.2,label='Whiskers = min, max'),
       Line2D([0],[0],marker='o',color='w',mfc='#555',ms=5,label='One cell')]
    if extra: h+=extra
    ax.legend(handles=h,frameon=False,loc=loc,bbox_to_anchor=anchor)
ZOOMNOTE='In zoom panels, ↑k = k points above the axis (all are shown in the full-scale panel).'
def foot(fig,text,y=.005): fig.text(.01,y,text,fontsize=FOOT,color='#444')
def save(fig,path,rect=(0,.04,1,1)):
    fig.tight_layout(rect=rect); fig.savefig(path,dpi=200); plt.close(fig)
