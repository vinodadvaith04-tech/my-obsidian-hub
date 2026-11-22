import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
g=pd.read_csv('mid_master_table.csv'); a=pd.read_csv('../k60_analysis_set3.csv'); a=a[a.grp=='-60 mV']
H=[-50,-45,-40,-35,-30]
def cell(d,h=None):
    return d.groupby('cellkey').agg(noise=('noise_lp3k_pA','mean'),drift=('drift_total_pA','mean'))
data={'-60':cell(a)}
for h in H: data[str(h)]=cell(g[g.holding_command_mV==h])
order=['-60']+[str(h) for h in H]
rows=[]
for k in order:
    for m in ['noise','drift']:
        v=data[k][m].values; q1,med,q3=np.percentile(v,[25,50,75]); rows.append((m,k,len(v),v.mean(),med,q1,q3,v.min(),v.max()))
st=pd.DataFrame(rows,columns=['measure','holding','n_cells','mean','median','q1','q3','min','max']).round(2); st.to_csv('mid_noise_drift_box_stats.csv',index=False)
cell_rows=[]
for k in order:
    for ck,r in data[k].iterrows(): cell_rows.append((k,ck,round(r.noise,3),round(r.drift,2)))
pd.DataFrame(cell_rows,columns=['holding','cellkey','noise_pA','drift_pA']).to_csv('mid_noise_drift_per_cell.csv',index=False)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(1,3,figsize=(13.2,4.8),gridspec_kw=dict(width_ratios=[1,1,1]))
def draw(ax,m,ylim,title,ylab,zero=False):
    rng=np.random.default_rng(7); lo,hi=ylim
    for k,key in enumerate(order):
        v=data[key][m].values; q1,med,q3=np.percentile(v,[25,50,75]); mn,mx,mu=v.min(),v.max(),v.mean(); w=.3
        c='#7a8a99' if key=='-60' else '#1f4e79' if m=='noise' else '#b5451b'
        ax.add_patch(plt.Rectangle((k-w,max(q1,lo)),2*w,min(q3,hi)-max(q1,lo),fc=c,alpha=.2,ec=c,lw=1.2,zorder=1))
        if lo<=med<=hi: ax.plot([k-w,k+w],[med,med],color=c,lw=2.4,zorder=3)
        ax.plot([k,k],[max(mn,lo),max(q1,lo)],color=c,lw=1.2,zorder=1); ax.plot([k,k],[min(q3,hi),min(mx,hi)],color=c,lw=1.2,zorder=1)
        if lo<=mn<=hi: ax.plot([k-w/2,k+w/2],[mn,mn],color=c,lw=1.2,zorder=1)
        if lo<=mx<=hi: ax.plot([k-w/2,k+w/2],[mx,mx],color=c,lw=1.2,zorder=1)
        if lo<=mu<=hi: ax.scatter([k],[mu],marker='D',s=40,fc='white',ec=c,lw=1.6,zorder=5)
        xs=k+rng.uniform(-.15,.15,len(v)); ins=(v>=lo)&(v<=hi)
        ax.scatter(xs[ins],v[ins],s=20,color=c,alpha=.85,lw=.4,ec='white',zorder=4)
        for x,y in zip(xs[~ins],v[~ins]):
            up=y>hi; ax.scatter([x],[hi-(hi-lo)*.015 if up else lo+(hi-lo)*.015],marker='^' if up else 'v',s=34,color=c,zorder=4)
        ax.text(k,-0.035,f'n={len(v)}',ha='center',va='top',fontsize=8,transform=ax.get_xaxis_transform(),color='#333')
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order); ax.tick_params(axis='x',pad=14)
    ax.set_xlim(-.6,len(order)-.4); ax.set_ylim(lo,hi); ax.set_xlabel('Command holding (mV)',labelpad=2); ax.set_ylabel(ylab); ax.set_title(title,loc='left',fontsize=9.5)
    if zero: ax.axhline(0,color='#999',lw=.7,zorder=0)
draw(axs[0],'noise',(0,9),'A  Noise (MAD-SD, 3 kHz low-pass)','Noise per cell (pA)')
draw(axs[1],'drift',(-800,300),'B  Drift over 2 min, full range','Drift per cell (pA over the recording)',True)
draw(axs[2],'drift',(-60,60),'C  Drift, zoom to ±60 pA','Drift per cell (pA over the recording)',True)
axs[1].legend(handles=[Line2D([0],[0],color='#555',lw=2.4,label='Median'),Line2D([0],[0],marker='D',color='w',mfc='white',mec='#555',mew=1.6,ms=6,label='Mean'),
  Line2D([0],[0],marker='s',color='w',mfc='#ccc',mec='#555',ms=8,label='Box = IQR'),Line2D([0],[0],color='#555',lw=1.2,label='Whiskers = min, max'),
  Line2D([0],[0],marker='o',color='w',mfc='#555',ms=5,label='One cell'),Line2D([0],[0],marker='^',color='w',mfc='#555',ms=6,label='Outside axis')],frameon=False,fontsize=7.5,loc='lower left',bbox_to_anchor=(0.0,0.02))
fig.text(.01,.005,'Grey = -60 mV cells from the first report. Each dot is one cell (mean of its recordings at that holding). Drift is signed: negative = current falling over the 2 min.',fontsize=7.5,color='#444')
fig.tight_layout(rect=(0,.03,1,1)); fig.savefig('fig_mid2.png',dpi=200)
print(st.to_string())
