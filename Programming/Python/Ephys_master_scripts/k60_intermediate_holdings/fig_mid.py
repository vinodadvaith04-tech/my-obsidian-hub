import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ch=pd.read_csv('mid_cell_holding.csv'); ref=pd.read_csv('../k60_event_freq_per_cell.csv')
H=[-50,-45,-40,-35,-30]
E={f'{h}':ch[ch.holding_command_mV==h].f_in.dropna().values for h in H}
I={f'{h}':ch[ch.holding_command_mV==h].f_out.dropna().values for h in H}
E['-60']=ref[ref.grp=='-60 mV'].f.values
I['0']=ref[ref.grp=='0 mV'].f.values; I['+5 to +15']=ref[ref.grp=='+5 to +15 mV'].f.values
orderE=['-60']+[str(h) for h in H]; orderI=[str(h) for h in H]+['0','+5 to +15']
YMAX=2.0
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(1,2,figsize=(11.2,4.8),gridspec_kw=dict(width_ratios=[6,7]))
rows=[]
def draw(ax,data,order,col,title,refkeys):
    rng=np.random.default_rng(5)
    for k,key in enumerate(order):
        v=np.asarray(data[key]); q1,med,q3=np.percentile(v,[25,50,75]); mn,mx,mu=v.min(),v.max(),v.mean(); w=.3
        c=col if key not in refkeys else '#7a8a99' if col=='#1f4e79' else '#b9987a'
        rows.append((title[:1],key,len(v),mu,med,q1,q3,mn,mx))
        ax.add_patch(plt.Rectangle((k-w,q1),2*w,min(q3,YMAX)-q1,fc=c,alpha=.2,ec=c,lw=1.2,zorder=1))
        ax.plot([k-w,k+w],[med,med],color=c,lw=2.4,zorder=3)
        ax.plot([k,k],[mn,q1],color=c,lw=1.2,zorder=1); ax.plot([k,k],[q3,min(mx,YMAX)],color=c,lw=1.2,zorder=1)
        ax.plot([k-w/2,k+w/2],[mn,mn],color=c,lw=1.2,zorder=1)
        if mx<=YMAX: ax.plot([k-w/2,k+w/2],[mx,mx],color=c,lw=1.2,zorder=1)
        if mu<=YMAX: ax.scatter([k],[mu],marker='D',s=40,fc='white',ec=c,lw=1.6,zorder=5)
        xs=k+rng.uniform(-.15,.15,len(v)); inr=v<=YMAX
        ax.scatter(xs[inr],v[inr],s=20,color=c,alpha=.85,lw=.4,ec='white',zorder=4)
        for x,y in zip(xs[~inr],v[~inr]):
            ax.scatter([x],[YMAX*0.985],marker='^',s=40,color=c,zorder=4); ax.text(x+.12,YMAX*0.97,f'{y:.1f} Hz',fontsize=7.5,va='top',color='#333')
        ax.text(k,-0.035,f'n={len(v)}',ha='center',va='top',fontsize=8,transform=ax.get_xaxis_transform(),color='#333')
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order); ax.tick_params(axis='x',pad=14)
    ax.set_xlim(-.6,len(order)-.4); ax.set_ylim(0,YMAX); ax.set_xlabel('Command holding (mV)',labelpad=2); ax.set_title(title,loc='left',fontsize=9.5)
    ax.set_ylabel('miniML event frequency per cell (Hz)')
draw(axs[0],E,orderE,'#1f4e79','A  Inward events (EPSC)',['-60'])
draw(axs[1],I,orderI,'#b5451b','B  Outward events (IPSC)',['0','+5 to +15'])
axs[1].legend(handles=[Line2D([0],[0],color='#555',lw=2.4,label='Median'),Line2D([0],[0],marker='D',color='w',mfc='white',mec='#555',mew=1.6,ms=6,label='Mean'),
  Line2D([0],[0],marker='s',color='w',mfc='#ccc',mec='#555',ms=8,label='Box = IQR'),Line2D([0],[0],color='#555',lw=1.2,label='Whiskers = min, max'),
  Line2D([0],[0],marker='o',color='w',mfc='#555',ms=5,label='One cell'),Line2D([0],[0],marker='^',color='w',mfc='#555',ms=6,label='Above axis')],frameon=False,fontsize=8,loc='upper left')
fig.text(.01,.005,'Greyed boxes: -60 mV (A) and 0, +5 to +15 mV (B) from the first report, for reference. Outward events at -60 mV and inward events at 0 mV were not analysed by miniML.',fontsize=7.5,color='#444')
fig.tight_layout(rect=(0,.03,1,1)); fig.savefig('fig_mid.png',dpi=200)
pd.DataFrame(rows,columns=['panel','holding','n_cells','mean','median','q1','q3','min','max']).round(3).to_csv('mid_event_freq_box_stats.csv',index=False)
print(pd.DataFrame(rows,columns=['panel','holding','n_cells','mean','median','q1','q3','min','max']).round(2))
