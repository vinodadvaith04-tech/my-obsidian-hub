import pandas as pd, numpy as np, json, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.lines import Line2D
G=pd.read_csv('audit_groups_windows.csv'); R=pd.read_csv('audit_group_summary.csv')
G['cell']=G.file.str.extract(r'^(\d{8}_[^_]+_[^_]+_[^_]+_slice\d+_C\d+)')[0]
Gr=G[G.why=='random']
cs=Gr[Gr.group=='I: Cs, 0 and +10 mV']; kc=Gr[Gr.group=='I: KCl + DNQX/AP5, -60 mV']; es=Gr[Gr.group=='E: Cs, -80 and -70 mV']
res={}
res['cells']={'Cs I':cs.cell.nunique(),'KCl I':kc.cell.nunique(),'Cs E':es.cell.nunique(),'Cs I recordings':cs.file.nunique()}
def mw(a,b): return float(stats.mannwhitneyu(a,b).pvalue)
res['p_truth']=mw(cs.truth_hz,kc.truth_hz); res['p_mini']=mw(cs.mini_markers_hz,kc.mini_markers_hz)
cr=cs[cs.truth>=3]; kr=kc[kc.truth>=3]; cr=cr.assign(rec=cr.mini_tp/cr.truth); kr=kr.assign(rec=kr.mini_tp/kr.truth)
res['p_recall']=mw(cr.rec,kr.rec); res['rec_cs_med']=float(cr.rec.median()); res['rec_kc_med']=float(kr.rec.median()); res['n_rec_win']=(len(cr),len(kr))
res['adj_factor_cs']=float(cs.truth_hz.mean()/cs.mini_markers_hz.mean()); res['adj_factor_kc']=float(kc.truth_hz.mean()/kc.mini_markers_hz.mean())
res['cs_adj']=3.225*res['adj_factor_cs']; res['kc_adj']=8.95*res['adj_factor_kc']; res['gap_adj']=res['kc_adj']/res['cs_adj']
json.dump(res,open('audit_stats.json','w'),indent=1); print(json.dumps(res,indent=1))
# ---- figure 5
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
CS='#4a7c59'; CUR='#b5451b'
def draw(ax,groups,ylim,ylabel,title):
    rng=np.random.default_rng(3); lo,hi=ylim
    for k,(lab,v,col,hatch) in enumerate(groups):
        v=np.asarray(v,float); v=v[~np.isnan(v)]; q1,med,q3=np.percentile(v,[25,50,75]); w=.3
        ax.add_patch(plt.Rectangle((k-w,q1),2*w,q3-q1,fc=col,alpha=.12 if hatch else .25,ec=col,lw=1.2,ls='--' if hatch else '-'))
        ax.plot([k-w,k+w],[med,med],color=col,lw=2.4); ax.plot([k,k],[v.min(),q1],color=col,lw=1.2); ax.plot([k,k],[q3,v.max()],color=col,lw=1.2)
        ax.plot([k-w/2,k+w/2],[v.min()]*2,color=col,lw=1.2); ax.plot([k-w/2,k+w/2],[v.max()]*2,color=col,lw=1.2)
        ax.scatter([k],[v.mean()],marker='D',s=40,fc='white',ec=col,lw=1.6,zorder=5)
        ax.scatter(k+rng.uniform(-.15,.15,len(v)),v,s=20,color=col,alpha=.85,ec='white',lw=.4,zorder=4,marker='o' if not hatch else 's')
        ax.text(k,-0.035,f'n={len(v)}',ha='center',va='top',fontsize=8,transform=ax.get_xaxis_transform(),color='#333')
    ax.set_xticks(range(len(groups))); ax.set_xticklabels([g[0] for g in groups],fontsize=8); ax.tick_params(axis='x',pad=14)
    ax.set_xlim(-.6,len(groups)-.4); ax.set_ylim(lo,hi); ax.set_ylabel(ylabel); ax.set_title(title,loc='left',fontsize=9.5)
fig,axs=plt.subplots(1,2,figsize=(12,4.9),gridspec_kw=dict(width_ratios=[1.25,1]))
draw(axs[0],[('Cs, 0/+10 mV\nminiML',cs.mini_markers_hz,CS,True),('Cs, 0/+10 mV\naudited',cs.truth_hz,CS,False),('KCl + DNQX/AP5\nminiML',kc.mini_markers_hz,CUR,True),('KCl + DNQX/AP5\naudited',kc.truth_hz,CUR,False)],(0,19),'IPSC frequency per 3 s window (Hz)','A  IPSC frequency in the audited windows')
draw(axs[1],[('Cs, 0/+10 mV',cr.rec,CS,False),('KCl + DNQX/AP5',kr.rec,CUR,False)],(0,1.05),'Share of audited IPSCs found by miniML','B  miniML detection rate per window')
axs[1].axhline(1,color='#999',lw=.8,ls='--',zorder=0)
axs[0].legend(handles=[Line2D([0],[0],marker='s',color='w',mfc='#ccc',mec='#555',ms=7,label='miniML count (dashed box)'),Line2D([0],[0],marker='o',color='w',mfc='#555',ms=6,label='Audited count (solid box)'),Line2D([0],[0],marker='D',color='w',mfc='white',mec='#555',mew=1.5,ms=6,label='Mean')],frameon=False,fontsize=7.5,loc='upper left')
fig.text(.01,.005,'Each dot is one 3 s window, random windows only. Audited = miniML markers you accepted + second-detector markers + events you clicked. Panel B: windows with at least 3 audited events.',fontsize=7.3,color='#444')
fig.tight_layout(rect=(0,.04,1,1)); fig.savefig('Figure5_audit_IPSC_detection.png',dpi=200); plt.close(fig)
