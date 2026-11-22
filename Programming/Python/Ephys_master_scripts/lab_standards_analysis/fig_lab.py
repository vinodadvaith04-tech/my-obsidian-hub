import numpy as np, pandas as pd, json, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
D=pd.read_csv('lab_cells_long.csv'); S=json.load(open('lab_stats.json'))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
PAST='#5b4a8a'; CUR='#b5451b'; CS='#4a7c59'
BOX=[]
def draw(ax,groups,ylim,ylabel,title,tag,log=False):
    """groups: list of (label, values, colour)"""
    rng=np.random.default_rng(11); lo,hi=ylim
    for k,(lab,v,col) in enumerate(groups):
        v=np.asarray(v,float); v=v[~np.isnan(v)]
        q1,med,q3=np.percentile(v,[25,50,75]); mn,mx,mu=v.min(),v.max(),v.mean(); w=.3
        BOX.append((tag,lab.replace('\n',' '),len(v),mu,med,q1,q3,mn,mx))
        cl=lambda y: min(max(y,lo),hi)
        ax.add_patch(plt.Rectangle((k-w,cl(q1)),2*w,max(cl(q3)-cl(q1),0),fc=col,alpha=.2,ec=col,lw=1.2,zorder=1))
        ax.plot([k-w,k+w],[med,med],color=col,lw=2.4,zorder=3)
        ax.plot([k,k],[cl(mn),cl(q1)],color=col,lw=1.2,zorder=1); ax.plot([k,k],[cl(q3),cl(mx)],color=col,lw=1.2,zorder=1)
        ax.plot([k-w/2,k+w/2],[mn,mn],color=col,lw=1.2,zorder=1); ax.plot([k-w/2,k+w/2],[mx,mx],color=col,lw=1.2,zorder=1)
        ax.scatter([k],[mu],marker='D',s=40,fc='white',ec=col,lw=1.6,zorder=5)
        xs=k+rng.uniform(-.15,.15,len(v)); ax.scatter(xs,v,s=20,color=col,alpha=.85,lw=.4,ec='white',zorder=4)
        ax.text(k,-0.035,f'n={len(v)}',ha='center',va='top',fontsize=8,transform=ax.get_xaxis_transform(),color='#333')
    ax.set_xticks(range(len(groups))); ax.set_xticklabels([g[0] for g in groups],fontsize=8); ax.tick_params(axis='x',pad=14)
    ax.set_xlim(-.6,len(groups)-.4); ax.set_ylim(lo,hi)
    if log: ax.set_yscale('log')
    ax.set_ylabel(ylabel); ax.set_title(title,loc='left',fontsize=9.5)
def legend(ax,loc='upper left',anchor=None):
    ax.legend(handles=[Line2D([0],[0],color='#555',lw=2.4,label='Median'),Line2D([0],[0],marker='D',color='w',mfc='white',mec='#555',mew=1.6,ms=6,label='Mean'),
      Line2D([0],[0],marker='s',color='w',mfc='#ccc',mec='#555',ms=8,label='Box = IQR'),Line2D([0],[0],color='#555',lw=1.2,label='Whiskers = min, max'),
      Line2D([0],[0],marker='o',color='w',mfc='#555',ms=5,label='One cell')],frameon=False,fontsize=7.5,loc=loc,bbox_to_anchor=anchor)
def V(src,meas): return D[(D.source==src)&(D.measure==meas)].freq.values
def VV(srcs,meas): return D[(D.source.isin(srcs))&(D.measure==meas)].freq.values
SAM="Sam's processed EPSC, control cells (his pipeline)"; AE='Advaith EPSC, Sam experiment (miniML)'; HI='Helena KCl IPSC (miniML)'; AI='Advaith KCl IPSC, Sam experiment (miniML)'
# ---------- Figure 1: three-way comparison
fig,axs=plt.subplots(1,2,figsize=(11.6,4.9))
draw(axs[0],[('Past lab\n(paper, DIV9)',V('Paper DIV9 control','E'),PAST),('Past lab\n(paper, DIV21)',V('Paper DIV21 control','E'),PAST),('Current lab\n(Sam + Advaith,\nK-Glu, -60 mV)',VV([SAM,AE],'E'),CUR),('Advaith Cs-Glu\n(-80 mV)',V('Cs-gluconate, E at -80 mV (miniML)','E'),CS)],(0,3.3),'EPSC frequency per cell (Hz)','A  EPSC frequency','E1')
draw(axs[1],[('Past lab\n(paper, DIV9)',V('Paper DIV9 control','I'),PAST),('Past lab\n(paper, DIV21)',V('Paper DIV21 control','I'),PAST),('Current lab\n(Helena + Advaith,\nKCl, -60 mV)',VV([HI,AI],'I'),CUR),('Advaith Cs-Glu\n(0 mV)',V('Cs-gluconate, I at 0 mV (miniML)','I'),CS)],(0,15),'IPSC frequency per cell (Hz)','B  IPSC frequency','I1')
legend(axs[1],'upper left',(0.0,1.0))
fig.text(.01,.005,'Each dot is one cell. Paper values are read off the published plots (control cells only). Different detectors, internals and holdings between groups: see caveats.',fontsize=7.5,color='#444')
fig.tight_layout(rect=(0,.04,1,1)); fig.savefig('Figure1_lab_standards_frequency.png',dpi=200); plt.close(fig)
# ---------- Figure 2: by source
fig,axs=plt.subplots(1,2,figsize=(12.6,4.9))
draw(axs[0],[('Paper\nDIV9',V('Paper DIV9 control','E'),PAST),('Paper\nDIV21',V('Paper DIV21 control','E'),PAST),("Sam's own\nprocessed data",V(SAM,'E'),CUR),("Advaith,\nSam's experiment",V(AE,'E'),CUR),('Advaith\nCs-Glu, -80',V('Cs-gluconate, E at -80 mV (miniML)','E'),CS)],(0,3.3),'EPSC frequency per cell (Hz)','A  EPSC frequency by data source','E2')
draw(axs[1],[('Paper\nDIV9',V('Paper DIV9 control','I'),PAST),('Paper\nDIV21',V('Paper DIV21 control','I'),PAST),("Helena\nKCl",V(HI,'I'),CUR),("Advaith KCl,\nSam's experiment",V(AI,'I'),CUR),('Advaith\nCs-Glu, 0',V('Cs-gluconate, I at 0 mV (miniML)','I'),CS)],(0,15),'IPSC frequency per cell (Hz)','B  IPSC frequency by data source','I2')
legend(axs[1],'upper left',(0.0,1.0))
fig.text(.01,.005,'Same data as Figure 1 with the current-lab sources shown separately. Each dot is one cell.',fontsize=7.5,color='#444')
fig.tight_layout(rect=(0,.04,1,1)); fig.savefig('Figure2_lab_standards_by_source.png',dpi=200); plt.close(fig)
# ---------- Figure 3: E/I individually measured
R=S['ratio']; cw=pd.read_csv('cs_E_I_same_cell.csv')
fig,ax=plt.subplots(figsize=(8.6,5.0))
items=[('Past lab\n(paper, DIV9)','Paper DIV9',PAST),('Past lab\n(paper, DIV21)','Paper DIV21',PAST),('Current lab\n(Sam + Advaith E,\nHelena + Advaith I)','Current lab (pooled)',CUR),('Advaith Cs-Glu\n(E at -80, I at 0,\nseparate cells)','Cs-gluconate (-80 / 0 mV, separate cells)',CS)]
for k,(lab,key,col) in enumerate(items):
    r=R[key]; ax.errorbar([k],[r['ratio_mean']],yerr=[[r['ratio_mean']-r['ci_mean'][0]],[r['ci_mean'][1]-r['ratio_mean']]],fmt='D',ms=8,color=col,mfc='white',mew=1.8,capsize=5,lw=1.6,zorder=4)
    ax.text(k+.12,r['ratio_mean'],f"{r['ratio_mean']:.2f}",va='center',fontsize=8.5,color='#333')
k=4; rng=np.random.default_rng(5); v=cw.ratio.values; q1,med,q3=np.percentile(v,[25,50,75]); w=.3
ax.add_patch(plt.Rectangle((k-w,q1),2*w,q3-q1,fc=CS,alpha=.2,ec=CS,lw=1.2)); ax.plot([k-w,k+w],[med,med],color=CS,lw=2.4)
ax.plot([k,k],[v.min(),q1],color=CS,lw=1.2); ax.plot([k,k],[q3,v.max()],color=CS,lw=1.2)
ax.scatter(k+rng.uniform(-.15,.15,len(v)),v,s=22,color=CS,alpha=.85,ec='white',lw=.4,zorder=4)
ax.text(k+.36,med,f'{med:.2f}',va='center',fontsize=8.5,color='#333')
ax.set_ylim(0,1.6); ax.set_yticks([0,.2,.4,.6,.8,1.0,1.2,1.4,1.6])
ax.set_xticks(range(5)); ax.set_xticklabels([i[0] for i in items]+['Advaith Cs-Glu\n(same cell,\nper cell)'],fontsize=8); ax.tick_params(axis='x',pad=4)
ax.set_xlim(-.6,4.7); ax.set_ylabel('E/I frequency ratio (EPSC Hz / IPSC Hz)'); ax.axhline(1,color='#999',lw=.8,ls='--',zorder=0)
ax.set_title('E/I ratio from separately measured EPSCs and IPSCs',loc='left',fontsize=9.5)
ax.axvline(3.5,color='#ccc',lw=.8)
fig.text(.01,.005,'Diamonds: ratio of group mean frequencies with 95% bootstrap interval (E and I come from different cells).\nRight: Cs, ratio within each cell (n=14): box = IQR, bar = median, dots = cells.',fontsize=7.3,color='#444')
fig.tight_layout(rect=(0,.07,1,1)); fig.savefig('Figure3_lab_standards_EI_individual.png',dpi=200); plt.close(fig)
# ---------- Figure 4: intermediate holding (Cs only)
ei=pd.read_csv('/home/claude/work/cs/cs_ei_per_cell_ref-80_0.csv')
fig,axs=plt.subplots(1,2,figsize=(11.6,4.9),gridspec_kw=dict(width_ratios=[1,1.15]))
t=ei[ei.h==-40]; ax=axs[0]
for _,r in t.iterrows(): ax.plot([0,1],[r.ratio_ref,r.ratio_h],color='#999',lw=.9,alpha=.8,zorder=1)
ax.scatter(np.zeros(len(t)),t.ratio_ref,s=28,color=CS,ec='white',lw=.4,zorder=3); ax.scatter(np.ones(len(t)),t.ratio_h,s=28,color=CS,ec='white',lw=.4,zorder=3)
for x,c in ((0,'ratio_ref'),(1,'ratio_h')): ax.plot([x-.18,x+.18],[t[c].median()]*2,color='k',lw=2.4,zorder=4)
ax.set_ylim(0,3); ax.axhline(1,color='#999',lw=.8,ls='--',zorder=0)
ax.set_xticks([0,1]); ax.set_xticklabels(['Measured separately\n(E at -80, I at 0)','Measured together\nat -40 mV']); ax.set_xlim(-.5,1.5)
ax.set_ylabel('E/I frequency ratio per cell'); ax.set_title(f'A  Same cells, two ways (n={len(t)})',loc='left',fontsize=9.5)
ax.text(.5,.04,'line = one cell; black bar = median',transform=ax.transAxes,ha='center',fontsize=7.5,color='#444')
mid=[-50,-40,-30]
draw(axs[1],[(f'{h} mV',ei[ei.h==h].fold.values,CS) for h in mid],(0,20),'Ratio at holding / separately measured ratio (fold)','B  How far each holding moves the ratio','F')
axs[1].set_yticks([0,2,4,6,8,10,12,14,16,18,20]); axs[1].axhline(1,color='#999',lw=.9,ls='--',zorder=0)
axs[1].axhspan(1/1.5,1.5,color=CS,alpha=.08,zorder=0); axs[1].set_xlabel('Command holding')
fig.text(.01,.005,'Only the Cs-gluconate recordings have an intermediate-holding E/I; the paper and the current-lab data were measured at one holding per cell. Shaded: ratio within x1/1.5 to x1.5 of the separate measurement.',fontsize=7.3,color='#444')
fig.tight_layout(rect=(0,.04,1,1)); fig.savefig('Figure4_lab_standards_EI_intermediate.png',dpi=200); plt.close(fig)
pd.DataFrame(BOX,columns=['panel','group','n_cells','mean','median','q1','q3','min','max']).round(4).to_csv('lab_box_stats.csv',index=False)
