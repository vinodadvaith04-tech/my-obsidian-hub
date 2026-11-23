import sys; sys.path.insert(0,'/home/claude/work/style')
import numpy as np, pandas as pd, json
from figstyle import *
setup(); L='/home/claude/work/lab/'; O='/home/claude/work/style/out/'
D=pd.read_csv(L+'lab_cells_long.csv'); S=json.load(open(L+'lab_stats.json'))
PAST=PU; CUR=RD; CS=GR; BOX=[]
def V(src,meas): return D[(D.source==src)&(D.measure==meas)].freq.values
def VV(srcs,meas): return D[(D.source.isin(srcs))&(D.measure==meas)].freq.values
SAM="Sam's processed EPSC, control cells (his pipeline)"; AE='Advaith EPSC, Sam experiment (miniML)'; HI='Helena KCl IPSC (miniML)'; AI='Advaith KCl IPSC, Sam experiment (miniML)'
yE='EPSC frequency per cell (Hz)'; yI='IPSC frequency per cell (Hz)'
# ---- Figure 1
gE=[('Past lab\n(paper, DIV9)',V('Paper DIV9 control','E'),PAST),('Past lab\n(paper, DIV21)',V('Paper DIV21 control','E'),PAST),('Current lab\n(Sam + Advaith,\nK-Glu, -60 mV)',VV([SAM,AE],'E'),CUR),('Advaith Cs-Glu\n(-80 mV)',V('Cs-gluconate, E at -80 mV (miniML)','E'),CS)]
gI=[('Past lab\n(paper, DIV9)',V('Paper DIV9 control','I'),PAST),('Past lab\n(paper, DIV21)',V('Paper DIV21 control','I'),PAST),('Current lab\n(Helena + Advaith,\nKCl, -60 mV)',VV([HI,AI],'I'),CUR),('Advaith Cs-Glu\n(0 mV)',V('Cs-gluconate, I at 0 mV (miniML)','I'),CS)]
fig,axs=plt.subplots(1,2,figsize=(11.6,5.0))
box_panel(axs[0],gE,full_ylim(gE),yE,'A  EPSC frequency, full scale','',stats=BOX,tag='E1',xtick_fs=8)
box_panel(axs[1],gI,full_ylim(gI),yI,'B  IPSC frequency, full scale','',stats=BOX,tag='I1',xtick_fs=8)
legend(axs[1],'upper left',(0.0,1.0))
foot(fig,'Each dot is one cell. Paper values are read off the published plots (control cells only). Different detectors, internals and holdings between groups: see caveats.'); save(fig,O+'Lab_Figure1.png')
# ---- Figure 2
gE=[('Paper\nDIV9',V('Paper DIV9 control','E'),PAST),('Paper\nDIV21',V('Paper DIV21 control','E'),PAST),("Sam's own\nprocessed data",V(SAM,'E'),CUR),("Advaith,\nSam's experiment",V(AE,'E'),CUR),('Advaith\nCs-Glu, -80',V('Cs-gluconate, E at -80 mV (miniML)','E'),CS)]
gI=[('Paper\nDIV9',V('Paper DIV9 control','I'),PAST),('Paper\nDIV21',V('Paper DIV21 control','I'),PAST),("Helena\nKCl",V(HI,'I'),CUR),("Advaith KCl,\nSam's experiment",V(AI,'I'),CUR),('Advaith\nCs-Glu, 0',V('Cs-gluconate, I at 0 mV (miniML)','I'),CS)]
fig,axs=plt.subplots(1,2,figsize=(12.6,5.0))
box_panel(axs[0],gE,full_ylim(gE),yE,'A  EPSC frequency by data source, full scale',stats=BOX,tag='E2',xtick_fs=8)
box_panel(axs[1],gI,full_ylim(gI),yI,'B  IPSC frequency by data source, full scale',stats=BOX,tag='I2',xtick_fs=8)
legend(axs[1],'upper left',(0.0,1.0))
foot(fig,'Same data as Figure 1 with the current-lab sources shown separately. Each dot is one cell.'); save(fig,O+'Lab_Figure2.png')
# ---- Figure 3
R=S['ratio']; cw=pd.read_csv(L+'cs_E_I_same_cell.csv')
fig,ax=plt.subplots(figsize=(8.8,5.2))
items=[('Past lab\n(paper, DIV9)','Paper DIV9',PAST),('Past lab\n(paper, DIV21)','Paper DIV21',PAST),('Current lab\n(Sam + Advaith E,\nHelena + Advaith I)','Current lab (pooled)',CUR),('Advaith Cs-Glu\n(E at -80, I at 0,\nseparate cells)','Cs-gluconate (-80 / 0 mV, separate cells)',CS)]
for k,(lab,key,col) in enumerate(items):
    r=R[key]; ax.errorbar([k],[r['ratio_mean']],yerr=[[r['ratio_mean']-r['ci_mean'][0]],[r['ci_mean'][1]-r['ratio_mean']]],fmt='D',ms=8,color=col,mfc='white',mew=1.8,capsize=5,lw=1.6,zorder=4)
    ax.text(k+.12,r['ratio_mean'],f"{r['ratio_mean']:.2f}",va='center',fontsize=8.5,color='#333')
rng=np.random.default_rng(5); v=cw.ratio.values; box(ax,4,v,CS,rng); med=np.median(v); ax.text(4+.36,med,f'{med:.2f}',va='center',fontsize=8.5,color='#333')
ax.set_ylim(0,max(1.6,v.max()*1.06)); ax.set_xticks(range(5)); ax.set_xticklabels([i[0] for i in items]+['Advaith Cs-Glu\n(same cell,\nper cell)'],fontsize=8); ax.tick_params(axis='x',pad=4)
ax.set_xlim(-.6,4.7); ax.set_ylabel('E/I frequency ratio (EPSC Hz / IPSC Hz)'); ax.axhline(1,color='#999',lw=.8,ls='--',zorder=0)
ax.set_title('E/I ratio from separately measured EPSCs and IPSCs',loc='left'); ax.axvline(3.5,color='#ccc',lw=.8)
foot(fig,'Diamonds: ratio of group mean frequencies with 95% bootstrap interval (E and I come from different cells).\nRight: Cs, ratio within each cell (n=14): box = IQR, bar = median, diamond = mean, dots = cells.',.004); save(fig,O+'Lab_Figure3.png',rect=(0,.07,1,1))
# ---- Figure 4
ei=pd.read_csv('/home/claude/work/cs/cs_ei_per_cell_ref-80_0.csv'); mid=[-50,-40,-30]
fig,axs=plt.subplots(1,3,figsize=(15.4,5.0),gridspec_kw=dict(width_ratios=[1,1.1,1.1]))
t=ei[ei.h==-40]; ax=axs[0]
for _,r in t.iterrows(): ax.plot([0,1],[r.ratio_ref,r.ratio_h],color='#999',lw=.9,alpha=.8,zorder=1)
ax.scatter(np.zeros(len(t)),t.ratio_ref,s=22,color=CS,ec='white',lw=.4,zorder=3); ax.scatter(np.ones(len(t)),t.ratio_h,s=22,color=CS,ec='white',lw=.4,zorder=3)
for x,c in ((0,'ratio_ref'),(1,'ratio_h')): ax.plot([x-.18,x+.18],[t[c].median()]*2,color='k',lw=2.4,zorder=4)
ax.set_ylim(0,max(3,t[['ratio_ref','ratio_h']].values.max()*1.06)); ax.axhline(1,color='#999',lw=.8,ls='--',zorder=0)
ax.set_xticks([0,1]); ax.set_xticklabels(['Measured separately\n(E at -80, I at 0)','Measured together\nat -40 mV']); ax.set_xlim(-.5,1.5)
ax.set_ylabel('E/I frequency ratio per cell'); ax.set_title(f'A  Same cells, two ways (n={len(t)})',loc='left')
ax.text(.5,.04,'line = one cell; black bar = median',transform=ax.transAxes,ha='center',fontsize=7.5,color='#444')
g=[(f'{h} mV',ei[ei.h==h].fold.values,CS) for h in mid]; yl='Ratio at holding / separately measured ratio (fold)'
box_panel(axs[1],g,full_ylim(g),yl,'B  How far each holding moves the ratio, full scale','Command holding',hline=1,band=(1/1.5,1.5,CS),stats=BOX,tag='F')
box_panel(axs[2],g,(0,6),yl,'C  Zoom to 0-6 fold','Command holding',hline=1,band=(1/1.5,1.5,CS),zoom=True)
legend(axs[1],'upper right',extra=[Line2D([0],[0],color='none',label='↑k = k points above the zoom axis')])
foot(fig,'Only the Cs-gluconate recordings have an intermediate-holding E/I; the paper and the current-lab data were measured at one holding per cell. Shaded: ratio within x1/1.5 to x1.5 of the separate measurement. C: same data as B on a narrower scale.'); save(fig,O+'Lab_Figure4.png')
pd.DataFrame(BOX,columns=['panel','group','n_cells','mean','median','q1','q3','min','max']).round(4).to_csv(O+'lab_box_stats_check.csv',index=False)
# ---- Figure 5
from scipy import stats as st_
G=pd.read_csv('/home/claude/work/lab2/audit_groups_windows.csv'); G['cell']=G.file.str.extract(r'^(\d{8}_[^_]+_[^_]+_[^_]+_slice\d+_C\d+)')[0]; Gr=G[G.why=='random']
cs=Gr[Gr.group=='I: Cs, 0 and +10 mV']; kc=Gr[Gr.group=='I: KCl + DNQX/AP5, -60 mV']
cr=cs[cs.truth>=3]; kr=kc[kc.truth>=3]; cr=cr.assign(rec=cr.mini_tp/cr.truth); kr=kr.assign(rec=kr.mini_tp/kr.truth)
gA=[('Cs, 0/+10 mV\nminiML',cs.mini_markers_hz,CS,.10,'--','s'),('Cs, 0/+10 mV\naudited',cs.truth_hz,CS,.25),('KCl + DNQX/AP5\nminiML',kc.mini_markers_hz,CUR,.10,'--','s'),('KCl + DNQX/AP5\naudited',kc.truth_hz,CUR,.25)]
gB=[('Cs, 0/+10 mV',cr.rec,CS),('KCl + DNQX/AP5',kr.rec,CUR)]
fig,axs=plt.subplots(1,2,figsize=(12,5.0),gridspec_kw=dict(width_ratios=[1.25,1]))
box_panel(axs[0],gA,full_ylim(gA),'IPSC frequency per 3 s window (Hz)','A  IPSC frequency in the audited windows',seed=3,xtick_fs=8)
box_panel(axs[1],gB,(0,1.05),'Share of audited IPSCs found by miniML','B  miniML detection rate per window',hline=1,seed=3,xtick_fs=8)
axs[0].legend(handles=[Line2D([0],[0],marker='s',color='w',mfc='#ccc',mec='#555',ms=7,label='miniML count (dashed box)'),Line2D([0],[0],marker='o',color='w',mfc='#555',ms=6,label='Audited count (solid box)'),Line2D([0],[0],marker='D',color='w',mfc='white',mec='#555',mew=1.5,ms=6,label='Mean')],frameon=False,loc='upper left')
foot(fig,'Each dot is one 3 s window, random windows only. Audited = miniML markers you accepted + second-detector markers + events you clicked. Panel B: windows with at least 3 audited events.'); save(fig,O+'Lab_Figure5.png')
print('ok', len(cs),len(kc),len(cr),len(kr))
