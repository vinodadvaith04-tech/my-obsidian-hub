"""Audit-based figures for the Lab standards report (uses the latest audit_labels.xlsx via audit_groups_windows.csv):
 1) precision ('accuracy') and recall of miniML per group, random windows only (unbiased), 95% cluster-bootstrap CI over files
 2) IPSC frequency from audited windows only (per cell, random 3-s windows)  -> stars Cs vs others
 3) projected IPSC frequency if miniML were perfect = per-cell miniML frequency x (audited true events / miniML markers) of the group
Old K-Glu no-blocker recordings and 'KCl blockers not confirmed' are excluded. Run: python3 fig_audit5.py"""
import sys; sys.path.insert(0,'/home/claude/work/style'); sys.path.insert(0,'/home/claude/work/lab5')
import numpy as np, pandas as pd, re, itertools
from scipy.stats import mannwhitneyu
from figstyle import *
setup(); W='/home/claude/work/lab5/'; O=W+'out/'; L='/home/claude/work/lab/'
rng=np.random.default_rng(7); B=5000
D=pd.read_csv(L+'lab_cells_long.csv'); G=pd.read_csv(W+'audit_groups_windows.csv')
PAST,CUR,CS=PU,RD,GR
def ck(f):
    b=f.rsplit('_',1)[0]; m=re.match(r'(\d{8})_.*_slice(\d+)_C(\d+)$',b)
    return f"{m.group(1)}_s{m.group(2)}_c{m.group(3)}" if m else b
G['cell']=G.file.map(lambda f: ck(f)); G.loc[G.group.str.contains('KCl|K-Glu'),'cell']=G.file.str.rsplit('_',n=1).str[0]
R=G[G.why=='random'].copy()
GI_CS='I: Cs, 0 and +10 mV'; GI_K='I: KCl + DNQX/AP5, -60 mV'; GE_CS='E: Cs, -80 and -70 mV'; GE_K='E: K-Glu + gabazine, -60 mV'
def boot_ratio(x,num,den,B=B):
    cl=[g for _,g in x.groupby('file')]; n=len(cl); out=[]
    for _ in range(B):
        s=pd.concat([cl[i] for i in rng.integers(0,n,n)]); out.append(s[num].sum()/s[den].sum())
    return np.percentile(out,[2.5,97.5])
rows=[]
for g in (GI_CS,GI_K,GE_CS,GE_K):
    x=R[R.group==g]
    rec=x.mini_tp.sum()/x.truth.sum(); pre=x.mini_tp.sum()/x.mini_markers.sum(); sc=x.truth.sum()/x.mini_markers.sum()
    # recall/precision CIs
    cl=[h for _,h in x.groupby('file')]; n=len(cl); rr=[];pp=[]
    for _ in range(B):
        s=pd.concat([cl[i] for i in rng.integers(0,n,n)]); rr.append(s.mini_tp.sum()/s.truth.sum()); pp.append(s.mini_tp.sum()/s.mini_markers.sum())
    rows.append(dict(group=g,windows=len(x),files=x.file.nunique(),true_events=int(x.truth.sum()),miniML_markers=int(x.mini_markers.sum()),found=int(x.mini_tp.sum()),
        recall=rec,recall_lo=np.percentile(rr,2.5),recall_hi=np.percentile(rr,97.5),precision=pre,prec_lo=np.percentile(pp,2.5),prec_hi=np.percentile(pp,97.5),
        scale_true_over_markers=sc,scale_lo=boot_ratio(x,'truth','mini_markers')[0],scale_hi=boot_ratio(x,'truth','mini_markers')[1],
        audited_hz=x.truth.sum()/x.dur.sum(),miniML_hz_same_windows=x.mini_markers.sum()/x.dur.sum()))
AU=pd.DataFrame(rows); AU.to_csv(W+'audit_precision_recall.csv',index=False); print(AU.round(3).to_string())
# ---------------- Figure A: precision & recall ----------------
fig,axs=plt.subplots(1,2,figsize=(12.5,5.2))
labs=['I: Cs-Glu\n(0, +10 mV)','I: KCl + DNQX/AP5\n(-60 mV)','E: Cs-Glu\n(-80, -70 mV)','E: K-Glu +\ngabazine (-60 mV)']
cols=[CS,CUR,CS,CUR]
for ax,(k,lo,hi,t) in zip(axs,[('recall','recall_lo','recall_hi','A  Recall: share of real events that miniML finds'),('precision','prec_lo','prec_hi','B  Precision ("accuracy"): share of miniML events that are real')]):
    for i,(r,c) in enumerate(zip(AU.itertuples(),cols)):
        v=getattr(r,k); l=getattr(r,lo); h=getattr(r,hi)
        ax.bar(i,v,width=.55,color=c,alpha=.25,edgecolor=c,lw=1.4)
        ax.errorbar(i,v,yerr=[[v-l],[h-v]],color=c,lw=2,capsize=6)
        ax.text(i,min(h+.03,1.06),f'{v*100:.0f}%',ha='center',va='bottom',fontsize=11,fontweight='bold',color='#222')
        ax.text(i,-.035,f'{r.windows} windows, {r.files} files\n{r.true_events} real events',ha='center',va='top',fontsize=7.5,color='#333',transform=ax.get_xaxis_transform())
    ax.set_xticks(range(4)); ax.set_xticklabels(labs,fontsize=8); ax.tick_params(axis='x',pad=34)
    ax.set_ylim(0,1.15); ax.set_yticks(np.arange(0,1.01,.2)); ax.set_yticklabels([f'{int(v*100)}%' for v in np.arange(0,1.01,.2)]); ax.set_xlim(-.6,3.6)
    ax.set_title(t,loc='left',fontsize=9.5); ax.axhline(1,color='#999',lw=.7,ls='--',zorder=0)
axs[0].set_ylabel('Recall'); axs[1].set_ylabel('Precision')
foot(fig,'Random 3-s audit windows only (unbiased); old K-Glu no-blocker recordings and recordings with unconfirmed blockers excluded. Bar = pooled over windows; whiskers = 95% bootstrap CI resampling files.\nRecall = real events found by miniML / all real events; precision = miniML events that are real / all miniML events. Real events = audited truth (miniML hits confirmed + events added by hand).',y=.005)
save(fig,O+'Lab_Fig_audit_precision_recall.png',rect=(0,.07,1,1))
# ---------------- per-cell tables ----------------
def holm(ps):
    ps=np.array(ps); o=np.argsort(ps); out=np.empty(len(ps)); run=0
    for r,i in enumerate(o): run=max(run,(len(ps)-r)*ps[i]); out[i]=min(1,run)
    return out
def star(p): return '***' if p<.001 else '**' if p<.01 else '*' if p<.05 else 'ns'
def vf(src): return D[(D.source==src)&(D.measure=='I')].set_index('cell').freq
P9=vf('Paper DIV9 control'); P21=vf('Paper DIV21 control'); HI=vf('Helena KCl IPSC (miniML)'); AI=vf('Advaith KCl IPSC, Sam experiment (miniML)'); CSF=vf('Cs-gluconate, I at 0 mV (miniML)')
KC=pd.concat([HI,AI]); 
# audited per-cell (random windows; Cs = 0 mV windows only so the cells match the frequency figure)
rc=R[(R.group==GI_CS)&(R.holding==0)].groupby('cell').agg(aud=('truth_hz','mean'),mini_same=('mini_markers_hz','mean'),w=('window','count'))
rk=R[R.group==GI_K].groupby('cell').agg(aud=('truth_hz','mean'),mini_same=('mini_markers_hz','mean'),w=('window','count'))
assert set(rc.index)==set(CSF.index) and set(rk.index)==set(KC.index)
def pairs_holm(g): 
    pr=[];pairs=list(itertools.combinations(range(4),2))
    for i,j in pairs: pr.append(mannwhitneyu(g[i][1],g[j][1]).pvalue)
    return pairs,pr,holm(pr)
def draw_brackets(ax,hi,step,pairs,ph):
    sig=sorted([(i,p) for (i,j),p in zip(pairs,ph) if j==3 and p<.05],key=lambda t:-t[0]); y=hi
    for i,p in sig:
        y+=step; ax.plot([i,i,3,3],[y-step*.3,y,y,y-step*.3],color='#333',lw=1,clip_on=False); ax.text((i+3)/2,y,star(p),ha='center',va='bottom',fontsize=11,fontweight='bold',color='#222')
    return len(sig)
XL=['Past lab\n(paper, DIV9)','Past lab\n(paper, DIV21)']
ROWS=[]
def fig_groups(g,tag,ylab,title,ax):
    pairs,pr,ph=pairs_holm(g); lo,hi=full_ylim(g); step=(hi-lo)*.085; n=sum(1 for (i,j),p in zip(pairs,ph) if j==3 and p<.05)
    box_panel(ax,g,(lo,hi+(n+1.1)*step),ylab,title,'',xtick_fs=8); draw_brackets(ax,hi,step,pairs,ph)
    for (i,j),a,b in zip(pairs,pr,ph): ROWS.append(dict(figure=tag,a=g[i][0].replace('\n',' '),b=g[j][0].replace('\n',' '),n_a=len(g[i][1]),n_b=len(g[j][1]),med_a=np.median(g[i][1]),med_b=np.median(g[j][1]),p=a,p_holm=b,stars=star(b),cs_pair=(j==3)))
# ---------------- Figure B: audited windows only ----------------
gB=[(XL[0],P9.values,PAST),(XL[1],P21.values,PAST),('Current lab KCl\n(audited windows,\n-60 mV)',rk.aud.values,CUR),('Advaith Cs-Glu\n(audited windows,\n0 mV)',rc.aud.values,CS)]
fig,axs=plt.subplots(1,2,figsize=(12.5,5.6),gridspec_kw=dict(width_ratios=[1.35,1]))
fig_groups(gB,'audited','IPSC frequency per cell (Hz)','A  IPSC frequency, audited windows only','' and axs[0]) if False else None
fig_groups(gB,'audited','IPSC frequency per cell (Hz)','A  IPSC frequency from audited windows only',axs[0])
ax=axs[1]; mx=max(rc.aud.max(),rc.mini_same.max(),rk.aud.max(),rk.mini_same.max())*1.06
ax.plot([0,mx],[0,mx],color='#999',ls='--',lw=1,zorder=0)
ax.scatter(rk.mini_same,rk.aud,s=26,color=CUR,alpha=.85,label=f'KCl cells (n={len(rk)})',edgecolor='w',lw=.5)
ax.scatter(rc.mini_same,rc.aud,s=26,color=CS,alpha=.85,label=f'Cs-Glu cells (n={len(rc)})',edgecolor='w',lw=.5)
ax.set_xlim(0,mx); ax.set_ylim(0,mx); ax.set_aspect('equal'); ax.set_xlabel('miniML IPSC frequency in the same windows (Hz)'); ax.set_ylabel('Audited (true) IPSC frequency (Hz)')
ax.set_title('B  Per cell: audited vs miniML in the same windows',loc='left'); ax.legend(frameon=False,loc='upper left'); ax.text(mx*.02,mx*.62,'dashed line:\nminiML = audit',fontsize=7.5,color='#666')
foot(fig,'Panel A: one dot per cell; paper = published control cells. Current-lab and Cs values = mean over that cell\'s random 3-s audit windows\n(KCl: 24 cells, 42 windows; Cs at 0 mV: 15 cells, 46 windows; 1-4 windows per cell, so each cell is a noisy estimate).\n'+'Stars: Mann-Whitney U, Holm over all 6 pairs; * p<0.05, ** p<0.01, *** p<0.001; bracket only when Cs differs significantly from that group.',y=.005)
save(fig,O+'Lab_Fig_audited_ipsc.png',rect=(0,.1,1,1))
# ---------------- Figure C: projected perfect miniML ----------------
fk=AU.set_index('group').scale_true_over_markers; fkl=AU.set_index('group').scale_lo; fkh=AU.set_index('group').scale_hi
projC=CSF*fk[GI_CS]; projK=KC*fk[GI_K]
gC=[(XL[0],P9.values,PAST),(XL[1],P21.values,PAST),(f'Current lab KCl\n(miniML x {fk[GI_K]:.2f},\n-60 mV)',projK.values,CUR),(f'Advaith Cs-Glu\n(miniML x {fk[GI_CS]:.2f},\n0 mV)',projC.values,CS)]
fig,axs=plt.subplots(1,2,figsize=(12.5,5.6),gridspec_kw=dict(width_ratios=[1.35,1]))
fig_groups(gC,'projected','IPSC frequency per cell (Hz)','A  IPSC frequency if miniML were perfect',axs[0])
ax=axs[0]
for k,(raw) in {2:KC.values,3:CSF.values}.items(): ax.plot([k-.36,k+.36],[np.median(raw)]*2,color='#222',lw=1.4,ls=':',zorder=5)
ax.plot([],[],color='#222',lw=1.4,ls=':',label='median as measured by miniML (before scaling)'); ax.legend(frameon=False,loc='upper left',fontsize=7.5)
# panel B: group medians, three scenarios, with CI
def boot_med(vals,fac_lo_hi_samples=None):
    out=[]
    for _ in range(B):
        s=rng.choice(vals,len(vals)); out.append(np.median(s)*(1 if fac_lo_hi_samples is None else rng.choice(fac_lo_hi_samples)))
    return np.percentile(out,[2.5,97.5])
def fac_samples(g):
    x=R[R.group==g]; cl=[h for _,h in x.groupby('file')]; n=len(cl); out=[]
    for _ in range(2000):
        s=pd.concat([cl[i] for i in rng.integers(0,n,n)]); out.append(s.truth.sum()/s.mini_markers.sum())
    return np.array(out)
fsC=fac_samples(GI_CS); fsK=fac_samples(GI_K)
sc=[('miniML as\nmeasured',np.median(KC),boot_med(KC.values),np.median(CSF),boot_med(CSF.values)),
    ('audited\nwindows only',rk.aud.median(),boot_med(rk.aud.values),rc.aud.median(),boot_med(rc.aud.values)),
    ('projected,\nperfect miniML',np.median(projK),boot_med(KC.values,fsK*1.0)if False else boot_med(KC.values,fsK),np.median(projC),boot_med(CSF.values,fsC))]
ax=axs[1]; SC=[]
for i,(n,mk,ck_,mc,cc) in enumerate(sc):
    ax.errorbar(i-.12,mk,yerr=[[mk-ck_[0]],[ck_[1]-mk]],fmt='o',color=CUR,ms=8,capsize=4,lw=1.8,label='Current lab (KCl, -60 mV)' if i==0 else None)
    ax.errorbar(i+.12,mc,yerr=[[mc-cc[0]],[cc[1]-mc]],fmt='s',color=CS,ms=8,capsize=4,lw=1.8,label='Advaith Cs-Glu (0 mV)' if i==0 else None)
    ax.text(i-.12,ck_[1]+.35,f'{mk:.1f}',ha='center',fontsize=8.5,color=CUR); ax.text(i+.12,cc[1]+.35,f'{mc:.1f}',ha='center',fontsize=8.5,color=CS)
    SC.append(dict(scenario=n.replace('\n',' '),kcl_median=mk,kcl_lo=ck_[0],kcl_hi=ck_[1],cs_median=mc,cs_lo=cc[0],cs_hi=cc[1],ratio=mk/mc))
ax.set_xticks(range(3)); ax.set_xticklabels([s[0]+f'\nKCl / Cs = {s[1]/s[3]:.1f}x' for s in sc],fontsize=8.5); ax.tick_params(axis='x',pad=4); ax.set_xlim(-.6,2.6); ax.set_ylim(0,max(s['kcl_hi'] for s in SC)*1.15)
ax.set_ylabel('Median IPSC frequency per cell (Hz)'); ax.set_title('B  Median per group under three views',loc='left'); ax.legend(frameon=False,loc='lower left',fontsize=7.5)
pd.DataFrame(SC).round(3).to_csv(W+'ipsc_scenario_medians.csv',index=False)
foot(fig,'Projection: each cell\'s miniML frequency x (real events / miniML events) from that group\'s random audit windows\n(Cs 72 windows, KCl 42 windows); paper values are not scaled. Whiskers in B: 95% bootstrap CI (cells and audit files resampled).\nStars: Mann-Whitney U, Holm over all 6 pairs; bracket only when Cs differs significantly from that group.',y=.005)
save(fig,O+'Lab_Fig_projected_ipsc.png',rect=(0,.1,1,1))
pd.DataFrame(ROWS).round(5).to_csv(W+'lab_audit_ipsc_stats.csv',index=False)
print(pd.DataFrame(ROWS).query('cs_pair')[['figure','a','n_a','n_b','med_a','med_b','p','p_holm','stars']].round(4).to_string()); print(pd.DataFrame(SC).round(2).to_string())
rc.to_csv(W+'audited_per_cell_cs.csv'); rk.to_csv(W+'audited_per_cell_kcl.csv')
