"""Lab-standards comparison figures with significance stars (Cs-gluconate vs the other groups).
Frequency (Fig 1 equivalent) + amplitude, rise time, charge, half-decay, half-width.
Stars: Mann-Whitney U, Holm-corrected within each measure/polarity over all 6 pairwise tests (4 groups) or the single test (2 groups).
Only Cs-vs-other comparisons that are significant (Holm p<0.05) get a bracket; * <0.05, ** <0.01, *** <0.001.
Run: python3 fig_lab5.py   (needs /home/claude/work/style/figstyle.py, lab/lab_cells_long.csv, lab/lab_stats.json, lab5/lab_measures_long.csv)"""
import sys; sys.path.insert(0,'/home/claude/work/style')
import numpy as np, pandas as pd, json, itertools, textwrap
from scipy.stats import mannwhitneyu
from figstyle import *
setup(); L='/home/claude/work/lab/'; W='/home/claude/work/lab5/'; O=W+'out/'
D=pd.read_csv(L+'lab_cells_long.csv'); S=json.load(open(L+'lab_stats.json')); M=pd.read_csv(W+'lab_measures_long.csv')
PAST,CUR,CS=PU,RD,GR
SAM="Sam's processed EPSC, control cells (his pipeline)"; AE='Advaith EPSC, Sam experiment (miniML)'; HI='Helena KCl IPSC (miniML)'; AI='Advaith KCl IPSC, Sam experiment (miniML)'
CSE='Cs-gluconate, E at -80 mV (miniML)'; CSI='Cs-gluconate, I at 0 mV (miniML)'; P9='Paper DIV9 control'; P21='Paper DIV21 control'
CURS={'E':[SAM,AE],'I':[HI,AI]}; CSS={'E':CSE,'I':CSI}
LAB4=['Past lab\n(paper, DIV9)','Past lab\n(paper, DIV21)',{'E':'Current lab\n(Sam + Advaith,\nK-Glu, -60 mV)','I':'Current lab\n(Helena + Advaith,\nKCl, -60 mV)'},{'E':'Advaith Cs-Glu\n(-80 mV)','I':'Advaith Cs-Glu\n(0 mV)'}]
def holm(ps):
    ps=np.array(ps); o=np.argsort(ps); out=np.empty(len(ps)); run=0
    for r,i in enumerate(o): run=max(run,(len(ps)-r)*ps[i]); out[i]=min(1,run)
    return out
def star(p): return '***' if p<.001 else '**' if p<.01 else '*' if p<.05 else 'ns'
def mv(var,meas,srcs): return M[(M.variable==var)&(M.measure==meas)&(M.source.isin(srcs))].value.values
def groups(var,meas,paper=True,cur=None):
    """returns list of (label, values, colour) ; index of Cs is last"""
    g=[]
    if paper: g+= [(LAB4[0],mv(var,meas,[P9]),PAST),(LAB4[1],mv(var,meas,[P21]),PAST)]
    lc=LAB4[2][meas]
    if cur==[AE]: lc="Current lab\n(Advaith only,\nK-Glu, -60 mV)"
    g+=[(lc,mv(var,meas,cur or CURS[meas]),CUR),(LAB4[3][meas],mv(var,meas,[CSS[meas]]),CS)]
    return g
def tests(g,labels=None):
    """MW for every pair, Holm over all pairs; returns dict {(i,Cs_idx): (p_raw,p_holm)} and full table"""
    pairs=list(itertools.combinations(range(len(g)),2))
    pr=[mannwhitneyu(g[i][1],g[j][1],alternative='two-sided').pvalue for i,j in pairs]; ph=holm(pr)
    return pairs,pr,ph
def brackets(ax,g,pairs,ph,hi,step,cs=None):
    cs=len(g)-1; sig=[(i,p) for (i,j),p in zip(pairs,ph) if j==cs and p<.05]
    sig.sort(key=lambda t:-t[0])          # nearest group first (lowest bracket)
    y=hi
    for i,p in sig:
        y+=step; ax.plot([i,i,cs,cs],[y-step*.3,y,y,y-step*.3],color='#333',lw=1,clip_on=False)
        ax.text((i+cs)/2,y,star(p),ha='center',va='bottom',fontsize=11,fontweight='bold',color='#222',clip_on=False)
    return hi+(len(sig)+1.1)*step if sig else hi
ROWS=[]
def panel(ax,g,pairs,ph,ylab,title,tag,zoomhi=None,xtick=8):
    lo,hi=full_ylim(g) if zoomhi is None else (0,zoomhi)
    step=(hi-lo)*.085
    n=sum(1 for (i,j),p in zip(pairs,ph) if j==len(g)-1 and p<.05)
    box_panel(ax,g,(lo,hi+(n+1.1)*step),ylab,title,'',stats=None,xtick_fs=xtick,zoom=zoomhi is not None)
    brackets(ax,g,pairs,ph,hi,step)
    return n
def pair_table(var,meas,g,pairs,pr,ph):
    for (i,j),a,b in zip(pairs,pr,ph):
        ROWS.append(dict(variable=var,measure=meas,a=g[i][0].replace('\n',' '),b=g[j][0].replace('\n',' '),n_a=len(g[i][1]),n_b=len(g[j][1]),
                         med_a=np.median(g[i][1]),med_b=np.median(g[j][1]),p=a,p_holm=b,stars=star(b),cs_pair=(j==len(g)-1)))
def two_panel(fname,var,ylabE,ylabI,namE,namI,paper=True,zoom=None,foot_txt='',cur=None,units=''):
    gs={m:groups(var,m,paper,cur and cur[m]) for m in 'EI'}
    st={m:tests(gs[m]) for m in 'EI'}
    for m in 'EI': pair_table(var,m,gs[m],*st[m])
    rows=2 if zoom else 1
    fig,axs=plt.subplots(rows,2,figsize=(12.5,5.4*rows),squeeze=False)
    for k,m in enumerate('EI'):
        pairs,pr,ph=st[m]
        panel(axs[0,k],gs[m],pairs,ph,ylabE if m=='E' else ylabI,('A  ' if m=='E' else 'B  ')+(namE if m=='E' else namI)+', full scale',var)
        if zoom:
            panel(axs[1,k],gs[m],pairs,ph,ylabE if m=='E' else ylabI,('C  ' if m=='E' else 'D  ')+(namE if m=='E' else namI)+', zoom',var,zoomhi=zoom[m])
    foot_txt='\n'.join(textwrap.fill(l,175) for l in foot_txt.split('\n'))
    foot(fig,foot_txt,y=.005); save(fig,O+fname,rect=(0,.095 if rows==1 else .04,1,1))
    return gs,st
FOOTSTAR='Stars: Mann-Whitney U, Holm-corrected over all pairwise tests in the panel; * p<0.05, ** p<0.01, *** p<0.001. A bracket is drawn only when Cs-gluconate differs significantly from that group; no bracket = not significant.'
# ---------- Figure 1: frequency (same as before + stars from lab_stats.json) ----------
def vf(src,m): return D[(D.source==src)&(D.measure==m)].freq.values
def vvf(srcs,m): return D[(D.source.isin(srcs))&(D.measure==m)].freq.values
gF={'E':[(LAB4[0],vf(P9,'E'),PAST),(LAB4[1],vf(P21,'E'),PAST),(LAB4[2]['E'],vvf(CURS['E'],'E'),CUR),(LAB4[3]['E'],vf(CSE,'E'),CS)],
    'I':[(LAB4[0],vf(P9,'I'),PAST),(LAB4[1],vf(P21,'I'),PAST),(LAB4[2]['I'],vvf(CURS['I'],'I'),CUR),(LAB4[3]['I'],vf(CSI,'I'),CS)]}
fig,axs=plt.subplots(1,2,figsize=(12.5,5.6))
for k,m in enumerate('EI'):
    g=gF[m]; pairs=list(itertools.combinations(range(4),2))
    T={(t['a'],t['b']):t['p_holm'] for t in S['tests'][m]}
    names=[t['a'] for t in S['tests'][m][:3]]; keys=[('Paper DIV9','Paper DIV21'),('Paper DIV9','Current lab (pooled)'),('Paper DIV9',[t['b'] for t in S['tests'][m] if 'Cs' in t['b']][0]),('Paper DIV21','Current lab (pooled)'),('Paper DIV21',[t['b'] for t in S['tests'][m] if 'Cs' in t['b']][0]),('Current lab (pooled)',[t['b'] for t in S['tests'][m] if 'Cs' in t['b']][0])]
    ph=[T[kk] for kk in keys]; pr=[next(t['p'] for t in S['tests'][m] if (t['a'],t['b'])==kk) for kk in keys]
    for (i,j),a,b in zip(pairs,pr,ph): ROWS.append(dict(variable='frequency',measure=m,a=g[i][0].replace('\n',' '),b=g[j][0].replace('\n',' '),n_a=len(g[i][1]),n_b=len(g[j][1]),med_a=np.median(g[i][1]),med_b=np.median(g[j][1]),p=a,p_holm=b,stars=star(b),cs_pair=(j==3)))
    panel(axs[k],g,pairs,ph,'EPSC frequency per cell (Hz)' if m=='E' else 'IPSC frequency per cell (Hz)',('A  EPSC' if m=='E' else 'B  IPSC')+' frequency, full scale','frequency')
foot(fig,'Each dot is one cell. Paper values are read off the published plots (control cells only). Different detectors, internals and holdings between groups: see caveats.\n'+FOOTSTAR,y=.005)
save(fig,O+'Lab_Figure1_stars.png',rect=(0,.07,1,1))
# ---------- amplitude, rise time (4 groups) ----------
two_panel('Lab_Fig_amplitude.png','amplitude','EPSC amplitude per cell (pA)','IPSC amplitude per cell (pA)','EPSC amplitude','IPSC amplitude',
  foot_txt='One dot per cell: Sam, Advaith, Helena and Cs = per-cell median of miniML event amplitudes; the paper dots are per-cell means (read off its plots, control cells), which sit above medians for skewed data.\n'+FOOTSTAR)
two_panel('Lab_Fig_risetime.png','risetime','EPSC 10-90% rise time per cell (ms)','IPSC 10-90% rise time per cell (ms)','EPSC rise time','IPSC rise time',
  foot_txt='One dot per cell: paper and Sam = per-cell mean as published; Advaith, Helena and Cs = per-cell mean of miniML events after dropping each cell top 1% (artefact guard), so all groups are means. Paper values read off the published plots (control cells).\n'+FOOTSTAR)
# ---------- charge (current lab vs Cs) ----------
two_panel('Lab_Fig_charge.png','charge','EPSC charge per event (pC)','IPSC charge per event (pC)','EPSC charge','IPSC charge',paper=False,
  foot_txt='No charge values are published in the paper. One dot per cell: per-cell mean charge (Sam as exported; Advaith, Helena, Cs after dropping each cell top 1% of events). Two groups only, so one Mann-Whitney test (no correction needed).\n'+FOOTSTAR)
# ---------- half-decay and half-width (miniML groups only) ----------
curHD={'E':[SAM,AE],'I':[HI,AI]}; curHW={'E':[AE],'I':[HI,AI]}
two_panel('Lab_Fig_halfdecay.png','halfdecay','EPSC half-decay time per cell (ms)','IPSC half-decay time per cell (ms)','EPSC half-decay','IPSC half-decay',paper=False,cur=curHD,
  foot_txt="Only miniML cells can be compared (the paper reports a decay tau, not half-decay). Current-lab EPSC = Sam's 63 cells (his per-cell mean half-decay) + Advaith's 3 cells; Advaith, Helena and Cs are per-cell means after dropping the top 1% of events, so all groups are means. One dot per cell.\n"+FOOTSTAR)
two_panel('Lab_Fig_halfwidth.png','halfwidth','EPSC half-width per cell (ms)','IPSC half-width per cell (ms)','EPSC half-width ("half time")','IPSC half-width ("half time")',paper=False,cur=curHW,
  foot_txt="Only miniML cells can be compared (no half-width in the paper or in Sam's exported sheet). Current-lab EPSC = Advaith's 3 cells only. One dot per cell (median of miniML events).\n"+FOOTSTAR)
pd.DataFrame(ROWS).round(5).to_csv(W+'lab_measure_stats.csv',index=False)
print(pd.DataFrame(ROWS).query('cs_pair').round(4)[['variable','measure','a','b','n_a','n_b','med_a','med_b','p','p_holm','stars']].to_string())
