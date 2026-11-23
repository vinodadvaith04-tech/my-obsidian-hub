import sys; sys.path.insert(0,'/home/claude/work/style')
import numpy as np, pandas as pd, json
from figstyle import *
from matplotlib.gridspec import GridSpec
setup(); C='/home/claude/work/cs/'; O='/home/claude/work/style/out/'
a=pd.read_csv(C+'cs_analysis_set.csv')
H=[-80,-70,-60,-50,-40,-30,0,10]
cm=lambda c: a.groupby(['cellkey','holding_command_mV'])[c].mean().unstack()
xl='Command holding (mV)'
def G(df,col,hs=H): return [(str(h),df[h].dropna().values if h in df.columns else np.array([]),col) for h in hs]
ZN='↑k / ↓k = k points beyond the zoom axis'
zextra=[Line2D([0],[0],color='none',label=ZN)]
# ================= Figure 1: noise, drift, holding current (row 1 full scale, row 2 zoom)
fig,axs=plt.subplots(2,3,figsize=(14,9.4))
spec=[('noise_lp3k_pA',BL,'Noise per cell (pA)','A  Noise (MAD-SD, 3 kHz low-pass)',(0,6.5),'D  Noise, zoom to 0-6.5 pA',False),
      ('drift_total_pA',RD,'Drift per cell (pA over the 2 min)','B  Drift over 2 min, full scale',(-150,100),'E  Drift, zoom to -150 to +100 pA',True),
      ('Ihold_median_pA',GR,'Holding current per cell (pA)','C  Holding current, full scale',(-450,350),'F  Holding current, zoom to -450 to +350 pA',True)]
for j,(c,col,yl,t1,zl,t2,z0) in enumerate(spec):
    g=G(cm(c),col); box_panel(axs[0,j],g,full_ylim(g,zero=not z0),yl,t1,xlabel=xl,zero=z0)
    box_panel(axs[1,j],g,zl,yl,t2,xlabel=xl,zero=z0,zoom=True)
legend(axs[0,0],'upper center',(0.55,1.0),extra=zextra)
foot(fig,'Each dot is one cell (mean of its recordings at that command holding). -20 mV and +15 mV (1 cell each) are in the tables but not drawn. Drift is signed: negative = current falling over the 2 min.\nUpper row: all cells, full scale. Lower row: the same data on a narrower scale.',.004)
save(fig,O+'Cs_Figure1.png',rect=(0,.045,1,1))
# ================= Figure 2: event frequency
fi,fo=cm('f_in'),cm('f_out'); yl='miniML frequency per cell (Hz)'
gi=G(fi,BL); go=G(fo,RD)
fig,axs=plt.subplots(2,2,figsize=(12.6,8.8)); ST=[]
box_panel(axs[0,0],gi,full_ylim(gi),yl,'A  Inward events (EPSC), full scale',xlabel=xl)
box_panel(axs[0,1],go,full_ylim(go),yl,'B  Outward events (IPSC), full scale',xlabel=xl)
box_panel(axs[1,0],gi,(0,.4),yl,'C  Inward events, zoom to 0-0.4 Hz',xlabel=xl,zoom=True)
box_panel(axs[1,1],go,(0,.4),yl,'D  Outward events, zoom to 0-0.4 Hz',xlabel=xl,zoom=True)
for ax in (axs[0,0],axs[1,1]): ax.axvspan(-.5,1.5,color=BL,alpha=.06,zorder=0)
for ax in (axs[0,1],axs[1,0]): ax.axvspan(5.5,7.5,color=RD,alpha=.06,zorder=0)
legend(axs[0,1],'upper left',(0.0,1.0),extra=zextra)
foot(fig,'Shaded: candidate isolating holdings. Blue: -80, -70 for EPSCs (A, and D = residual outward events). Red: 0, +10 for IPSCs (B, and C = residual inward events).\nEach dot is one cell (mean of its recordings at that holding). -20 and +15 mV (1 cell each) are not drawn.',.004)
save(fig,O+'Cs_Figure2.png',rect=(0,.05,1,1))
# ================= Figure 3: E/I preservation (linear axes; row 1 full scale, row 2 zoom)
res=json.load(open(C+'cs_tables_part2.json'))['primary']['ei']
pc=pd.read_csv(C+'cs_ei_per_cell_ref-80_0.csv'); mid=[-60,-50,-40,-30]
def retention(): 
    r=[]
    for h in mid:
        t=pc[pc.h==h]; r.append((h,len(t),np.median(t.rE)*100,*np.percentile(t.rE*100,[25,75]),np.median(t.rI)*100,*np.percentile(t.rI*100,[25,75])))
    return pd.DataFrame(r,columns=['h','n','rE','rE25','rE75','rI','rI25','rI75'])
ret=retention()
def ret_panel(ax,ylim,title,zoom):
    x=np.arange(len(mid)); lo,hi=ylim
    ax.errorbar(x-.04,ret.rE,yerr=[ret.rE-ret.rE25,ret.rE75-ret.rE],color=BL,marker='o',lw=1.8,capsize=3,label='EPSC frequency (inward) / that at -80')
    ax.errorbar(x+.04,ret.rI,yerr=[ret.rI-ret.rI25,ret.rI75-ret.rI],color=RD,marker='s',lw=1.8,capsize=3,label='IPSC frequency (outward) / that at 0')
    for k,r in ret.iterrows(): ax.text(k,-0.035,f'n={int(r.n)}',ha='center',va='top',fontsize=8,transform=ax.get_xaxis_transform(),color='#333')
    ax.set_xticks(x); ax.set_xticklabels(mid); ax.tick_params(axis='x',pad=14); ax.set_xlabel(xl,labelpad=2); ax.set_xlim(-.4,len(mid)-.6)
    ax.set_ylim(lo,hi); ax.axhline(100,color='#999',lw=.9,ls='--',zorder=0); ax.set_ylabel('Share of the isolated-holding frequency still detected (%)'); ax.set_title(title,loc='left')
ref=[('-80 | 0',res['-80/0']['cross'],PU),('-70 | 0',res['-70/0']['cross'],PU),('-80 | +10',res['-80/+10' if '-80/+10' in res else '-80/10']['cross'],PU),('-70 | +10',res['-70/+10' if '-70/+10' in res else '-70/10']['cross'],PU)]
fig=plt.figure(figsize=(14.6,9.6)); gs=GridSpec(2,3,figure=fig,width_ratios=[1,1,.9]); ax=[[fig.add_subplot(gs[i,j]) for j in range(3)] for i in range(2)]
gf=[(str(h),pc[pc.h==h].fold.values,PU) for h in mid]
box_panel(ax[0][0],gf,full_ylim(gf),'Within-cell E/I ratio at h / reference E/I ratio (fold)','A  E/I frequency ratio relative to the separately measured one, full scale',xlabel=xl,hline=1,band=(1/1.5,1.5,PU))
box_panel(ax[1][0],gf,(0,20),'Within-cell E/I ratio at h / reference E/I ratio (fold)','D  E/I ratio, zoom to 0-20 fold',xlabel=xl,hline=1,band=(1/1.5,1.5,PU),zoom=True)
ret_panel(ax[0][1],(0,max(ret.rE75)*1.06),'B  How much of each event type is still detected, full scale',False)
ret_panel(ax[1][1],(0,60),'E  Detection, zoom to 0-60%',True)
ax[0][1].legend(frameon=False,loc='upper right')
box_panel(ax[0][2],ref,(-65,0),'Command holding where the cell\'s E/I ratio\nequals the reference (mV, interpolated)','C  Holding that reproduces the separate E/I ratio',xlabel='Reference holdings (E | I)',zero=False)
ax[1][2].axis('off'); legend(ax[1][2],'center',(0.5,0.5),extra=zextra)
for t in ax[0][0].texts: pass
ax[1][1].text(0.02,.97,'\u2191 EPSC at -60 mV: median 138% (off the axis)',transform=ax[1][1].transAxes,fontsize=7.5,va='top',color=BL)
foot(fig,'Each dot is one cell. Ratio = inward/outward event frequency. Reference = inward frequency at -80 mV (E) and outward frequency at 0 mV (I) in the same cell. Grey band in A, D: within x1.5.\nB, E: median and IQR of cells. Lower row (D, E): the same data as the upper row on a narrower scale; there is no zoom for C because all values lie between -60 and 0 mV.',.004)
save(fig,O+'Cs_Figure3.png',rect=(0,.055,1,1))
print(ret.round(1)); print(list(res.keys()))
