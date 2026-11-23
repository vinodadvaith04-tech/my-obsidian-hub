import sys; sys.path.insert(0,'/home/claude/work/style')
import numpy as np, pandas as pd, json
from figstyle import *
setup(); C='/home/claude/work/cs/'; O='/home/claude/work/style/out/'
a=pd.read_csv(C+'cs_analysis_set.csv'); ei=pd.read_csv(C+'cs_ei_per_cell_ref-80_0.csv'); res=json.load(open(C+'cs_tables_part2.json'))['primary']
cm=lambda c: a.groupby(['cellkey','holding_command_mV'])[c].mean().unstack()
fi,fo=cm('f_in'),cm('f_out'); GN='#2e7d4f'
def L(h): return '0' if h==0 else ('+' if h>0 else '−')+str(abs(h))
H=[-80,-70,-60,-50,-40,-30,0,10]; xl='Command holding (mV)'
# ---- Option A
fig,ax=plt.subplots(figsize=(7.4,4.9))
for f,c,n in ((fi,BL,'Inward events (EPSC)'),(fo,RD,'Outward events (IPSC)')):
    m=np.array([f[h].dropna().median() for h in H]); lo=np.array([f[h].dropna().quantile(.25) for h in H]); hi=np.array([f[h].dropna().quantile(.75) for h in H]); x=np.arange(len(H))
    ax.fill_between(x,lo,hi,color=c,alpha=.13,lw=0); ax.plot(x,m,'-o',color=c,lw=2.2,ms=5,label=n)
ax.set_xticks(range(len(H))); ax.set_xticklabels([L(h) for h in H]); ax.set_xlabel(xl); ax.set_ylabel('Event frequency per cell (Hz)')
ax.axvline(4,color=GN,lw=1.3,ls='--'); ax.text(4.08,4.6,'−40 mV',color=GN,fontsize=9,va='top')
ax.annotate('EPSCs isolated',xy=(0,0.77),xytext=(-0.2,2.6),color=BL,fontsize=9,arrowprops=dict(arrowstyle='-',color=BL,lw=.8))
ax.annotate('IPSCs isolated',xy=(6.5,3.0),xytext=(5.0,5.4),color=RD,fontsize=9,arrowprops=dict(arrowstyle='-',color=RD,lw=.8))
ax.set_ylim(0,6.3); ax.legend(frameon=False,loc='upper left',bbox_to_anchor=(0.0,1.0)); ax.set_title('A  Both event types across holdings (median and IQR of cells)',loc='left')
foot(fig,'Cell-level medians; band = IQR. n = 7 to 16 cells per holding.'); save(fig,O+'Cs_OptA.png')
# ---- Option B
P={r['h']:r for r in res['purity']}; ref=(fi[-80]/(fi[-80]+fo[0])).dropna(); refm=ref.median()*100
fig,ax=plt.subplots(figsize=(7.4,5.2)); x=np.arange(len(H)); e=[P[h]['pctE'] for h in H]; i=[100-v for v in e]
ax.bar(x,e,color=BL,width=.7,label='Inward (EPSC)'); ax.bar(x,i,bottom=e,color=RD,width=.7,label='Outward (IPSC)')
for k,h in enumerate(H): ax.text(k,e[k]/2 if e[k]>12 else e[k]+3,f'{e[k]:.0f}%',ha='center',va='center',color='white' if e[k]>12 else BL,fontsize=9,fontweight='bold')
ax.axhline(refm,color='black',lw=1.6,ls='--',label=f'Inward share expected from the separate measurements ({refm:.0f}%)')
ax.set_xticks(x); ax.set_xticklabels([L(h) for h in H]); ax.set_xlabel(xl); ax.set_ylabel('Share of events (%)'); ax.set_ylim(0,100)
h_,l_=ax.get_legend_handles_labels(); fig.legend(h_,l_,frameon=False,loc='upper left',bbox_to_anchor=(0.01,0.955),ncol=1)
fig.suptitle('B  What share of events is inward?',x=0.01,ha='left',fontsize=9.5,y=0.99)
foot(fig,'Median over cells of each cell\'s inward share. Dashed line: E at −80 mV and I at 0 mV, same cell (n=%d).'%len(ref)); save(fig,O+'Cs_OptB.png',rect=(0,.04,1,.86))
# ---- Option C: fold per cell, full + zoom (linear)
hs=[-50,-40,-30]; g=[(L(h),ei[ei.h==h].fold.values,GN if h==-40 else GY) for h in hs]
fig,axs=plt.subplots(1,2,figsize=(10.4,5.0)); yl='E/I ratio here ÷ E/I ratio measured separately (fold)'
box_panel(axs[0],g,full_ylim(g),yl,'C  How far off is the ratio at each holding? Full scale',xlabel=xl,hline=1,band=(1/1.5,1.5,GN))
box_panel(axs[1],g,(0,6),yl,'D  Zoom to 0-6 fold',xlabel=xl,hline=1,band=(1/1.5,1.5,GN),zoom=True)
legend(axs[1],'upper right',extra=[Line2D([0],[0],color='none',label='↑k = k points above the zoom axis')])
foot(fig,'One dot = one cell; band = within ×1.5 of 1 (ratio preserved). Separate = E at −80 mV, I at 0 mV, same cell. Right: same data on a narrower scale.'); save(fig,O+'Cs_OptC.png')
# ---- Option D
fig,ax=plt.subplots(figsize=(7.0,4.9))
for col,c,n,off in (('rE',BL,'EPSC frequency vs. that at −80 mV',-.04),('rI',RD,'IPSC frequency vs. that at 0 mV',.04)):
    m=[ei[ei.h==h][col].median()*100 for h in hs]; q1=[ei[ei.h==h][col].quantile(.25)*100 for h in hs]; q3=[ei[ei.h==h][col].quantile(.75)*100 for h in hs]
    xx=np.arange(3)+off; ax.errorbar(xx,m,yerr=[np.array(m)-q1,np.array(q3)-m],color=c,marker='o',ms=7,lw=2.2,capsize=4,label=n)
    for k,(xi,mi) in enumerate(zip(xx,m)):
        dy={('rE',0):0,('rE',1):14,('rE',2):13,('rI',0):-11,('rI',1):-11,('rI',2):0}[(col,k)]; dx={'rE':-.1,'rI':.1}[col]
        if col=='rE' and k==1: dx=-.13
        if col=='rI' and k==1: dx=.13
        ax.text(xi+dx,mi+dy,f'{mi:.0f}%',color=c,ha='right' if dx<0 else 'left',va='center',fontsize=9,fontweight='bold')
ax.axhline(100,color='black',lw=1.1,ls='--'); ax.text(-.45,103,'same as at the isolating holding',ha='left',fontsize=8.5)
ax.set_xticks(range(3)); ax.set_xticklabels([L(h) for h in hs]); ax.set_xlim(-.5,2.5); ax.set_ylim(0,135); ax.set_xlabel(xl); ax.set_ylabel('Events still detected (% of isolated-holding frequency)')
ax.legend(frameon=False,loc='upper right',bbox_to_anchor=(1,1.0)); ax.set_title('D  How many events are lost at an intermediate holding?',loc='left')
foot(fig,'Median and IQR of cells (n = 14, 12, 10).'); save(fig,O+'Cs_OptD.png')
# ---- Option E
c=np.array(res['ei']['-80/0']['cross']); fig,ax=plt.subplots(figsize=(4.8,4.9))
st=[]; box_panel(ax,[('Cs-Glu cells',c,PU)],(-55,2),'Command holding where the cell\'s E/I ratio\nmatches the separate one (mV)','E  One number per cell',stats=st,xtick_fs=9)
m=np.median(c); q1,q3=np.percentile(c,[25,75]); ax.text(.36,m,f'median {m:.0f} mV\n(IQR {q1:.0f} to {q3:.0f})'.replace('-','−'),va='center',fontsize=9,color=PU)
ax.set_xlim(-.6,1.2); foot(fig,'Interpolated between the tried holdings; one dot = one cell.'); save(fig,O+'Cs_OptE.png')
print('ok')
