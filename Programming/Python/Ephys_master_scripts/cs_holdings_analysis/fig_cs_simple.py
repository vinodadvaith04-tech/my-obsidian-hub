import numpy as np, pandas as pd, matplotlib, json
matplotlib.use('Agg'); import matplotlib.pyplot as plt
a=pd.read_csv('cs_analysis_set.csv'); ei=pd.read_csv('cs_ei_per_cell_ref-80_0.csv'); res=json.load(open('cs_tables_part2.json'))['primary']
cm=lambda c: a.groupby(['cellkey','holding_command_mV'])[c].mean().unstack()
fi,fo=cm('f_in'),cm('f_out')
BL='#1f4e79'; RD='#b5451b'; PU='#5b4a8a'; GN='#2e7d4f'; GY='#8a8a8a'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.titlesize':12})
def L(h): return '0' if h==0 else ('+' if h>0 else '−')+str(abs(h))
H=[-80,-70,-60,-50,-40,-30,0,10]
S={}
# ---------- Option A: median E and I frequency vs holding
fig,ax=plt.subplots(figsize=(7.4,4.8))
for f,c,n in ((fi,BL,'Inward events (EPSC)'),(fo,RD,'Outward events (IPSC)')):
    m=[];lo=[];hi=[]
    for h in H:
        v=f[h].dropna(); m.append(v.median()); lo.append(v.quantile(.25)); hi.append(v.quantile(.75))
    m,lo,hi=map(np.array,(m,lo,hi)); x=np.arange(len(H))
    ax.fill_between(x,lo,hi,color=c,alpha=.13,lw=0); ax.plot(x,m,'-o',color=c,lw=2.6,ms=6,label=n)
ax.set_xticks(range(len(H))); ax.set_xticklabels([L(h) for h in H]); ax.set_xlabel('Command holding (mV)'); ax.set_ylabel('Event frequency per cell (Hz)')
ax.axvline(4,color=GN,lw=1.4,ls='--'); ax.text(4.08,3.55,'−40 mV',color=GN,fontsize=10.5,va='top')
ax.annotate('EPSCs isolated',xy=(0,0.77),xytext=(-0.2,2.1),color=BL,fontsize=10.5,arrowprops=dict(arrowstyle='-',color=BL,lw=.8))
ax.annotate('IPSCs isolated',xy=(6.5,3.0),xytext=(5.0,4.4),color=RD,fontsize=10.5,arrowprops=dict(arrowstyle='-',color=RD,lw=.8))
ax.set_ylim(0,5.2); ax.legend(frameon=False,loc='upper center',bbox_to_anchor=(0.45,1.0)); ax.set_title('A  Both event types across holdings (median and IQR of cells)',loc='left')
fig.text(.01,.005,'Cell-level medians; band = IQR. n = 7 to 16 cells per holding.',fontsize=8.5,color='#444'); fig.tight_layout(rect=(0,.03,1,1)); fig.savefig('FigureOpt_A_E_and_I_lines.png',dpi=200); plt.close(fig)
# ---------- Option B: share of events inward vs outward, stacked, with reference
pur=res['purity']; P={r['h']:r for r in pur}
ref=(fi[-80]/(fi[-80]+fo[0])).dropna(); refm=ref.median()*100
fig,ax=plt.subplots(figsize=(7.4,5.2))
x=np.arange(len(H)); e=[P[h]['pctE'] for h in H]; i=[100-v for v in e]
ax.bar(x,e,color=BL,width=.7,label='Inward (EPSC)'); ax.bar(x,i,bottom=e,color=RD,width=.7,label='Outward (IPSC)')
for k,h in enumerate(H): ax.text(k,e[k]/2 if e[k]>12 else e[k]+3,f'{e[k]:.0f}%',ha='center',va='center',color='white' if e[k]>12 else BL,fontsize=10,fontweight='bold')
ax.axhline(refm,color='black',lw=1.8,ls='--',label=f'Inward share expected from the separate measurements ({refm:.0f}%)')
ax.set_xticks(x); ax.set_xticklabels([L(h) for h in H]); ax.set_xlabel('Command holding (mV)'); ax.set_ylabel('Share of events (%)'); ax.set_ylim(0,100)
h_,l_=ax.get_legend_handles_labels(); fig.legend(h_,l_,frameon=False,loc='upper left',bbox_to_anchor=(0.01,0.955),ncol=1,fontsize=9.5)
fig.suptitle('B  What share of events is inward?',x=0.01,ha='left',fontsize=12,y=0.99)
fig.text(.01,.005,'Median over cells of each cell\'s inward share. Dashed line: E at \u221280 mV and I at 0 mV, same cell (n=%d).'%len(ref),fontsize=8.5,color='#444')
fig.tight_layout(rect=(0,.03,1,.86)); fig.savefig('FigureOpt_B_event_share.png',dpi=200); plt.close(fig); S['ref_share']=refm
# ---------- Option C: fold per cell at -50,-40,-30 on log axis
fig,ax=plt.subplots(figsize=(6.8,5.0)); rng=np.random.default_rng(2); hs=[-50,-40,-30]
for k,h in enumerate(hs):
    v=ei[ei.h==h].fold.values; c=GN if h==-40 else GY
    ax.scatter(k+rng.uniform(-.13,.13,len(v)),v,s=40,color=c,alpha=.9,ec='white',lw=.6,zorder=3)
    ax.plot([k-.28,k+.28],[np.median(v)]*2,color=GN if h==-40 else '#444',lw=3,zorder=4)
    ax.text(k,-0.13,f'median {np.median(v):.2f}\u00d7  (n={len(v)})',ha='center',va='top',fontsize=9.5,color='#333',transform=ax.get_xaxis_transform())
ax.axhline(1,color='black',lw=1.4,ls='--'); ax.axhspan(1/1.5,1.5,color=GN,alpha=.08,lw=0); ax.set_yscale('log'); ax.set_ylim(0.04,40)
ax.set_yticks([.0625,.125,.25,.5,1,2,4,8,16,32]); ax.set_yticklabels(['1/16','1/8','1/4','1/2','1 (same)','2','4','8','16','32'])
ax.set_xticks(range(3)); ax.set_xticklabels([L(h) for h in hs]); ax.set_xlim(-.6,2.6); ax.set_xlabel('Command holding (mV)',labelpad=38); ax.set_ylabel('E/I ratio here \u00f7 E/I ratio measured separately')
ax.text(2.55,1.12,'ratio preserved',ha='right',fontsize=9.5,color='#444'); ax.text(-.55,28,'too much EPSC',fontsize=9.5,color=BL); ax.text(-.55,.065,'too much IPSC',fontsize=9.5,color=RD)
ax.set_title('C  How far off is the ratio at each holding?',loc='left')
fig.text(.01,.005,'One dot = one cell; bar = median; band = within \u00d71.5. Separate = E at \u221280 mV, I at 0 mV, same cell.',fontsize=8.5,color='#444')
fig.tight_layout(rect=(0,.03,1,1)); fig.savefig('FigureOpt_C_fold_dots.png',dpi=200); plt.close(fig)
# ---------- Option D: detection retained
fig,ax=plt.subplots(figsize=(7.0,4.8)); hs=[-50,-40,-30]
for col,c,n,off in (('rE',BL,'EPSC frequency vs. that at \u221280 mV',-.04),('rI',RD,'IPSC frequency vs. that at 0 mV',.04)):
    ei['rr']=ei['rE'] if col=='rE' else ei['rI']
    m=[ei[ei.h==h].rr.median()*100 for h in hs]; q1=[ei[ei.h==h].rr.quantile(.25)*100 for h in hs]; q3=[ei[ei.h==h].rr.quantile(.75)*100 for h in hs]
    xx=np.arange(3)+off; ax.errorbar(xx,m,yerr=[np.array(m)-q1,np.array(q3)-m],color=c,marker='o',ms=8,lw=2.6,capsize=4,label=n)
    for k,(xi,mi) in enumerate(zip(xx,m)):
        dy={('rE',0):0,('rE',1):14,('rE',2):13,('rI',0):-11,('rI',1):-11,('rI',2):0}[(col,k)]; dx={'rE':-.1,'rI':.1}[col]
        if col=='rE' and k==1: dx=-.13
        if col=='rI' and k==1: dx=.13
        ax.text(xi+dx,mi+dy,f'{mi:.0f}%',color=c,ha='right' if dx<0 else 'left',va='center',fontsize=10.5,fontweight='bold')
ax.axhline(100,color='black',lw=1.2,ls='--'); ax.text(2.45,106,'same as at the isolating holding',ha='right',fontsize=9.5)
ax.set_xticks(range(3)); ax.set_xticklabels([L(h) for h in hs]); ax.set_xlim(-.5,2.5); ax.set_ylim(0,110); ax.set_xlabel('Command holding (mV)'); ax.set_ylabel('Events still detected (% of isolated-holding frequency)')
ax.legend(frameon=False,loc='upper right',bbox_to_anchor=(1,.93),fontsize=10); ax.set_title('D  How many events are lost at an intermediate holding?',loc='left')
fig.text(.01,.005,'Median and IQR of cells (n = 14, 12, 10).',fontsize=8.5,color='#444'); fig.tight_layout(rect=(0,.03,1,1)); fig.savefig('FigureOpt_D_detection_retained.png',dpi=200); plt.close(fig)
# ---------- Option E: matching holding dot strip
c=np.array(res['ei']['-80/0']['cross']); fig,ax=plt.subplots(figsize=(4.8,4.8)); rng=np.random.default_rng(4)
ax.scatter(rng.uniform(-.12,.12,len(c)),c,s=55,color=PU,alpha=.9,ec='white',lw=.7,zorder=3); m=np.median(c); q1,q3=np.percentile(c,[25,75])
ax.add_patch(plt.Rectangle((-.3,q1),.6,q3-q1,fc=GN,alpha=.15,ec='none')); ax.plot([-.3,.3],[m,m],color=GN,lw=3.2)
ax.text(.34,m,f'median {m:.0f} mV\n(IQR {q1:.0f} to {q3:.0f})'.replace('-','−'),va='center',fontsize=10.5,color=GN)
ax.set_xlim(-.55,1.0); ax.set_ylim(-52,0); ax.set_xticks([]); ax.set_ylabel('Command holding where the cell\'s E/I ratio\nmatches the separate one (mV)'); ax.spines['bottom'].set_visible(False)
ax.set_yticks([-50,-40,-30,-20,-10,0]); ax.set_yticklabels(['−50','−40','−30','−20','−10','0'])
ax.set_title('E  One number per cell',loc='left'); ax.text(0,-51,f'n={len(c)} cells',ha='center',fontsize=9.5,color='#333')
fig.text(.01,.005,'Interpolated between the tried holdings; one dot = one cell.',fontsize=8.5,color='#444'); fig.tight_layout(rect=(0,.03,1,1)); fig.savefig('FigureOpt_E_matching_holding.png',dpi=200); plt.close(fig)
json.dump(S,open('simple_fig_vals.json','w')); print(S)
