"""Cs report, choice of the intermediate holding (26 Nov 2026): six more views of the same question (Options F-K).
Same per-cell set as Figure 3 and Options A-E (reference = inward at -80 mV and outward at 0 mV in the same cell)."""
import sys; sys.path.insert(0,'/home/claude/work/style'); sys.path.insert(0,'/home/claude/work/mid_opts')
from cs_mid_data import *
from figstyle import *
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D
setup(); O='/home/claude/work/mid_opts/'; GN='#2e7d4f'; xl='Command holding (mV)'
def stars(p): return '***' if p<.001 else '**' if p<.01 else '*' if p<.05 else 'ns'
def wp(x,y=None):
    try: return stats.wilcoxon(x) .pvalue if y is None else stats.wilcoxon(x,y).pvalue
    except Exception: return np.nan
import textwrap
def WF(fig,t,w=150): foot(fig,'\n'.join(textwrap.fill(p,w) for p in t.split('\n')),.004)
def L(h): return '0' if h==0 else ('+' if h>0 else '−')+str(abs(h))
HS=[-50,-40,-30]; piv=lambda c: pc.pivot(index='cell',columns='h',values=c)
fold,dsh,aerr,labs=piv('fold'),piv('dshare'),piv('abserr'),None
# ============ F: heatmap, one row per cell
order=base.sort_values('ratio_ref').index.tolist(); cols=[-60,-50,-40,-30]
fig,ax=plt.subplots(figsize=(8.2,6.4)); CB,CR,CG='#8fb0d3','#e3a58a','#6fae85'
for i,c in enumerate(order):
    for j,h in enumerate(cols):
        v=fold.at[c,h] if (c in fold.index and h in fold.columns) else np.nan
        if np.isnan(v): ax.add_patch(Rectangle((j-.5,i-.5),1,1,fc='white',ec='#ddd',lw=.6,hatch='////',zorder=1)); continue
        col=CG if 1/1.5<v<1.5 else CB if v>=1.5 else CR
        ax.add_patch(Rectangle((j-.5,i-.5),1,1,fc=col,ec='white',lw=1.4,zorder=1))
        ax.text(j,i,('×%.0f'%v if v>=10 else '×%.1f'%v if v>=1 else '×%.2f'%v).replace('×','× '),ha='center',va='center',fontsize=8.2,color='#111')
ax.set_xlim(-.5,len(cols)-.5); ax.set_ylim(len(order)-.5,-.5); ax.set_xticks(range(len(cols))); ax.set_xticklabels([L(h) for h in cols]); ax.xaxis.tick_top(); ax.xaxis.set_label_position('top')
ax.set_xlabel(xl); ax.set_yticks(range(len(order))); ax.set_yticklabels([short(c) for c in order],fontsize=8.2)
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0)
for j,h in enumerate(cols):
    v=fold[h].dropna() if h in fold.columns else []; k=int(((v>1/1.5)&(v<1.5)).sum()); ax.text(j,len(order)-.25,f'{k} of {len(v)}',ha='center',va='top',fontsize=9,fontweight='bold',color=GN)
ax.text(-.5,len(order)-.25,'within × 1.5:',ha='right',va='top',fontsize=8.5,color=GN)
fig.legend(handles=[Patch(fc=CG,label='within × 1.5 of the separate E/I ratio'),Patch(fc=CB,label='too inward (EPSC-heavy), more than × 1.5'),Patch(fc=CR,label='too outward (IPSC-heavy), less than × 0.67'),Patch(fc='white',ec='#aaa',hatch='////',label='not recorded at this holding')],frameon=False,loc='lower center',ncol=2,fontsize=8,bbox_to_anchor=(.5,.045))
ax.set_title('F  Every cell at every intermediate holding',loc='left',pad=32)
WF(fig,'Number in each square = E/I frequency ratio in that cell at that holding divided by the ratio from the same cell at −80 mV (inward) and 0 mV (outward). Rows sorted by the cell\'s separate E/I ratio.'); save(fig,O+'Cs_OptF.png',rect=(0,.1 if 'Cs_OptF' in ('Cs_OptF','Cs_OptH') else .09,1,1))
# ============ G: within fold bars + cell votes
fig,axs=plt.subplots(1,2,figsize=(10.6,4.9),gridspec_kw=dict(width_ratios=[1.25,1]))
ax=axs[0]; x=np.arange(len(cols)); w=.36
for off,lim,col,lab in ((-w/2,1.5,GN,'within × 1.5'),(w/2,2.0,'#a9cfb8','within × 2')):
    for j,h in enumerate(cols):
        v=fold[h].dropna(); k=int(((v>1/lim)&(v<lim)).sum()); n=len(v)
        ax.bar(j+off,100*k/n,width=w*.95,color=col,label=lab if j==0 else None); ax.text(j+off,100*k/n+2,f'{k}/{n}',ha='center',fontsize=8.5,color='#222')
ax.set_xticks(x); ax.set_xticklabels([L(h) for h in cols]); ax.set_xlabel(xl); ax.set_ylabel('Cells whose E/I ratio is preserved (%)'); ax.set_ylim(0,100); ax.legend(frameon=False,loc='upper left',fontsize=8.5)
ax.set_title('G1  Share of cells that keep the ratio',loc='left'); ax.text(0,-9,'(3 cells)',ha='center',fontsize=7.5,color='#666')
# votes: cells with -50,-40,-30 all recorded
full=[c for c in base.index if all((c in fold.index) and (h in fold.columns) and not np.isnan(fold.at[c,h]) for h in HS)]
votes={h:0 for h in HS}
for c in full: votes[min(HS,key=lambda h: abs(np.log2(fold.at[c,h])))]+=1
ax=axs[1]; ax.bar(range(3),[votes[h] for h in HS],color=[GY if h!=-40 else GN for h in HS],width=.6)
for j,h in enumerate(HS): ax.text(j,votes[h]+.15,str(votes[h]),ha='center',fontsize=10,fontweight='bold')
ax.set_xticks(range(3)); ax.set_xticklabels([L(h) for h in HS]); ax.set_xlabel(xl); ax.set_ylabel('Number of cells'); ax.set_ylim(0,max(votes.values())+1.5); ax.set_yticks(range(0,max(votes.values())+2))
ax.set_title(f'G2  Which holding was closest, cell by cell (n = {len(full)})',loc='left')
WF(fig,'G1: ratio in the cell at that holding ÷ separately measured ratio (same cell). G2: only cells recorded at −50, −40 and −30 mV; each cell votes for the holding whose ratio is closest to its separate ratio (fold distance).'); save(fig,O+'Cs_OptG.png',rect=(0,.1 if 'Cs_OptG' in ('Cs_OptF','Cs_OptH') else .09,1,1))
print('votes',votes,len(full),'cells',[short(c) for c in full])
# ============ H: agreement scatter
fig,axs=plt.subplots(2,3,figsize=(12.4,8.4))
for r,(ymax,zm) in enumerate(((9,False),(3,True))):
    for j,h in enumerate(HS):
        ax=axs[r,j]; t=pc[pc.h==h]; xs=np.linspace(0,1.55,50)
        ax.fill_between(xs,xs/1.5,xs*1.5,color=GN,alpha=.13,lw=0); ax.plot(xs,xs,color=GN,lw=1.3)
        up=int((t.ratio_h>ymax).sum()); ax.scatter(t.ratio_ref,t.ratio_h,s=26,color=BL if h!=-40 else GN,edgecolor='white',lw=.5,zorder=3)
        ax.set_xlim(0,1.55); ax.set_ylim(0,ymax); ax.set_xlabel('E/I ratio measured separately\n(inward at −80 mV ÷ outward at 0 mV)'); ax.set_ylabel('E/I ratio at the intermediate holding')
        k=int(((t.fold>1/1.5)&(t.fold<1.5)).sum()); ax.set_title(f'{L(h)} mV{"" if r==0 else ", zoom"}   ({k} of {len(t)} within × 1.5)',loc='left',fontsize=9.5)
        if zm and up: ax.text(.97,.97,f'↑{up}',transform=ax.transAxes,ha='right',va='top',fontsize=8.5,color=BL if h!=-40 else GN,fontweight='bold')
fig.legend(handles=[Line2D([0],[0],color=GN,lw=1.3,label='same ratio as measured separately'),Patch(fc=GN,alpha=.13,label='within × 1.5')],frameon=False,loc='lower center',ncol=2,fontsize=8.5,bbox_to_anchor=(.5,.045))
WF(fig,'One dot = one cell. A dot on the line means the single holding reproduces what two separate holdings gave. Lower row: same data on a narrower scale (↑k = k cells above the axis).'); save(fig,O+'Cs_OptH.png',rect=(0,.1 if 'Cs_OptH' in ('Cs_OptF','Cs_OptH') else .09,1,1))
# ============ I: difference in inward share + absolute error
fig,axs=plt.subplots(1,2,figsize=(10.8,5.2)); gr=[(L(h),dsh[h].dropna().values,GN if h==-40 else GY) for h in HS]
box_panel(axs[0],gr,(-100,100),'Inward share of events at this holding minus the\nseparately measured share (percentage points)','I1  Over or under the separate measurement?',xlabel=xl,zero=True)
for j,h in enumerate(HS):
    v=dsh[h].dropna(); p=wp(v.values); axs[0].text(j,95,stars(p)+f'\np = {p:.3f}' if p>=.001 else stars(p)+'\np < 0.001',ha='center',va='top',fontsize=8.5,color='#222' if p<.05 else '#777')
gr2=[(L(h),aerr[h].dropna().values,GN if h==-40 else GY) for h in HS]
box_panel(axs[1],gr2,(0,100),'Distance from the separate measurement\n(percentage points of inward share)','I2  How far off is each holding?',xlabel=xl)
def brk(ax,j1,j2,y,p):
    ax.plot([j1,j1,j2,j2],[y-1.5,y,y,y-1.5],color='#222',lw=.9); ax.text((j1+j2)/2,y+.5,stars(p),ha='center',va='bottom',fontsize=10)
res={}
for (h1,h2,y) in ((-40,-50,90),(-40,-30,80)):
    t=pd.DataFrame({'a':aerr[h1],'b':aerr[h2]}).dropna(); p=wp(t.a.values,t.b.values); res[(h1,h2)]=(len(t),t.a.median(),t.b.median(),p)
    brk(axs[1],HS.index(h1),HS.index(h2),y,p)
WF(fig,'One dot = one cell. I1: Wilcoxon signed-rank test of each holding against 0. I2: paired Wilcoxon between −40 mV and its neighbours (cells with both). Inward share = inward ÷ (inward + outward) frequency. * p < 0.05, ** p < 0.01, *** p < 0.001, ns; no correction for multiple comparisons.'); save(fig,O+'Cs_OptI.png',rect=(0,.1 if 'Cs_OptI' in ('Cs_OptF','Cs_OptH') else .09,1,1))
print('abs err paired',res); print({h:(len(dsh[h].dropna()),dsh[h].median(),wp(dsh[h].dropna().values)) for h in HS})
# ============ J: per-cell inward-share curves + ECDF of matching holding
fig,axs=plt.subplots(1,2,figsize=(11.6,5.2),gridspec_kw=dict(width_ratios=[1.35,1])); ax=axs[0]
Hh=[-80,-70,-60,-50,-40,-30,0]; shr=(fi/(fi+fo))*100
for c in base.index:
    s=shr.loc[c].reindex(Hh).dropna(); ax.plot(s.index,s.values,color='#9aa3ad',lw=.9,alpha=.85,marker='o',ms=2.5,mec='none',zorder=1)
    if c in cross: ax.scatter(cross[c],base.at[c,'share_ref']*100,marker='D',s=36,color=GN,edgecolor='white',lw=.6,zorder=4)
med=shr.loc[base.index].reindex(columns=Hh).median(); ax.plot(Hh,med.values,color=BL,lw=2.6,marker='o',ms=5,zorder=3)
ax.axvline(-40,color=GN,lw=1.2,ls='--',zorder=0); ax.set_xlim(-84,4); ax.set_ylim(0,100); ax.set_xticks(Hh); ax.set_xticklabels([L(h) for h in Hh]); ax.set_xlabel(xl); ax.set_ylabel('Inward events (% of all events) in the cell')
ax.set_title('J1  Each cell\'s inward share across holdings',loc='left')
ax.legend(handles=[Line2D([0],[0],color='#9aa3ad',lw=1,marker='o',ms=3,label='One cell'),Line2D([0],[0],color=BL,lw=2.6,marker='o',ms=5,label='Median of cells'),Line2D([0],[0],color='none',marker='D',mfc=GN,mec='white',ms=7,label='Where the curve crosses the cell\'s separate share')],frameon=False,loc='upper right',fontsize=8)
cv=np.sort(np.array(list(cross.values()))); ax=axs[1]; ax.step(np.r_[-60,cv,0],np.r_[0,np.arange(1,len(cv)+1)/len(cv)*100,100],where='post',color=GN,lw=2.2)
ax.scatter(cv,np.arange(1,len(cv)+1)/len(cv)*100,s=22,color=GN,zorder=3); m=np.median(cv); q1,q3=np.percentile(cv,[25,75])
ax.axvspan(q1,q3,color=GN,alpha=.1,lw=0); ax.axvline(m,color=GN,lw=1.2,ls='--'); ax.axhline(50,color='#999',lw=.8,ls=':')
ax.text(m+1,10,f'median {m:.0f} mV\nIQR {q1:.0f} to {q3:.0f}'.replace('-','−'),fontsize=9,color=GN); ax.set_xlim(-50,0); ax.set_ylim(0,100); ax.set_xlabel('Command holding where the cell matches its separate E/I ratio (mV)'); ax.set_ylabel('Cells that have matched at or before this holding (%)')
ax.set_title(f'J2  Cumulative distribution (n = {len(cv)} cells)',loc='left')
WF(fig,'J1: the diamond is where each cell\'s curve passes through the inward share it had from the separate −80 mV and 0 mV recordings (interpolated; same rule as Option E). J2: ordered matching holdings, one dot per cell.'); save(fig,O+'Cs_OptJ.png',rect=(0,.1 if 'Cs_OptJ' in ('Cs_OptF','Cs_OptH') else .09,1,1))
print('cross q',m,q1,q3,len(cv),'<= -40:',(cv<=-39.5).sum())
# ============ K: bootstrap
rng=np.random.default_rng(7); cells=list(base.index); B=10000; win={h:0 for h in HS}; meds={h:[] for h in HS}; ok=0
A=aerr.reindex(cells)[HS].values
for b in range(B):
    idx=rng.integers(0,len(cells),len(cells)); m_=[]
    for k in range(3):
        v=A[idx,k]; v=v[~np.isnan(v)]; m_.append(np.median(v) if len(v)>=3 else np.nan); meds[HS[k]].append(m_[-1])
    if not np.isnan(m_).any(): win[HS[int(np.argmin(m_))]]+=1; ok+=1
fig,axs=plt.subplots(1,2,figsize=(10.4,4.9))
ax=axs[0]; obs=[np.nanmedian(aerr[h].dropna()) for h in HS]; lo=[np.nanpercentile(meds[h],2.5) for h in HS]; hi=[np.nanpercentile(meds[h],97.5) for h in HS]
for j,h in enumerate(HS):
    c=GN if h==-40 else GY; ax.errorbar(j,obs[j],yerr=[[obs[j]-lo[j]],[hi[j]-obs[j]]],color=c,marker='o',ms=8,lw=2,capsize=5); ax.text(j+.12,obs[j],f'{obs[j]:.0f}',va='center',fontsize=9,color=c)
ax.set_xticks(range(3)); ax.set_xticklabels([L(h) for h in HS]); ax.set_xlim(-.6,2.6); ax.set_ylim(0,60); ax.set_xlabel(xl); ax.set_ylabel('Median distance from the separate measurement\n(percentage points of inward share)'); ax.set_title('K1  Typical error, with 95% bootstrap interval',loc='left')
ax=axs[1]; pr=[100*win[h]/ok for h in HS]; ax.bar(range(3),pr,color=[GN if h==-40 else GY for h in HS],width=.6)
for j,v in enumerate(pr): ax.text(j,v+1.5,f'{v:.0f}%',ha='center',fontsize=10,fontweight='bold')
ax.set_xticks(range(3)); ax.set_xticklabels([L(h) for h in HS]); ax.set_xlabel(xl); ax.set_ylabel('Resamples in which this holding is the closest (%)'); ax.set_ylim(0,100); ax.set_title(f'K2  How often is each holding the best? ({B:,} resamples)',loc='left')
WF(fig,'Cells (n = 14) were drawn with replacement; in each draw the median distance from the separate measurement was computed for each holding, using the drawn cells recorded at that holding (at least 3). Draws with fewer than 3 cells at a holding are skipped.'); save(fig,O+'Cs_OptK.png',rect=(0,.1 if 'Cs_OptK' in ('Cs_OptF','Cs_OptH') else .09,1,1))
print('bootstrap win',{h:round(100*win[h]/ok,1) for h in HS},'ok',ok,'CI',list(zip(HS,np.round(obs,1),np.round(lo,1),np.round(hi,1))))
# robustness: fold-based distance
print('median |log2 fold|',{h:round(float(np.nanmedian(pc[pc.h==h].absl2)),2) for h in HS})
