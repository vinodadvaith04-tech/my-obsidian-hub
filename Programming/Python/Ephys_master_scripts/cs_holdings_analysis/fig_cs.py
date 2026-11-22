import numpy as np, pandas as pd, matplotlib, json
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
a=pd.read_csv('cs_analysis_set.csv')
H=[-80,-70,-60,-50,-40,-30,0,10]
cm=lambda c: a.groupby(['cellkey','holding_command_mV'])[c].mean().unstack()
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
BOX=[]
def lab(h): return f'{h:+d}'.replace('+0','0') if h!=0 else '0'
def draw(ax,data,order,col,ylim,ylabel,title,tag,log=False,unit='',zero=False,hline=None):
    rng=np.random.default_rng(11); lo,hi=ylim
    for k,key in enumerate(order):
        v=np.asarray(data[key],float); v=v[~np.isnan(v)]
        if len(v)==0: continue
        q1,med,q3=np.percentile(v,[25,50,75]); mn,mx,mu=v.min(),v.max(),v.mean(); w=.3
        BOX.append((tag,key,len(v),mu,med,q1,q3,mn,mx))
        cl=lambda y: min(max(y,lo),hi)
        ax.add_patch(plt.Rectangle((k-w,cl(q1)),2*w,max(cl(q3)-cl(q1),0),fc=col,alpha=.2,ec=col,lw=1.2,zorder=1))
        if lo<=med<=hi: ax.plot([k-w,k+w],[med,med],color=col,lw=2.4,zorder=3)
        ax.plot([k,k],[cl(mn),cl(q1)],color=col,lw=1.2,zorder=1); ax.plot([k,k],[cl(q3),cl(mx)],color=col,lw=1.2,zorder=1)
        if lo<=mn<=hi: ax.plot([k-w/2,k+w/2],[mn,mn],color=col,lw=1.2,zorder=1)
        if lo<=mx<=hi: ax.plot([k-w/2,k+w/2],[mx,mx],color=col,lw=1.2,zorder=1)
        if lo<=mu<=hi: ax.scatter([k],[mu],marker='D',s=40,fc='white',ec=col,lw=1.6,zorder=5)
        xs=k+rng.uniform(-.15,.15,len(v)); ins=(v>=lo)&(v<=hi)
        ax.scatter(xs[ins],v[ins],s=20,color=col,alpha=.85,lw=.4,ec='white',zorder=4)
        for x,y in zip(xs[~ins],v[~ins]):
            up=y>hi
            if log: ypos=hi*.93 if up else lo*1.08
            else: ypos=hi-(hi-lo)*.015 if up else lo+(hi-lo)*.015
            ax.scatter([x],[ypos],marker='^' if up else 'v',s=34,color=col,zorder=4)
        ax.text(k,-0.035,f'n={len(v)}',ha='center',va='top',fontsize=8,transform=ax.get_xaxis_transform(),color='#333')
    ax.set_xticks(range(len(order))); ax.set_xticklabels([str(o) for o in order]); ax.tick_params(axis='x',pad=14)
    ax.set_xlim(-.6,len(order)-.4); ax.set_ylim(lo,hi)
    if log: ax.set_yscale('log')
    ax.set_xlabel('Command holding (mV)',labelpad=2); ax.set_ylabel(ylabel); ax.set_title(title,loc='left',fontsize=9.5)
    if zero: ax.axhline(0,color='#999',lw=.7,zorder=0)
    if hline is not None: ax.axhline(hline,color='#999',lw=.9,ls='--',zorder=0)
def legend(ax,loc='upper left',anchor=None):
    ax.legend(handles=[Line2D([0],[0],color='#555',lw=2.4,label='Median'),Line2D([0],[0],marker='D',color='w',mfc='white',mec='#555',mew=1.6,ms=6,label='Mean'),
      Line2D([0],[0],marker='s',color='w',mfc='#ccc',mec='#555',ms=8,label='Box = IQR'),Line2D([0],[0],color='#555',lw=1.2,label='Whiskers = min, max'),
      Line2D([0],[0],marker='o',color='w',mfc='#555',ms=5,label='One cell'),Line2D([0],[0],marker='^',color='w',mfc='#555',ms=6,label='Outside axis')],frameon=False,fontsize=7.5,loc=loc,bbox_to_anchor=anchor)
BL='#1f4e79'; RD='#b5451b'; GR='#4a7c59'
# ---------------- Figure 1: noise, drift, holding current
def cell_vals(c): return {h:cm(c)[h].dropna().values if h in cm(c).columns else np.array([]) for h in H}
fig,axs=plt.subplots(1,3,figsize=(14,4.8))
draw(axs[0],cell_vals('noise_lp3k_pA'),H,BL,(0,9),'Noise per cell (pA)','A  Noise (MAD-SD, 3 kHz low-pass)','noise')
draw(axs[1],cell_vals('drift_total_pA'),H,RD,(-150,100),'Drift per cell (pA over the 2 min)','B  Drift over 2 min','drift',zero=True)
draw(axs[2],cell_vals('Ihold_median_pA'),H,GR,(-450,250),'Holding current per cell (pA)','C  Holding current','ihold',zero=True)
legend(axs[0],'upper center',(0.55,1.0))
fig.text(.01,.005,'Each dot is one cell (mean of its recordings at that command holding). -20 mV and +15 mV (1 cell each) are in the tables but not drawn. Drift is signed: negative = current falling over the 2 min.',fontsize=7.5,color='#444')
fig.tight_layout(rect=(0,.03,1,1)); fig.savefig('Figure1_cs_noise_drift_holding_current.png',dpi=200); plt.close(fig)
# ---------------- Figure 2: event frequency per cell, E and I
fi,fo=cm('f_in'),cm('f_out')
D_in={h:fi[h].dropna().values for h in H}; D_out={h:fo[h].dropna().values for h in H}
fig,axs=plt.subplots(2,2,figsize=(12.6,8.6))
draw(axs[0,0],D_in,H,BL,(0,3.3),'miniML frequency per cell (Hz)','A  Inward events (EPSC), full range','in_full')
draw(axs[0,1],D_out,H,RD,(0,12),'miniML frequency per cell (Hz)','B  Outward events (IPSC), full range','out_full')
draw(axs[1,0],D_in,H,BL,(0,0.4),'miniML frequency per cell (Hz)','C  Inward events, zoom to 0-0.4 Hz','in_zoom')
draw(axs[1,1],D_out,H,RD,(0,0.4),'miniML frequency per cell (Hz)','D  Outward events, zoom to 0-0.4 Hz','out_zoom')
legend(axs[0,1],'upper left',(0.0,1.0))
for ax in (axs[0,0],axs[1,1]): ax.axvspan(-.5,1.5,color=BL,alpha=.06,zorder=0)
for ax in (axs[0,1],axs[1,0]): ax.axvspan(5.5,7.5,color=RD,alpha=.06,zorder=0)
fig.text(.01,.005,'Shaded: candidate isolating holdings. Blue: -80, -70 for EPSCs (A, and D = residual outward events). Red: 0, +10 for IPSCs (B, and C = residual inward events).\nEach dot is one cell (mean of its recordings at that holding). -20 and +15 mV (1 cell each) are not drawn.',fontsize=7.5,color='#444')
fig.tight_layout(rect=(0,.035,1,1)); fig.savefig('Figure2_cs_event_frequency_boxplots.png',dpi=200); plt.close(fig)
pd.DataFrame(BOX,columns=['panel','holding','n_cells','mean','median','q1','q3','min','max']).round(4).to_csv('cs_box_stats.csv',index=False)
# ---------------- Figure 3: E/I preservation
res=json.load(open('cs_tables_part2.json'))['primary']['ei']
qi,qo=cm('q_in'),cm('q_out')
mid=[-60,-50,-40,-30]
def per_cell(e,i):
    rows=[]
    for c in fi.index:
        E0,I0=fi.at[c,e] if e in fi else np.nan, fo.at[c,i] if i in fo else np.nan
        if not(E0>0 and I0>0): continue
        for h in mid:
            x,y=fi.at[c,h] if h in fi else np.nan, fo.at[c,h] if h in fo else np.nan
            if x>0 and y>0: rows.append(dict(cell=c,h=h,E0=E0,I0=I0,fi=x,fo=y,rE=x/E0,rI=y/I0,fold=(x/E0)/(y/I0),ratio_h=x/y,ratio_ref=E0/I0))
    return pd.DataFrame(rows)
pc=per_cell(-80,0); pc.to_csv('cs_ei_per_cell_ref-80_0.csv',index=False)
fig,axs=plt.subplots(1,3,figsize=(14.6,4.9),gridspec_kw=dict(width_ratios=[1.0,1.0,.9]))
draw(axs[0],{h:pc[pc.h==h].fold.values for h in mid},mid,'#5b4a8a',(0.05,40),'Within-cell E/I ratio at h  /  reference E/I ratio (fold)','A  E/I frequency ratio relative to the separately measured one','fold',log=True,hline=1)
axs[0].set_yticks([.0625,.125,.25,.5,1,2,4,8,16,32]); axs[0].set_yticklabels(['1/16','1/8','1/4','1/2','1','2','4','8','16','32'])
axs[0].axhspan(1/1.5,1.5,color='#5b4a8a',alpha=.08,zorder=0); axs[0].set_ylim(0.0625,40)
axs[0].text(-0.55,1.03,'ratio preserved (x1/1.5 to x1.5)',fontsize=7.5,ha='left',va='bottom',color='#444')
# panel B: retention
ret=[]
for h in mid:
    t=pc[pc.h==h]; ret.append((h,len(t),np.median(t.rE)*100,*np.percentile(t.rE*100,[25,75]),np.median(t.rI)*100,*np.percentile(t.rI*100,[25,75])))
ret=pd.DataFrame(ret,columns=['h','n','rE','rE25','rE75','rI','rI25','rI75'])
ax=axs[1]; x=np.arange(len(mid))
ax.errorbar(x-.04,ret.rE,yerr=[ret.rE-ret.rE25,ret.rE75-ret.rE],color=BL,marker='o',lw=1.8,capsize=3,label='EPSC frequency (inward) / that at -80')
ax.errorbar(x+.04,ret.rI,yerr=[ret.rI-ret.rI25,ret.rI75-ret.rI],color=RD,marker='s',lw=1.8,capsize=3,label='IPSC frequency (outward) / that at 0')
for k,r in ret.iterrows(): ax.text(k,-0.035,f'n={int(r.n)}',ha='center',va='top',fontsize=8,transform=ax.get_xaxis_transform(),color='#333')
ax.set_xticks(x); ax.set_xticklabels(mid); ax.tick_params(axis='x',pad=14); ax.set_xlabel('Command holding (mV)',labelpad=2)
ax.set_ylabel('Share of the isolated-holding frequency still detected (%)'); ax.set_yscale('log'); ax.set_ylim(2,300); ax.axhline(100,color='#999',lw=.9,ls='--')
ax.set_yticks([3,10,30,100,300]); ax.set_yticklabels(['3','10','30','100','300']); ax.legend(frameon=False,fontsize=7.5,loc='upper right')
ax.set_title('B  How much of each event type is still detected',loc='left',fontsize=9.5)
# panel C: crossing holdings
refs=[('-80 / 0','E at -80, I at 0'),('-70 / 0','E at -70, I at 0'),('-80/10','E at -80, I at +10'),('-70/10','E at -70, I at +10')]
rng=np.random.default_rng(3); ax=axs[2]
keys=['-80/0','-70/0','-80/10','-70/10']; names=['-80 | 0','-70 | 0','-80 | +10','-70 | +10']
for k,key in enumerate(keys):
    c=np.array(res[key]['cross']); q1,med,q3=np.percentile(c,[25,50,75]); w=.3
    BOX.append(('cross',names[k],len(c),c.mean(),med,q1,q3,c.min(),c.max()))
    ax.add_patch(plt.Rectangle((k-w,q1),2*w,q3-q1,fc='#5b4a8a',alpha=.2,ec='#5b4a8a',lw=1.2)); ax.plot([k-w,k+w],[med,med],color='#5b4a8a',lw=2.4)
    ax.plot([k,k],[c.min(),q1],color='#5b4a8a',lw=1.2); ax.plot([k,k],[q3,c.max()],color='#5b4a8a',lw=1.2)
    ax.plot([k-w/2,k+w/2],[c.min()]*2,color='#5b4a8a',lw=1.2); ax.plot([k-w/2,k+w/2],[c.max()]*2,color='#5b4a8a',lw=1.2)
    ax.scatter([k],[c.mean()],marker='D',s=40,fc='white',ec='#5b4a8a',lw=1.6,zorder=5)
    ax.scatter(k+rng.uniform(-.15,.15,len(c)),c,s=20,color='#5b4a8a',alpha=.85,lw=.4,ec='white',zorder=4)
    ax.text(k,-0.035,f'n={len(c)}',ha='center',va='top',fontsize=8,transform=ax.get_xaxis_transform(),color='#333')
ax.set_xticks(range(4)); ax.set_xticklabels(names); ax.tick_params(axis='x',pad=14); ax.set_xlim(-.6,3.6); ax.set_ylim(-65,0)
ax.set_xlabel('Reference holdings (E | I)',labelpad=2); ax.set_ylabel('Command holding where the cell\'s E/I ratio\nequals the reference (mV, interpolated)')
ax.set_title('C  Holding that reproduces the separate E/I ratio',loc='left',fontsize=9.5)
legend(axs[0],'lower left',(0.0,0.0))
fig.text(.01,.005,'Each dot is one cell. Ratio = inward/outward event frequency. Reference = inward frequency at -80 mV (E) and outward frequency at 0 mV (I) in the same cell. Grey band in A: within x1.5.',fontsize=7.5,color='#444')
fig.tight_layout(rect=(0,.03,1,1)); fig.savefig('Figure3_cs_EI_ratio_vs_holding.png',dpi=200); plt.close(fig)
pd.DataFrame(BOX,columns=['panel','holding','n_cells','mean','median','q1','q3','min','max']).round(4).to_csv('cs_box_stats.csv',index=False)
print(ret.round(1)); print(pd.DataFrame(BOX).tail(5))
