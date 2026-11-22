import pandas as pd, numpy as np, json, itertools
from scipy import stats
D=pd.read_csv('lab_cells_long.csv')
rng=np.random.default_rng(2026)
def grp(era_src,measure):
    return D[(D.measure==measure)&(D.source.isin(era_src))].freq.values
G={
 'E':{'Paper DIV9':['Paper DIV9 control'],'Paper DIV21':['Paper DIV21 control'],
      'Current lab (pooled)':["Sam's processed EPSC, control cells (his pipeline)",'Advaith EPSC, Sam experiment (miniML)'],
      'Cs-gluconate (-80 mV)':['Cs-gluconate, E at -80 mV (miniML)']},
 'I':{'Paper DIV9':['Paper DIV9 control'],'Paper DIV21':['Paper DIV21 control'],
      'Current lab (pooled)':['Helena KCl IPSC (miniML)','Advaith KCl IPSC, Sam experiment (miniML)'],
      'Cs-gluconate (0 mV)':['Cs-gluconate, I at 0 mV (miniML)']}}
out={'desc':{},'tests':{},'ratio':{},'cs_within':{},'cs_mid':{}}
def desc(v):
    v=np.asarray(v,float); q=np.percentile(v,[25,50,75])
    return dict(n=len(v),median=q[1],q1=q[0],q3=q[2],mean=v.mean(),sem=v.std(ddof=1)/np.sqrt(len(v)),min=v.min(),max=v.max())
for meas in 'EI':
    for g,s in G[meas].items(): out['desc'][f'{meas}|{g}']=desc(grp(s,meas))
# sources separately too
for (era,src,meas),g in D.groupby(['era','source','measure']): out['desc'][f'{meas}|src|{src}']=desc(g.freq.values)
# tests
for meas in 'EI':
    keys=list(G[meas]); res=[]
    for a,b in itertools.combinations(keys,2):
        x,y=grp(G[meas][a],meas),grp(G[meas][b],meas)
        u,p=stats.mannwhitneyu(x,y,alternative='two-sided')
        res.append(dict(a=a,b=b,p=p,fold=np.median(x)/np.median(y)))
    # Holm
    ps=np.array([r['p'] for r in res]); o=np.argsort(ps); adj=np.empty(len(ps)); run=0
    for rank,i in enumerate(o): run=max(run,(len(ps)-rank)*ps[i]); adj[i]=min(run,1)
    for r,a in zip(res,adj): r['p_holm']=a
    out['tests'][meas]=res
# E/I individually measured (group level; E and I from different cells): ratio of means and medians, independent bootstrap
def boot_ratio(e,i,B=10000):
    rm,rmed=[],[]
    for _ in range(B):
        ee=rng.choice(e,len(e)); ii=rng.choice(i,len(i))
        rm.append(ee.mean()/ii.mean()); rmed.append(np.median(ee)/np.median(ii))
    return dict(ratio_mean=e.mean()/i.mean(),ci_mean=np.percentile(rm,[2.5,97.5]).tolist(),ratio_med=np.median(e)/np.median(i),ci_med=np.percentile(rmed,[2.5,97.5]).tolist(),nE=len(e),nI=len(i))
pairs={'Paper DIV9':('Paper DIV9','Paper DIV9'),'Paper DIV21':('Paper DIV21','Paper DIV21'),'Current lab (pooled)':('Current lab (pooled)','Current lab (pooled)'),'Cs-gluconate (-80 / 0 mV, separate cells)':('Cs-gluconate (-80 mV)','Cs-gluconate (0 mV)')}
for k,(e,i) in pairs.items(): out['ratio'][k]=boot_ratio(grp(G['E'][e],'E'),grp(G['I'][i],'I'))
# Current lab, miniML only (excluding Sam's own pipeline) as sensitivity: Advaith EPSC (n=3) vs KCl pooled
# Cs same-cell E/I
cs=pd.read_csv('/home/claude/work/cs/cs_cell_holding.csv')
E=cs[cs.holding_command_mV==-80].set_index('cellkey').f_in; I=cs[cs.holding_command_mV==0].set_index('cellkey').f_out
both=pd.concat([E,I],axis=1,keys=['E','I']).dropna(); both=both[(both.E>0)&(both.I>0)]; both['ratio']=both.E/both.I
r=both.ratio.values; bs=[np.median(rng.choice(r,len(r))) for _ in range(10000)]
out['cs_within']=dict(n=len(r),median=float(np.median(r)),q1=float(np.percentile(r,25)),q3=float(np.percentile(r,75)),mean=float(r.mean()),ci_med=np.percentile(bs,[2.5,97.5]).tolist(),min=float(r.min()),max=float(r.max()))
both.to_csv('cs_E_I_same_cell.csv')
# Cs intermediate: per-cell E/I (inward/outward) at -40 and reference
ei=pd.read_csv('/home/claude/work/cs/cs_ei_per_cell_ref-80_0.csv')
mid={}
for h in (-50,-40,-30):
    t=ei[ei.h==h]
    if len(t)<3: continue
    w=stats.wilcoxon(np.log(t.ratio_h),np.log(t.ratio_ref)) if len(t)>=6 else None
    mid[h]=dict(n=len(t),ratio_h_med=float(t.ratio_h.median()),ratio_h_q1=float(t.ratio_h.quantile(.25)),ratio_h_q3=float(t.ratio_h.quantile(.75)),
        ratio_ref_med=float(t.ratio_ref.median()),ratio_ref_q1=float(t.ratio_ref.quantile(.25)),ratio_ref_q3=float(t.ratio_ref.quantile(.75)),
        fold_med=float(t.fold.median()),fold_q1=float(t.fold.quantile(.25)),fold_q3=float(t.fold.quantile(.75)),p=float(w.pvalue) if w else None,
        n_within_1p5=int(((t.fold>=1/1.5)&(t.fold<=1.5)).sum()),
        group_ratio_means=float(t.fi.mean()/t.fo.mean()))
out['cs_mid']=mid
# sensitivity: second detector
sens={}
for meas in 'EI':
    for src,g in D[D.measure==meas].groupby('source'):
        v=g.freq_d2.dropna().values
        if len(v): sens[f'{meas}|{src}']=dict(n=len(v),median=float(np.median(v)),q1=float(np.percentile(v,25)),q3=float(np.percentile(v,75)),mean=float(v.mean()),n_minimal=len(g))
out['d2']=sens
json.dump(out,open('lab_stats.json','w'),indent=1,default=float)
for k,v in out['desc'].items(): print(k,{a:round(b,3) for a,b in v.items()})
print(json.dumps(out['tests'],indent=0,default=float)[:1800]); print(json.dumps(out['ratio'],default=float)); print(out['cs_within']); print(json.dumps(out['cs_mid'],default=float)); print(json.dumps(out['d2'],default=float))
