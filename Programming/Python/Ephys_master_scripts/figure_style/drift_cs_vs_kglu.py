"""Cs-gluconate vs K-gluconate drift, matched by command holding. Drift = |change of the fitted linear trend of 1-s block medians over the
2-min recording| (drift_total_pA, as in the reports). Unit = cell (mean of its recordings at that holding)."""
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
import statsmodels.formula.api as smf
W='/home/claude/work/'
c=pd.read_csv(W+'cs/cs_analysis_set.csv'); c['internal']='Cs'
m=pd.read_csv(W+'k60/mid/mid_master_table.csv'); a=pd.read_csv(W+'k60/k60_analysis_set3.csv')
a['cellkey']=a.date.astype(str)+'_s'+a.slice.astype(str)+'_c'+a.cell.astype(str)
k=pd.concat([m,a],ignore_index=True); k['internal']='K-Glu'; k=k.drop_duplicates('new_name')
cols=['new_name','cellkey','holding_command_mV','internal','drift_total_pA','Ihold_median_pA','Ra_b']
d=pd.concat([c[cols],k[cols]],ignore_index=True)
d['ad']=d.drift_total_pA.abs(); d['aI']=d.Ihold_median_pA.abs(); d['h']=d.holding_command_mV; d.loc[d.h.between(5,15),'h']=10
cell=d.groupby(['internal','cellkey','h']).agg(ad=('ad','mean'),aI=('aI','mean'),n=('ad','size')).reset_index()
cell.to_csv('drift_cs_vs_kglu_per_cell.csv',index=False)
H=[-60,-50,-40,-30,0]; rng=np.random.default_rng(1); rows=[]
for h in H:
    x=cell[(cell.h==h)&(cell.internal=='K-Glu')].ad.values; y=cell[(cell.h==h)&(cell.internal=='Cs')].ad.values
    bs=[np.median(rng.choice(y,len(y)))/np.median(rng.choice(x,len(x))) for _ in range(5000)]
    rows.append(dict(h=h,nK=len(x),nCs=len(y),medK=np.median(x),medCs=np.median(y),ratio=np.median(y)/np.median(x),lo=np.percentile(bs,2.5),hi=np.percentile(bs,97.5),p=mannwhitneyu(y,x).pvalue,
                     n25K=(x>25).sum(),n25Cs=(y>25).sum()))
R=pd.DataFrame(rows); print(R.round(3).to_string()); R.to_csv('drift_cs_vs_kglu_by_holding.csv',index=False)
r=d[d.h.isin(H)&d.ad.notna()].copy(); r['lg']=np.log(r.ad+1); r['hc']=r.h.astype(str); r['Cs']=(r.internal=='Cs').astype(int)
mm=smf.mixedlm('lg~Cs+C(hc)',r,groups=r.cellkey).fit(); b,se=mm.params['Cs'],mm.bse['Cs']
print('mixed log(|drift|+1): Cs vs K-Glu ratio %.2f (%.2f-%.2f) p=%.3f'%(np.exp(b),np.exp(b-1.96*se),np.exp(b+1.96*se),mm.pvalues['Cs']))
r2=r[r.aI.notna()].copy(); r2['lI']=np.log(r2.aI+1); m2=smf.mixedlm('lg~Cs+C(hc)+lI',r2,groups=r2.cellkey).fit(); b,se=m2.params['Cs'],m2.bse['Cs']
print(' + |holding current| adjusted: ratio %.2f (%.2f-%.2f) p=%.3f'%(np.exp(b),np.exp(b-1.96*se),np.exp(b+1.96*se),m2.pvalues['Cs']))
print(cell[cell.h.isin(H)].groupby(['internal','h']).ad.agg(['max','mean']).round(1).T.to_string())
# sensitivity: drop Cs recordings flagged leaky or Ra not ok (the Cs report's strict set)
bad=set(c[(c.leaky==True)|(c.ra_ok==False)].new_name)
r3=r[~r.new_name.isin(bad)]
mm=smf.mixedlm('lg~Cs+C(hc)',r3,groups=r3.cellkey).fit(); b,se=mm.params['Cs'],mm.bse['Cs']
print('strict Cs (no leaky / Ra not ok; %d recordings dropped): ratio %.2f (%.2f-%.2f) p=%.3f'%(len(bad),np.exp(b),np.exp(b-1.96*se),np.exp(b+1.96*se),mm.pvalues['Cs']))
cs3=r3.groupby(['internal','cellkey','h']).ad.mean().reset_index()
for h in H:
    x=cs3[(cs3.h==h)&(cs3.internal=='K-Glu')].ad; y=cs3[(cs3.h==h)&(cs3.internal=='Cs')].ad
    print(h,len(x),len(y),round(x.median(),1),round(y.median(),1),'p=%.3f'%mannwhitneyu(y,x).pvalue)
