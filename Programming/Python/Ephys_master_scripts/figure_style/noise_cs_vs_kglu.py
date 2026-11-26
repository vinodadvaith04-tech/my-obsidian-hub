"""Cs-gluconate vs K-gluconate baseline noise, matched by command holding. Noise = peak-to-peak = 6 x MAD-SD (3 kHz low-pass),
whole recording (lab standard decided 26 Nov 2026). Unit = cell (mean of its recordings at that holding)."""
import sys; sys.path.insert(0,'/home/claude/work/style')
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
import statsmodels.formula.api as smf
W='/home/claude/work/'
c=pd.read_csv(W+'cs/cs_analysis_set.csv'); c['internal']='Cs'
m=pd.read_csv(W+'k60/mid/mid_master_table.csv'); a=pd.read_csv(W+'k60/k60_analysis_set3.csv')
a['cellkey']=a.date.astype(str)+'_s'+a.slice.astype(str)+'_c'+a.cell.astype(str)
k=pd.concat([m,a],ignore_index=True); k['internal']='K-Glu'
k=k.drop_duplicates('new_name')
cols=['new_name','cellkey','holding_command_mV','internal','noise_lp3k_pA','noise_raw_pA','Ra_b']
d=pd.concat([c[cols],k[cols]],ignore_index=True)
d['p2p']=6*d.noise_lp3k_pA; d['p2p_raw']=6*d.noise_raw_pA; d['h']=d.holding_command_mV
d.loc[d.h.between(5,15),'h']=10   # +5..+15 pooled as in the reports
cell=d.groupby(['internal','cellkey','h']).agg(p2p=('p2p','mean'),p2p_raw=('p2p_raw','mean'),Ra=('Ra_b','mean'),n=('p2p','size')).reset_index()
cell.to_csv('noise_cs_vs_kglu_per_cell.csv',index=False)
H=[-60,-50,-40,-30,0]; rng=np.random.default_rng(1)
rows=[]
for h in H:
    x=cell[(cell.h==h)&(cell.internal=='K-Glu')].p2p.values; y=cell[(cell.h==h)&(cell.internal=='Cs')].p2p.values
    p=mannwhitneyu(y,x).pvalue
    bs=[np.median(rng.choice(y,len(y)))/np.median(rng.choice(x,len(x))) for _ in range(5000)]
    rows.append(dict(h=h,nK=len(x),nCs=len(y),medK=np.median(x),medCs=np.median(y),ratio=np.median(y)/np.median(x),lo=np.percentile(bs,2.5),hi=np.percentile(bs,97.5),p=p))
R=pd.DataFrame(rows); print(R.round(3).to_string()); R.to_csv('noise_cs_vs_kglu_by_holding.csv',index=False)
# pooled: recording-level log noise, holding as categorical, random intercept per cell
r=d[d.h.isin(H)&d.p2p.notna()].copy(); r['lg']=np.log(r.p2p); r['hc']=r.h.astype(str); r['Cs']=(r.internal=='Cs').astype(int)
mm=smf.mixedlm('lg~Cs+C(hc)',r,groups=r.cellkey).fit(); b,se=mm.params['Cs'],mm.bse['Cs']
print('mixed model Cs vs K-Glu: ratio %.2f (95%% CI %.2f-%.2f) p=%.3f'%(np.exp(b),np.exp(b-1.96*se),np.exp(b+1.96*se),mm.pvalues['Cs']))
r2=r[r.Ra_b.notna()].copy(); r2['lRa']=np.log(r2.Ra_b)
m2=smf.mixedlm('lg~Cs+C(hc)+lRa',r2,groups=r2.cellkey).fit(); b,se=m2.params['Cs'],m2.bse['Cs']
print('  + Ra adjusted: ratio %.2f (%.2f-%.2f) p=%.3f; Ra coef %.2f'%(np.exp(b),np.exp(b-1.96*se),np.exp(b+1.96*se),m2.pvalues['Cs'],m2.params['lRa']))
r['lgraw']=np.log(r.p2p_raw); m3=smf.mixedlm('lgraw~Cs+C(hc)',r,groups=r.cellkey).fit(); b,se=m3.params['Cs'],m3.bse['Cs']
print('  raw 20 kHz trace: ratio %.2f (%.2f-%.2f) p=%.3f'%(np.exp(b),np.exp(b-1.96*se),np.exp(b+1.96*se),m3.pvalues['Cs']))
print(d.groupby('internal').cellkey.nunique(), d[d.h.isin(H)].groupby(['internal']).size())
print(cell.groupby('internal').Ra.median())
