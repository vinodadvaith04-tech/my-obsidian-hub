import sys; sys.path.insert(0,'.')
import pandas as pd, numpy as np
pd.set_option('display.width',250)
t=pd.read_pickle('tru.pkl'); k=pd.read_pickle('key.pkl')
def f1(p,r): return 2*p*r/(p+r)
# overall + F1 (all windows) per group
a=pd.read_csv('summary_all.csv'); r=pd.read_csv('summary_random.csv')
for d in (a,r):
    d['F1_M']=f1(d.prec_M,d.rec_M); d['F1_D2']=f1(d.prec_D2,d.rec_D2)
print(a.round(3)[['group','truth','rec_M','prec_M','F1_M','rec_D2','prec_D2','F1_D2','rec_either']]); print(r.round(3)[['group','truth','rec_M','prec_M','F1_M','rec_D2','prec_D2','F1_D2','rec_either']])
# overall all groups
for nm,tt,kk in (('all',t,k),('random',t[t.why_chosen=='random'],k[k.why_chosen=='random'])):
    rm=tt.found_by_miniML.mean(); rd=tt.found_by_old.mean()
    pm=1-kk[kk.detector=='miniML'].false.mean(); pd_=1-kk[kk.detector!='miniML'].false.mean()
    print(nm,'overall M rec %.3f prec %.3f F1 %.3f | D2 rec %.3f prec %.3f F1 %.3f | either %.3f'%(rm,pm,f1(pm,rm),rd,pd_,f1(pd_,rd),(tt.found_by_miniML|tt.found_by_old).mean()))
# overlap categories
t['cat']=np.select([t.found_by_miniML&t.found_by_old,t.found_by_miniML&~t.found_by_old,~t.found_by_miniML&t.found_by_old],['both','miniML only','D2 only'],'neither (clicked)')
print(pd.crosstab(t.grp,t.cat,normalize='index').round(3)); print(pd.crosstab(t.grp,t.cat))
# holding classes for Cs up
t['hc']=np.where(t.holding_mV>=0,'>=0 mV',np.where(t.holding_mV<=-55,'<=-55 mV','-50..-30 mV'))
print(t.groupby(['grp','hc']).agg(n=('hc','size'),M=('found_by_miniML','mean'),D2=('found_by_old','mean')).round(2))
# drug
print(t.groupby(['grp','drug']).agg(n=('hc','size'),M=('found_by_miniML','mean'),D2=('found_by_old','mean')).round(2))
# network probability categories for misses
ps=pd.read_pickle('predscore2.pkl') if False else pd.read_pickle('/home/claude/work/audit2/predscore2.pkl')
ps['grp']=np.where(ps.internal=='KCl','KCl',ps.internal+' '+ps.pol)
m=ps[~ps.found_M].copy()
m['band']=pd.cut(m.pmax,[-1,0.1,0.3,0.5,1.01],labels=['<0.1','0.1-0.3','0.3-0.5','>=0.5'])
x=pd.crosstab(m.grp,m.band); print(x); print(pd.crosstab(m.grp,m.band,normalize='index').round(2))
print('overall',m.band.value_counts(normalize=True).round(3).to_dict(), len(m), 'events with preds',len(ps))
x.to_csv('pmax_bands.csv')
