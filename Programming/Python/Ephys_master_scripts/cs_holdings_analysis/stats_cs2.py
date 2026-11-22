import pandas as pd, numpy as np, json
from scipy import stats
pd.set_option('display.width',250,'display.max_columns',80,'display.max_rows',400)
g=pd.read_csv('cs_master_table.csv')
a=g[(g.use_for_miniML=='Y')&g.holding_command_mV.notna()&(g.duration_s>=30)].copy()
a['ra_ok']=a.Ra_b<35
a['leaky']=a.Ihold_median_pA.abs()>1000
strict=a[a.ra_ok&~a.leaky].copy()
print('primary',len(a),a.cellkey.nunique(),'strict',len(strict),strict.cellkey.nunique())
print('excluded in strict:',a[~(a.ra_ok&~a.leaky)][['new_name','holding_command_mV','Ra_b','Ihold_median_pA']].to_string())
H=[-80,-70,-60,-50,-40,-30,-20,0,10,15]
def med(v): 
    v=np.asarray(v,float); v=v[~np.isnan(v)]
    return (np.nan,)*3 if len(v)==0 else tuple(np.percentile(v,[50,25,75]))
def cellmean(d,col): return d.groupby(['cellkey','holding_command_mV'])[col].mean().unstack()
out={}
# ---- Table 1
t1=[]
for h in H:
    d=a[a.holding_command_mV==h]
    t1.append(dict(h=h,days=d.date.nunique(),slices=d[['date','slice']].drop_duplicates().shape[0],cells=d.cellkey.nunique(),recs=len(d),minutes=round(d.valid_in.sum()/60 if 'valid_in' in d else d.valid_in.sum(),1) if False else round(d.duration_s.sum()/60,1),
       cells_per=round(len(d)/d.cellkey.nunique(),1)))
out['t1']=t1
print(pd.DataFrame(t1))
# ---- Table 2 seal test
t2=[]
for h in H:
    d=a[a.holding_command_mV==h]; r=dict(h=h,cells=d.cellkey.nunique())
    for c,n in (('Ra_b','Ra'),('Rm_b','Rm'),('Cm_b','Cm')):
        m=d.groupby('cellkey')[c].mean(); r[n]=med(m); r[n+'_mean']=(m.mean(),m.std()/np.sqrt(m.notna().sum()))
    r['ra35_rec']=int((d.Ra_b>=35).sum()); r['n_rec']=len(d); r['ra_max']=d.Ra_b.max()
    t2.append(r)
out['t2']=t2
for r in t2: print(r['h'],r['cells'],[round(x,1) for x in r['Ra']],[round(x) for x in r['Rm']],[round(x) for x in r['Cm']],r['ra35_rec'],r['n_rec'],round(r['ra_max'],1))
# ---- Table 3 noise drift ihold
t3=[]
for h in H:
    d=a[a.holding_command_mV==h]; r=dict(h=h,cells=d.cellkey.nunique())
    for c,n in (('noise_lp3k_pA','noise'),('drift_total_pA','drift'),('Ihold_median_pA','ihold'),('wander_sd_pA','wander')):
        m=d.groupby('cellkey')[c].mean(); r[n]=med(m); r[n+'_mean']=(m.mean(),m.std()/np.sqrt(m.notna().sum()))
    r['transients_per_rec']=med(d.groupby('cellkey').n_transients.mean())
    r['trans_rec_with']=int((d.n_transients>0).sum()); r['n_rec']=len(d)
    t3.append(r)
out['t3']=t3
for r in t3: print(r['h'],r['cells'],[round(x,1) for x in r['noise']],[round(x,1) for x in r['drift']],[round(x) for x in r['ihold']],[round(x,1) for x in r['transients_per_rec']])
# ---- Table 4 events
t4=[]
for h in H:
    d=a[a.holding_command_mV==h]; r=dict(h=h,cells=d.cellkey.nunique())
    for c in ('f_in','f_out','amp_in','amp_out','q_in','q_out','d2_freq_inward_Hz','d2_freq_outward_Hz'):
        m=d.groupby('cellkey')[c].mean(); r[c]=med(m); r[c+'_mean']=(m.mean(),m.std()/np.sqrt(m.notna().sum()) if m.notna().sum()>1 else np.nan); r[c+'_n']=int(m.notna().sum())
    r['n_ev_in']=int(d.n_in.sum()); r['n_ev_out']=int(d.n_out.sum())
    t4.append(r)
out['t4']=t4
for r in t4: print(r['h'],r['cells'],'fin',[round(x,2) for x in r['f_in']],'fout',[round(x,2) for x in r['f_out']],'ain',round(r['amp_in'][0],1),'aout',round(r['amp_out'][0],1),'d2',round(r['d2_freq_inward_Hz'][0],2),round(r['d2_freq_outward_Hz'][0],2),r['n_ev_in'],r['n_ev_out'])
json.dump(out,open('cs_tables_part1.json','w'),default=lambda o:float(o) if isinstance(o,(np.floating,np.integer)) else str(o))
a.to_csv('cs_analysis_set.csv',index=False); strict.to_csv('cs_analysis_set_strict.csv',index=False)
