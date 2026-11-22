import pandas as pd, numpy as np
from scipy import stats
pd.set_option('display.width',260,'display.max_rows',400,'display.max_columns',60)
a=pd.read_csv('k60_analysis_set2.csv')
a['abs_ihold']=a.Ihold_median_pA.abs()
a['noise']=a.noise_lp3k_pA
# comparator Cs -60
cs=pd.read_csv('k60_gapfree_features_cs60.csv'); cs['abs_drift']=cs.drift_total_pA.abs(); cs['abs_ihold']=cs.Ihold_median_pA.abs(); cs['noise']=cs.noise_lp3k_pA
cs['date']=cs.date.astype(int)
tcs=pd.read_csv('k60_transients_cs60.csv'); tcs=tcs[(tcs.amp_pA>=200)&(tcs.fwhm_ms<=3)]
print('Cs comparator',len(cs),cs.groupby(['session']).size().to_dict())
print(cs[['abs_ihold','abs_drift','noise','wander_sd_pA']].describe().round(1))
print(cs.groupby('session')[['noise','abs_drift']].median().round(1))
# thresholds
def pct(x,thr): return f"{(x>thr).sum()}/{len(x)} ({(x>thr).mean()*100:.0f}%)"
for gr,d in a.groupby('grp'):
    print(gr,len(d),'|drift|>10',pct(d.abs_drift,10),'>25',pct(d.abs_drift,25),'>50',pct(d.abs_drift,50),'| noise>4',pct(d.noise,4),'>5',pct(d.noise,5),'| |Ihold|>100',pct(d.abs_ihold,100),'>200',pct(d.abs_ihold,200),' transient files',(d.n_transients>0).sum(),' wander>10',pct(d.wander_sd_pA,10))
print('Cs60 |drift|>10',pct(cs.abs_drift,10),'>25',pct(cs.abs_drift,25),'>50',pct(cs.abs_drift,50),'noise>4',pct(cs.noise,4),'>5',pct(cs.noise,5),'|Ih|>100',pct(cs.abs_ihold,100))
# cell-level means
m60=a[a.grp=='-60 mV'].groupby('cellkey')[['abs_ihold','Ihold_median_pA','abs_drift','drift_total_pA','wander_sd_pA','noise','transient_rate_per_min','event_n','event_f','event_amp']].mean()
m0=a[a.grp=='0 mV'].groupby('cellkey')[['abs_ihold','Ihold_median_pA','abs_drift','drift_total_pA','wander_sd_pA','noise','transient_rate_per_min','event_n','event_f','event_amp']].mean()
both=m0.index.intersection(m60.index); print('paired cells',list(both))
print('paired (cell means, -60 vs 0)')
res=[]
for c in m0.columns:
    x=m60.loc[both,c]; y=m0.loc[both,c]
    ok=x.notna()&y.notna()
    if ok.sum()>=4:
        w=stats.wilcoxon(x[ok],y[ok]).pvalue
    else: w=np.nan
    res.append((c,ok.sum(),x[ok].mean(),y[ok].mean(),w))
print(pd.DataFrame(res,columns=['metric','n','mean_-60','mean_0','p_wilcoxon']).round(3))
# unpaired all -60 cells vs all 0 cells (cell level)
res=[]
for c in m0.columns:
    x=m60[c].dropna(); y=m0[c].dropna()
    res.append((c,len(x),len(y),x.median(),y.median(),stats.mannwhitneyu(x,y).pvalue))
print(pd.DataFrame(res,columns=['metric','n60','n0','med60','med0','p_MWU']).round(3))
# vs Cs comparator (recording level MWU)
for c in ['abs_ihold','abs_drift','noise','wander_sd_pA']:
    x=a[a.grp=='-60 mV'][c]; print(c,'KGlu-60 vs Cs-60 MWU p=',round(stats.mannwhitneyu(x,cs[c]).pvalue,4),'medians',round(x.median(),1),round(cs[c].median(),1))
# confirmed-only sensitivity for -60
conf=a[(a.grp=='-60 mV')&(a.label_status=='confirmed (journal/stated)')]
print('confirmed only n',len(conf),conf[['abs_ihold','abs_drift','noise','event_n','event_f']].median().round(2).to_dict())
print('all -60 n',(a.grp=='-60 mV').sum(),a[a.grp=='-60 mV'][['abs_ihold','abs_drift','noise','event_n','event_f']].median().round(2).to_dict())
# events
for gr,d in a.groupby('grp'):
    print(gr,'events total',d.event_n.sum(),'recs',len(d),'recs<5ev',(d.event_n<5).sum(),'recs0ev',(d.event_n==0).sum(),'rec>=20',(d.event_n>=20).sum(),'dur min',round(d.duration_s.sum()/60,1),' freq pooled',round(d.event_n.sum()/d.duration_s.sum(),3))
    print('   d2 inward total',d.d2_n_inward.sum(),'outward total',d.d2_n_outward.sum())
