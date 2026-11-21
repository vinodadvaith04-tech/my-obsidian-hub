import pickle, numpy as np, pandas as pd, glob
t=pd.read_excel('/home/claude/work/r3/audit_latest/audit_labels.xlsx','audit_truth_events')
out=[]
for f in glob.glob('exp_[01].pkl'): out+=pickle.load(open(f,'rb'))
TOL=0.015; rows=[]
for o in out:
    pol='down' if o['dir']=='negative' else 'up'
    T=t[(t.window==o['window'])&(t.polarity==pol)].time_s.values
    ev=o['t']; ev=ev[(ev>=o['start'])&(ev<=o['end'])] if len(ev) else ev
    hit=sum(1 for x in T if len(ev) and np.min(np.abs(ev-x))<=TOL)
    fp=sum(1 for e in ev if len(T)==0 or np.min(np.abs(T-e))>TOL)
    rows.append(dict(window=o['window'],dir=o['dir'],pol=pol,cfg=o['cfg'],internal=o['internal'],hold=o['hold'],truth=len(T),hit=hit,det=len(ev),fp=fp))
df=pd.DataFrame(rows); df.to_pickle('exp_eval.pkl')
def agg(g):
    tr=g.truth.sum(); h=g.hit.sum(); d=g.det.sum(); fp=g.fp.sum()
    p=(d-fp)/d if d else np.nan; r=h/tr if tr else np.nan
    return pd.Series(dict(truth=tr,det=d,recall=round(r,3),precision=round(p,3),F1=round(2*p*r/(p+r),3) if d and tr else np.nan,fp=fp))
if __name__=='__main__':
    print(df.window.nunique(),'windows', len(df)//df.cfg.nunique(),'window-directions')
    order=['base','w150','w300','w1200','w2400','thr.7','thr.3','thr.2','thr.1','rel.1','rel.5','pw2','pw10','lp2k','invert','w300thr.3','w1200thr.3','w300thr.2','w1200thr.2','sens']
    for name,sub in (('ALL',df),('Cs-gluconate',df[df.internal=='Cs-gluconate']),('KGlu',df[df.internal=='KGlu']),('KCl',df[df.internal=='KCl']),
                     ('Cs down (EPSC-like)',df[(df.internal=='Cs-gluconate')&(df.pol=='down')]),('Cs up (IPSC-like)',df[(df.internal=='Cs-gluconate')&(df.pol=='up')])):
        r=sub.groupby('cfg').apply(agg).reindex(order); print('\n==',name); print(r.to_string())
