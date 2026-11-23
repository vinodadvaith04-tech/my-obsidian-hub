import pickle, glob, numpy as np, pandas as pd
t=pd.read_excel('../r3/audit/audit_labels.xlsx','audit_truth_events')
rows=[]
for fn in sorted(glob.glob('preds/W*.pkl')):
    d=pickle.load(open(fn,'rb')); wid=d['meta']['id']; g=t[t.window==wid]
    for r in g.itertuples():
        d_='negative' if r.polarity=='down' else 'positive'
        if d_ not in d['res']: continue
        R=d['res'][d_]; p=R['pred']; dt=R['dt']; off=R['off']
        i=int((r.time_s-off)/dt)
        a,b=max(0,i-int(0.05/dt)),min(len(p),i+int(0.05/dt))
        if b<=a: continue
        seg=p[a:b]
        rows.append(dict(window=wid,time_s=r.time_s,pol=r.polarity,found_M=r.found_by_miniML,found_O=r.found_by_old,clicked=r.clicked_by_you,pmax=float(seg.max()),internal=r.internal,hold=r.holding_mV,amp=abs(r.amp_pA),batch=r.batch))
df=pd.DataFrame(rows); df.to_pickle('predscore2.pkl'); print(len(df))
df['cls']=pd.cut(df.pmax,[-1,0.1,0.3,0.5,0.8,1.01],labels=['<0.1','0.1-0.3','0.3-0.5','0.5-0.8','>=0.8'])
print(pd.crosstab([df.pol,df.found_M],df.cls))
miss=df[~df.found_M]
print('\nMISSED events:',len(miss))
for pol,g in miss.groupby('pol'):
    print(pol,len(g),'pmax<0.5: %.0f%%  pmax>=0.5: %.0f%%'%(100*(g.pmax<0.5).mean(),100*(g.pmax>=0.5).mean()),' pmax<0.1: %.0f%%'%(100*(g.pmax<0.1).mean()))
print('\nmissed by internal/pol'); print(miss.groupby(['internal','pol']).apply(lambda g: pd.Series(dict(n=len(g),blind=(g.pmax<0.5).mean().round(2),median_pmax=g.pmax.median().round(2)))))
print('\nround 2 only'); m2=miss[miss.batch.str.startswith('4')]; print(m2.groupby('pol').apply(lambda g: pd.Series(dict(n=len(g),blind=(g.pmax<0.5).mean().round(2)))))
