import sys; sys.path.insert(0,'/home/claude/work/audit2')
import mlrun as m, pickle, numpy as np, pandas as pd, time, os
from scipy.signal import decimate
part,nparts=int(sys.argv[1]),int(sys.argv[2])
m.fl=pd.read_pickle('master_sheets.pkl')['files'].set_index('new_name')
t=pd.read_excel('../r3/audit/audit_labels.xlsx','audit_truth_events')
w=pd.read_excel('../r3/audit/audit_labels.xlsx','audit_windows'); w=w[w.batch.str.startswith('4')]
lab=pd.read_excel('../r3/audit/audit_labels.xlsx','audit_labels'); skip=set(lab[lab.note.astype(str).str.strip()=='skip'].window)
# windows with Cs/KGlu and >=4 truth events, up-heavy first
g=t.groupby('window').agg(n=('time_s','size'),nup=('polarity',lambda s:(s=='up').sum()),internal=('internal','first'))
g=g[(g.n>=4)&(g.internal!='KCl')].sort_values('nup',ascending=False)
ids=[i for i in g.index if i in set(w.window) and i not in skip][:40]
rows=[]; t0=time.time()
for k,wid in enumerate(ids):
    if k%nparts!=part: continue
    r=w[w.window==wid].iloc[0]; nm=r.file.replace('_Adv6_1B12_','_Adv3_1B12_')
    pre=pickle.load(open(f'preds/{wid}.pkl','rb'))
    tr=m.load(nm,r.start_s-3,r.end_s+3); tt=t[t.window==wid]
    for fac in (1,2,4):
        if fac==1: cur=None
        else:
            d=decimate(np.asarray(tr.data,dtype=np.float64),fac,zero_phase=True)
            from miniML import MiniTrace
            tr2=MiniTrace(data=d,sampling_interval=tr.sampling*fac,y_unit='pA'); tr2._off=tr._off
        for pol,dname in (('up','positive'),('down','negative')):
            if fac==1: ev=pre['res'][dname]['events']['t'].values
            else:
                e,_=m.detect(tr2,dname,thr=0.5); ev=e['t'].values
            ev=ev[(ev>=r.start_s)&(ev<=r.end_s)]
            T=tt[tt.polarity==pol].time_s.values
            hit=sum(1 for x in T if len(ev) and np.min(np.abs(ev-x))<=0.015)
            fp=sum(1 for e_ in ev if len(T)==0 or np.min(np.abs(T-e_))>0.015)
            rows.append(dict(window=wid,fac=fac,pol=pol,internal=tt.internal.iloc[0],n_truth=len(T),hit=hit,n_det=len(ev),fp=fp))
    print(k,wid,round(time.time()-t0),flush=True)
pd.DataFrame(rows).to_pickle(f'resamp_{part}.pkl')
