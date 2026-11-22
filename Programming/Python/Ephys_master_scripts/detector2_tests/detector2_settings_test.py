import sys, copy, numpy as np, pandas as pd, pyabf, pickle, time
sys.path.insert(0,'/home/claude/work/r3/scripts'); import psc_detection as P
from multiprocessing import Pool
t=pd.read_excel('/home/claude/work/r3/audit_latest/audit_labels.xlsx','audit_truth_events')
W=pd.read_pickle('/home/claude/work/r3/exp_windows.pkl'); W=W[W.ok]
W=W[W.window.isin(set(t.window))]
rng=np.random.default_rng(3); sel=[]
for internal,n in (('Cs-gluconate',64),('KGlu',28),('KCl',28)):
    g=W[W.internal==internal]; sel+=list(g.window.values[rng.permutation(len(g))[:n]])
W=W.set_index('window').loc[sorted(sel)]
def cfg(**kw):
    c=copy.deepcopy(P.CFG)
    for k,v in kw.items(): c[k]=v
    return c
base=P.CFG
T2=dict(inward=list(base['TEMPLATES']['inward'])+[(1.0,5.0,35.0)], outward=list(base['TEMPLATES']['outward'])+[(2.5,12.0,80.0),(6.0,40.0,160.0)])
C={'base':cfg(),
   'K3':cfg(K_SD=dict(inward=3.0,outward=3.5),ABS_MIN_PA=3.0),
   'K2.5':cfg(K_SD=dict(inward=2.5,outward=3.0),ABS_MIN_PA=3.0),
   'UF4':cfg(UNSTABLE_FACTOR=4.0),
   'crit2.5':cfg(CAND_CRIT=2.5),
   'tmpl+':cfg(TEMPLATES=T2),
   'combo':cfg(K_SD=dict(inward=3.0,outward=3.5),ABS_MIN_PA=3.0,UNSTABLE_FACTOR=4.0,CAND_CRIT=2.5,TEMPLATES=T2)}
def run(wid):
    r=W.loc[wid]; ab=pyabf.ABF(r.path); ab.setSweep(0); I=np.asarray(ab.sweepY,float); fs=float(ab.dataRate)
    s0=max(0,r.start_s-5); i0=int(s0*fs); i1=min(len(I),int((r.end_s+5)*fs)); seg=I[i0:i1]
    out=[]
    for name,c in C.items():
        try: ok,rej,info=P.detect(seg,fs,cfg=c)
        except Exception as e: print('ERR',wid,name,e,flush=True); continue
        for pol in ('inward','outward'):
            d='down' if pol=='inward' else 'up'
            if ('negative' if pol=='inward' else 'positive') not in r.directions_shown: continue
            T=t[(t.window==wid)&(t.polarity==d)].time_s.values
            ev=np.array([x['peak_s']+s0 for x in ok if x['polarity']==pol]); ev=ev[(ev>=r.start_s)&(ev<=r.end_s)] if len(ev) else ev
            hit=sum(1 for x in T if len(ev) and np.min(np.abs(ev-x))<=0.015)
            fp=sum(1 for e in ev if len(T)==0 or np.min(np.abs(T-e))>0.015)
            out.append(dict(window=wid,cfg=name,pol=d,internal=r.internal,truth=len(T),hit=hit,det=len(ev),fp=fp))
    return out
if __name__=='__main__':
    t0=time.time()
    with Pool(2) as p:
        res=[]
        for i,x in enumerate(p.imap(run,list(W.index),chunksize=2)):
            res.append(x)
            if i%10==0: print(i,round(time.time()-t0),flush=True)
    pd.DataFrame([y for x in res for y in x]).to_pickle('exp_d2.pkl'); print('DONE',round(time.time()-t0))
