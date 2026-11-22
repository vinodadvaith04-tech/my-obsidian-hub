import sys, numpy as np, pandas as pd, pyabf, pickle, time
sys.path.insert(0,'/home/claude/work/r3/scripts'); import psc_detection as P
from multiprocessing import Pool
t=pd.read_excel('/home/claude/work/r3/audit_latest/audit_labels.xlsx','audit_truth_events')
W=pd.read_pickle('/home/claude/work/r3/exp_windows.pkl'); W=W[W.ok].set_index('window')
miss=t[(~t.found_by_old)&t.window.isin(W.index)]
wins=sorted(miss.window.unique())
CFG=P.CFG; TOL=0.015
def run(wid):
    r=W.loc[wid]; ab=pyabf.ABF(r.path); ab.setSweep(0); I=np.asarray(ab.sweepY,float); fs=float(ab.dataRate)
    s0=max(0,r.start_s-5); i0=int(s0*fs); i1=min(len(I),int((r.end_s+5)*fs)); seg=I[i0:i1]
    ok,rej,info,aux=P.detect(seg,fs,return_aux=True)
    out=[]
    for e in miss[miss.window==wid].itertuples():
        pol='inward' if e.polarity=='down' else 'outward'; sgn=-1.0 if pol=='inward' else 1.0
        tt=e.time_s-s0
        okn=[x for x in ok if x['polarity']==pol and abs(x['peak_s']-tt)<=TOL]
        lim=info['limit_inward_pA' if pol=='inward' else 'limit_outward_pA']
        d=dict(window=wid,time_s=e.time_s,pol=e.polarity,internal=e.internal,hold=e.holding_mV,amp=e.amp_pA,limit=lim,noise=info['noise_lp_pA'],
               found_miniML=e.found_by_miniML,clicked=e.clicked_by_you)
        if okn: d['cat']='found in re-run'; out.append(d); continue
        rj=[x for x in rej if x['polarity']==pol and x.get('peak_s') is not None and abs(x['peak_s']-tt)<=TOL]
        if rj:
            x=min(rj,key=lambda x:abs(x['peak_s']-tt)); d['cat']='rejected: '+x.get('reason','?'); d['rise']=x.get('rise1090_ms'); d['td']=x.get('td'); d['resid']=x.get('resid_ratio'); d['fitA']=x.get('A'); out.append(d); continue
        cands=P.candidates(aux['y'],aux['ylp'],fs,CFG,pol,valid=~aux['mask'])
        near=[c for c in cands if tt-0.040<=c['i0']/fs<=tt+0.003]
        if not near:
            if aux['mask'][min(len(aux['mask'])-1,int(tt*fs))]: d['cat']='in masked/unstable stretch'
            else: d['cat']='no stage-1 candidate (match score < 3)'
        else:
            c=max(near,key=lambda c:c['crit']); rp=P.rough_peak(sgn*aux['ylp'],fs,c['i0'],pol,CFG)
            if rp is None: d['cat']='dropped: rough peak on a slope'
            elif rp[1]<CFG['PREGATE']*lim: d['cat']='dropped: below 0.6 x amplitude limit'
            else: d['cat']='dropped: merged with neighbour / duplicate'
            d['crit']=c['crit']
        out.append(d)
    return out
if __name__=='__main__':
    t0=time.time()
    with Pool(2) as p: res=p.map(run,wins,chunksize=4)
    df=pd.DataFrame([x for r in res for x in r]); df.to_pickle('diag.pkl'); print(len(df),'missed events in',len(wins),'windows',round(time.time()-t0),'s')
