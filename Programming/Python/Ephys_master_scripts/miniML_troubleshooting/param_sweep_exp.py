import sys, os, time, pickle, warnings
os.environ['TF_USE_LEGACY_KERAS']='1'; os.environ['TF_CPP_MIN_LOG_LEVEL']='3'; warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/claude/work/miniML/core')
import numpy as np, pandas as pd
if not hasattr(np,'trapz'): np.trapz=np.trapezoid
import tensorflow as tf
from scipy.signal import butter, sosfiltfilt
from miniML import MiniTrace, EventDetection
part,nparts=int(sys.argv[1]),int(sys.argv[2])
MODEL=tf.keras.models.load_model('/home/claude/work/miniML/models/GC_lstm_model.h5',compile=False)
W=pd.read_pickle('/home/claude/work/r3/exp_windows.pkl'); W=W[W.ok]
t=pd.read_excel('/home/claude/work/r3/audit_latest/audit_labels.xlsx','audit_truth_events')
ntr=t.groupby('window').size(); W=W[W.window.map(ntr).fillna(0)>=2]
rng=np.random.default_rng(7)
sel=[]
for internal,n in (('Cs-gluconate',34),('KGlu',18),('KCl',8)):
    g=W[W.internal==internal]; sel+=list(g.window.values[rng.permutation(len(g))[:n]])
sel=sorted(sel); W=W.set_index('window').loc[sel]
if part==0: pd.Series(sel).to_csv('sel_windows.csv',index=False)
# config: name -> (win, thr, peak_w, rel, prefilter, invert)
C={'base':(600,.5,5,.25,0,0),
 'w150':(150,.5,5,.25,0,0),'w300':(300,.5,5,.25,0,0),'w1200':(1200,.5,5,.25,0,0),'w2400':(2400,.5,5,.25,0,0),
 'thr.7':(600,.7,5,.25,0,0),'thr.3':(600,.3,5,.25,0,0),'thr.2':(600,.2,5,.25,0,0),'thr.1':(600,.1,5,.25,0,0),
 'rel.1':(600,.5,5,.10,0,0),'rel.5':(600,.5,5,.50,0,0),'pw2':(600,.5,2,.25,0,0),'pw10':(600,.5,10,.25,0,0),
 'lp2k':(600,.5,5,.25,1,0),'invert':(600,.5,5,.25,0,1),
 'w300thr.3':(300,.3,5,.25,0,0),'w1200thr.3':(1200,.3,5,.25,0,0),'w300thr.2':(300,.2,5,.25,0,0),'w1200thr.2':(1200,.2,5,.25,0,0),
 'sens':(600,.2,2,.10,0,0)}
CACHE={}
orig=EventDetection._EventDetection__predict
def patched(self):
    key=(self._k,self.window_size,self.event_direction)
    if key in CACHE: self.prediction=CACHE[key]
    else: orig(self); CACHE[key]=self.prediction
EventDetection._EventDetection__predict=patched
sos=butter(4,2000,'low',fs=10000,output='sos')
out=[]; t0=time.time()
for k,(wid,r) in enumerate(W.iterrows()):
    if k%nparts!=part: continue
    tr0=MiniTrace.from_axon_file(filename=r.path,channel=0,scaling=1.0,unit='pA')
    fs=1/tr0.sampling; i0=max(0,int((r.start_s-3)*fs)); i1=min(len(tr0.data),int((r.end_s+3)*fs))
    raw=np.asarray(tr0.data,dtype=np.float64)[i0:i1].copy(); off=i0*tr0.sampling
    CACHE.clear()
    for d in r.directions_shown.split(','):
        for name,(win,thr,pw,rel,pf,inv) in C.items():
            data=raw
            if pf: data=sosfiltfilt(sos,data)
            dd=d
            if inv: data=-data; dd='positive' if d=='negative' else 'negative'
            tr=MiniTrace(data=data.copy(),sampling_interval=tr0.sampling,y_unit='pA')
            det=EventDetection(data=tr,model=MODEL,model_threshold=thr,window_size=win,batch_size=512,event_direction=dd,verbose=0)
            det._k=(wid,pf,inv)
            try:
                det.detect_events(eval=True,peak_w=pw,rel_prom_cutoff=rel,convolve_win=20,gradient_convolve_win=40)
                if hasattr(det,'event_peak_times'):
                    tt=np.asarray(det.event_peak_times)+off; sc=np.asarray(det.event_scores)
                else: tt=np.array([]); sc=np.array([])
            except Exception as e:
                tt=np.array([]); sc=np.array([]); print('ERR',wid,name,repr(e)[:100],flush=True)
            out.append(dict(window=wid,dir=d,cfg=name,t=tt,score=sc,start=r.start_s,end=r.end_s,internal=r.internal,hold=r.holding_mV))
    print(k,wid,round(time.time()-t0),flush=True)
    pickle.dump(out,open(f'exp_{part}.pkl','wb'))
pickle.dump(out,open(f'exp_{part}.pkl','wb')); print('DONE',flush=True)
