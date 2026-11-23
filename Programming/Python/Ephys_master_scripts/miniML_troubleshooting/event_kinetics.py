import pickle, numpy as np, pandas as pd
from scipy.signal import find_peaks
W=pickle.load(open('traces.pkl','rb'))
df=pd.read_pickle('truth_feat.pkl'); ps=pd.read_pickle('predscore.pkl')
df=df.merge(ps[['window','time_s','p_narrow']],on=['window','time_s'],how='left')
rows=[]
for r in df.itertuples():
    w=W[r.window]; y=w['y']; dt=w['meta']['dt']; t0=w['meta']['t0']
    s=-1 if r.pol=='down' else 1
    i=int((r.time_s-t0)/dt)
    a=max(0,i-80); b=min(len(y),i+120)
    seg=s*y[a:b]
    k=int(np.argmax(seg[60 if i>=80 else i-a-0:  (i-a)+40])) + (60 if i>=80 else 0) if False else None
    # refine peak within ±3 ms of marker
    lo=max(0,i-30); hi=min(len(y),i+30)
    j=lo+int(np.argmax(s*y[lo:hi])); pk=s*y[j]
    base=np.median(s*y[max(0,j-400):max(1,j-150)]) if j>150 else np.nan
    amp=pk-base
    rise=np.nan; dec=np.nan
    if amp>0:
        pre=s*y[max(0,j-300):j+1]-base
        i10=np.where(pre>=0.1*amp)[0]; i90=np.where(pre>=0.9*amp)[0]
        # last crossing before peak
        def last_cross(thr):
            idx=np.where(pre<thr)[0]
            return idx[-1] if len(idx) else 0
        rise=(len(pre)-1-last_cross(0.1*amp)-(len(pre)-1-last_cross(0.9*amp)))*dt*1000
        post=s*y[j:j+1500]-base
        h=np.where(post<=0.5*amp)[0]
        dec=h[0]*dt*1000 if len(h) else np.nan
    rows.append((rise,dec,amp))
df['rise_ms']=[x[0] for x in rows]; df['hdec_ms']=[x[1] for x in rows]; df['amp2']=[x[2] for x in rows]
df['blind']=(~df.found_M)&(df.p_narrow<0.1)
df.to_pickle('truth_feat2.pkl')
pd.set_option('display.width',200)
df['grp']=np.where(df.internal=='KCl','KCl',df.internal+' '+df.pol)
print(df.groupby(['found_M']).agg(rise=('rise_ms','median'),hdec=('hdec_ms','median'),amp=('amp2','median'),snr=('snr','median')).round(2))
print(df.groupby(['grp','found_M']).agg(n=('amp','size'),rise=('rise_ms','median'),hdec=('hdec_ms','median'),amp=('amp2','median'),snr=('snr','median')).round(2))
