import json, numpy as np, pandas as pd
raw=json.load(open('paper_raw.json'))
cal={'DIV9_sEPSC':[(1.5,172.0),(1.0,193.6),(0.5,215.3),(0,237.0)],
 'DIV9_sIPSC':[(8,275.5),(4,308.0),(0,340.5)],
 'DIV21_sEPSC':[(2.0,223.1),(1.5,244.6),(1.0,266.0),(0,287.5)],
 'DIV21_sIPSC':[(20,334.3),(10,366.5),(0,398.7)]}
rows=[]
for k,v in raw.items():
    c=np.array([[a,b] for a,b,_ in v['circ']])
    xs=np.sort(c[:,0]); 
    gaps=np.diff(xs); i=gaps.argmax(); split=(xs[i]+xs[i+1])/2
    ys=np.array([p[1] for p in cal[k]]); vals=np.array([p[0] for p in cal[k]])
    m,b=np.polyfit(ys,vals,1)
    print(k,'split',round(split,1),'gap',round(gaps[i],1),'resid',np.abs(m*ys+b-vals).max().round(3))
    for x,y in c:
        rows.append(dict(panel=k,group='CTL' if x<split else 'VU',x=x,y=y,value=m*y+b))
df=pd.DataFrame(rows)
print(df.groupby(['panel','group']).value.agg(['count','mean','median','std','min','max']).round(2))
df.to_csv('paper_dots.csv',index=False)
