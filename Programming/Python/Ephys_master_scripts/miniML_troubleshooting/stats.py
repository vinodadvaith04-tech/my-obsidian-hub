import sys; sys.path.insert(0,'.')
from load import *
import json
# false flags
fl=set()
for r in lab.itertuples():
    for col,pre in (('miniML_false_ids','M'),('old_false_ids','O')):
        v=getattr(r,col)
        if isinstance(v,str):
            for x in v.replace(',',';').split(';'):
                x=x.strip()
                if x: fl.add((r.window,pre+x.split('.')[0]))
key['false']=[(w,i) in fl for w,i in zip(key.window,key.id)]
key=key.merge(win[['window','internal','why_chosen','holding_mV','drug']],on='window')
key['grp']=np.where(key.internal=='KCl','KCl',key.internal+' '+key.polarity)
tru['grp']=np.where(tru.internal=='KCl','KCl',tru.internal+' '+tru.polarity)
def summarize(t,k):
    rows=[]
    for g in ['Cs-gluconate up','Cs-gluconate down','KGlu down','KGlu up','KCl']:
        a=t[t.grp==g]; b=k[k.grp==g]
        n=len(a)
        if n==0: continue
        m=b[b.detector=='miniML']; o=b[b.detector!='miniML']
        rows.append(dict(group=g,windows=a.window.nunique(),truth=n,
            rec_M=a.found_by_miniML.mean(),rec_D2=a.found_by_old.mean(),rec_either=(a.found_by_miniML|a.found_by_old).mean(),
            rec_clickonly=(a.clicked_by_you&~a.found_by_miniML&~a.found_by_old).mean(),
            prec_M=1-m.false.mean(),prec_D2=1-o.false.mean(),nM=len(m),nD2=len(o)))
    return pd.DataFrame(rows)
allw=summarize(tru,key)
rnd=summarize(tru[tru.why_chosen=='random'],key[key.why_chosen=='random'])
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
print('ALL windows'); print(allw.round(3).to_string())
print('RANDOM windows'); print(rnd.round(3).to_string())
allw.to_csv('summary_all.csv',index=False); rnd.to_csv('summary_random.csv',index=False)
key.to_pickle('key.pkl'); tru.to_pickle('tru.pkl')
