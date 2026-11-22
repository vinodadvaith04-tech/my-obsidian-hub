import pandas as pd, numpy as np
A='/mnt/user-data/uploads/Ephys pilot & practice/audit/audit_labels.xlsx'
w=pd.read_excel(A,sheet_name='audit_labels'); wi=pd.read_excel(A,sheet_name='audit_windows')[['window','end_s','stratum','directions_shown']]
t=pd.read_excel(A,sheet_name='audit_truth_events'); ek=pd.read_excel(A,sheet_name='audit_events_key')
w=w.merge(wi,on='window'); w['dur']=w.end_s-w.start_s
print(w.dur.describe().round(2).to_dict())
# per-window per-polarity counts
rows=[]
for _,r in w.iterrows():
    tt=t[t.window==r.window]; ee=ek[ek.window==r.window]
    for pol in ('up','down'):
        T=tt[tt.polarity==pol]; M=ee[(ee.detector=='miniML')&(ee.polarity==pol)]
        O=ee[(ee.detector!='miniML')&(ee.polarity==pol)]
        rows.append(dict(window=r.window,pol=pol,internal=r.internal,drug=r.drug,holding=r.holding_mV,why=r.why_chosen,batch=r.batch,file=r.file,dur=r.dur,
            truth=len(T),mini_tp=int(T.found_by_miniML.sum()),mini_markers=len(M),alt_tp=int(T.found_by_old.sum()),alt_markers=len(O),clicked=int(T.clicked_by_you.sum())))
P=pd.DataFrame(rows); P.to_csv('audit_window_pol.csv',index=False)
print(P[(P.pol=='up')&(P.internal=='Cs-gluconate')].groupby(['holding','why'])[['truth','mini_tp','mini_markers','alt_tp','dur']].sum().assign(n=P[(P.pol=='up')&(P.internal=='Cs-gluconate')].groupby(['holding','why']).size()))
