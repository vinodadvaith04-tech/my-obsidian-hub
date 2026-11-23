import pandas as pd, numpy as np
P='/mnt/user-data/uploads/Ephys pilot & practice/audit/audit_labels.xlsx'
X=pd.ExcelFile(P)
lab=X.parse('audit_labels'); win=X.parse('audit_windows'); key=X.parse('audit_events_key'); tru=X.parse('audit_truth_events')
win['hold_class']=np.where(win.holding_mV>=0,'positive (>=0 mV)',np.where(win.holding_mV<=-55,'negative (<=-55 mV)','intermediate'))
tru=tru.merge(win[['window','why_chosen','stratum']],on='window',how='left')
