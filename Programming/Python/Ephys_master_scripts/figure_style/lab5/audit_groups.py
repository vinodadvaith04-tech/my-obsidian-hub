import pandas as pd, numpy as np
P=pd.read_csv('/home/claude/work/lab5/audit_window_pol.csv')
def grp(r):
    if r.internal=='Cs-gluconate':
        if r.pol=='up' and r.holding in (0,10): return 'I: Cs, 0 and +10 mV'
        if r.pol=='down' and r.holding in (-80,-70): return 'E: Cs, -80 and -70 mV'
    if r.internal=='KCl' and r.drug=='DNQX+AP5' and r.pol=='down': return 'I: KCl + DNQX/AP5, -60 mV'
    if r.internal=='KCl' and r.drug=='blockers not confirmed' and r.pol=='down': return 'KCl, blockers not confirmed (not used)'
    if r.internal=='KGlu' and r.drug=='gabazine' and r.pol=='down': return 'E: K-Glu + gabazine, -60 mV'
    return None
P['group']=P.apply(grp,axis=1); G=P.dropna(subset=['group']).copy()
for c in ('truth','mini_tp','mini_markers','alt_tp','alt_markers'): G[c+'_hz']=G[c]/G.dur
G.to_csv('/home/claude/work/lab5/audit_groups_windows.csv',index=False)
out=[]
for sel,name in (('random','random windows only'),('all','all audited windows')):
    g=G if sel=='all' else G[G.why=='random']
    for k,x in g.groupby('group'):
        out.append(dict(sel=name,group=k,windows=len(x),files=x.file.nunique(),truth_hz_median=x.truth_hz.median(),truth_hz_q1=x.truth_hz.quantile(.25),truth_hz_q3=x.truth_hz.quantile(.75),truth_hz_mean=x.truth_hz.mean(),
          mini_all_hz_mean=x.mini_markers_hz.mean(),mini_tp_hz_mean=x.mini_tp_hz.mean(),alt_tp_hz_mean=x.alt_tp_hz.mean(),
          recall_mini=x.mini_tp.sum()/max(x.truth.sum(),1),precision_mini=x.mini_tp.sum()/max(x.mini_markers.sum(),1),recall_alt=x.alt_tp.sum()/max(x.truth.sum(),1),n_truth=int(x.truth.sum())))
R=pd.DataFrame(out); R.to_csv('/home/claude/work/lab5/audit_group_summary.csv',index=False)
pd.set_option('display.width',250); pd.set_option('display.max_columns',30)
print(R.round(2).to_string())
