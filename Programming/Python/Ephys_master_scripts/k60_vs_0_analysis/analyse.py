import pandas as pd, numpy as np, re, json
from scipy import stats
pd.set_option('display.width',260,'display.max_rows',400,'display.max_columns',60)
ff=pd.read_csv('k60_files_slice.csv'); ft=pd.read_csv('k60_gapfree_features.csv'); cc=pd.read_csv('k60_capcheck.csv')
mm=pd.read_csv('k60_miniML_slice.csv'); al=pd.read_csv('k60_alt_slice.csv'); ev=pd.read_csv('k60_miniML_events_slice.csv'); tr=pd.read_csv('k60_transients.csv')
g=ft.merge(ff[['new_name','holding_flag','sex','T.Code']],on='new_name')
g['t']=pd.to_datetime(g.start_time)
g['basis']=np.select([g.holding_flag.fillna('').str.startswith('default'),g.holding_flag.fillna('').str.startswith('stated')],['assumed','stated'],'journal')
# ---- cap checks: clamp voltage, nearest before/after
cc['t']=pd.to_datetime(cc.time)
cc['V']=np.where(cc.v_pre.abs()>200,cc.v_pre/20,np.where(cc.ch1.str.endswith('/pA'),np.nan,cc.v_pre))
for col in ['Ra_MOhm','Rm_MOhm','Cm_pF','Ibase_pA','V']:
    pass
def near(r,side):
    s=cc[(cc.date==r.date)&(cc.slice==r.slice)&(cc.cell==r.cell)]
    s=s[s.t<r.t] if side=='b' else s[s.t>r.t]
    if side=='b': s=s.sort_values('t').tail(2)
    else: s=s.sort_values('t').head(2)
    return s
rows=[]
for r in g.itertuples():
    d=dict(new_name=r.new_name)
    for side in ('b','a'):
        s=near(r,side)
        s_ok=s[s.Ra_MOhm.notna()]
        d[f'V_{side}']=s.V.median() if len(s) else np.nan
        d[f'gap_{side}_min']=((r.t-s.t.iloc[-1]).total_seconds()/60 if side=='b' else (s.t.iloc[0]-r.t).total_seconds()/60) if len(s) else np.nan
        for c,n in [('Ra_MOhm','Ra'),('Rm_MOhm','Rm'),('Cm_pF','Cm'),('Ibase_pA','Ib')]:
            d[f'{n}_{side}']=s_ok[c].median() if len(s_ok) else np.nan
    rows.append(d)
g=g.merge(pd.DataFrame(rows),on='new_name')
# ---- label status
def status(r):
    tgt=r.holding_command_mV
    vs=[v for v in (r.V_b,r.V_a) if pd.notna(v)]
    if r.basis!='assumed': return 'confirmed (journal/stated)'
    if tgt==-60:
        if vs and all(abs(v-(-60))>5 for v in vs): return 'contradicted by cap check'
        if any(abs(v+60)<=5 for v in vs): return 'assumed, cap checks at -60'
    return 'assumed, unverified'
g['label_status']=g.apply(status,axis=1)
g['grp']=np.select([g.holding_command_mV==-60,g.holding_command_mV==0],['-60 mV','0 mV'],'+5 to +15 mV')
# ---- transients episodes
t=tr[(tr.amp_pA>=200)&(tr.fwhm_ms<=3)].sort_values(['new_name','t_s'])
ep=[]
for nm,gg in t.groupby('new_name'):
    last=-9
    for r in gg.itertuples():
        if r.t_s-last>0.1: ep.append(dict(new_name=nm,t_s=r.t_s,sign=r.sign,amp=r.amp_pA,fwhm=r.fwhm_ms))
        elif r.amp_pA>ep[-1]['amp']: ep[-1].update(amp=r.amp_pA,fwhm=r.fwhm_ms,sign=r.sign)
        last=r.t_s
ep=pd.DataFrame(ep)
g['n_transients']=g.new_name.map(ep.groupby('new_name').size()).fillna(0).astype(int)
g['transient_rate_per_min']=g.n_transients/(g.duration_s/60)
g['max_transient_pA']=g.new_name.map(ep.groupby('new_name').amp.max())
# ---- miniML / detector 2
mm['dirn']=mm.direction
neg=mm[mm.direction=='negative'].set_index('new_name'); pos=mm[mm.direction=='positive'].set_index('new_name')
for nm,d in (('neg',neg),('pos',pos)):
    g[f'n_{nm}']=g.new_name.map(d.n_events); g[f'f_{nm}']=g.new_name.map(d.frequency_Hz)
    g[f'amp_{nm}']=g.new_name.map(d.amplitude_abs_median); g[f'valid_{nm}']=g.new_name.map(d.valid_time_s)
    g[f'cut_{nm}']=g.new_name.map(d.cut_time_s); g[f'noise_ml_{nm}']=g.new_name.map(d.noise_mad_pA)
al=al.set_index('new_name')
for c in ['n_inward','n_outward','freq_inward_Hz','freq_outward_Hz','noise_lp_pA','analysed_s']:
    g['d2_'+c]=g.new_name.map(al[c])
g.to_csv('k60_master_table.csv',index=False); ep.to_csv('k60_episodes.csv',index=False)
print(g.groupby(['grp','label_status']).size())
print(g.groupby('grp')[['cut_neg','cut_pos','valid_neg','valid_pos']].describe().T.round(1).head(0))
print(g[['n_neg','n_pos']].notna().groupby(g.grp).sum())
