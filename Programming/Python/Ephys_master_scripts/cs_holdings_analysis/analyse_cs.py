import pandas as pd, numpy as np
pd.set_option('display.width',250,'display.max_columns',80,'display.max_rows',300)
f=pd.read_pickle('/tmp/files_new.pkl'); mm=pd.read_pickle('/tmp/miniML_per_file_new.pkl'); al=pd.read_pickle('/tmp/alt_det_per_file_new.pkl')
ft=pd.read_csv('k60_gapfree_features_cs.csv'); cc=pd.read_csv('k60_capcheck_cs.csv'); tr=pd.read_csv('k60_transients_cs.csv')
s=f[(f.internal=='Cs-gluconate')&(f.drug.isna())&(f.protocol_type=='sPSC gap-free')].rename(columns={'date_header':'date'})
s['cellkey']=s.date.astype(int).astype(str)+'_s'+s.slice.astype(int).astype(str)+'_c'+s.cell.astype(int).astype(str)
g=ft.drop(columns=['session','date','slice','cell','holding_command_mV','duration_s','start_time']).merge(
  s[['new_name','cellkey','date','slice','cell','session','holding_command_mV','holding_flag','analysis_note','use_for_miniML','duration_s','start_time','T.Code','sex']],on='new_name')
g['t']=pd.to_datetime(g.date.astype(int).astype(str)+' '+g.start_time.astype(str),format='%Y%m%d %H:%M:%S')
cc['t']=pd.to_datetime(cc.time); cc['cellkey']=cc.date.astype(str)+'_s'+cc.slice.astype(str)+'_c'+cc.cell.astype(str)
rows=[]
for r in g.itertuples():
    d=dict(new_name=r.new_name); x=cc[(cc.cellkey==r.cellkey)&cc.Ra_MOhm.notna()]
    for side in 'ba':
        y=x[x.t<r.t].sort_values('t').tail(1) if side=='b' else x[x.t>r.t].sort_values('t').head(1)
        if len(y):
            y=y.iloc[0]; d.update({f'Ra_{side}':y.Ra_MOhm,f'Rm_{side}':y.Rm_MOhm,f'Cm_{side}':y.Cm_pF,f'gap_{side}_min':abs((r.t-y.t).total_seconds())/60})
    rows.append(d)
g=g.merge(pd.DataFrame(rows),on='new_name',how='left')
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
for nm,d in (('in','negative'),('out','positive')):
    x=mm[mm.direction==d].set_index('new_name')
    g['n_'+nm]=g.new_name.map(x.n_events); g['f_'+nm]=g.new_name.map(x.frequency_Hz); g['amp_'+nm]=g.new_name.map(x.amplitude_abs_median)
    g['q_'+nm]=g.new_name.map(x.charge_mean); g['valid_'+nm]=g.new_name.map(x.valid_time_s); g['noise_ml_'+nm]=g.new_name.map(x.noise_mad_pA)
a2=al.set_index('new_name')
for c in ['n_inward','n_outward','freq_inward_Hz','freq_outward_Hz']: g['d2_'+c]=g.new_name.map(a2[c])
g['abs_ihold']=g.Ihold_median_pA.abs(); g['abs_drift']=g.drift_total_pA.abs()
g.to_csv('cs_master_table.csv',index=False)
print(len(g)); print(g.groupby('holding_command_mV',dropna=False).agg(n=('new_name','size'),cells=('cellkey','nunique'),Ra=('Ra_b','median'),noise=('noise_lp3k_pA','median'),drift=('abs_drift','median'),ihold=('Ihold_median_pA','median'),fin=('f_in','median'),fout=('f_out','median'),ain=('amp_in','median'),aout=('amp_out','median'),tr=('n_transients','sum')).round(2))
