import pandas as pd, numpy as np, json
from scipy import stats
pd.set_option('display.width',250,'display.max_columns',80,'display.max_rows',300)
s=pd.read_pickle('/tmp/mid_sel.pkl')          # selected files rows from master (files sheet)
mm=pd.read_pickle('/tmp/miniML_per_file_new.pkl'); al=pd.read_pickle('/tmp/alt_det_per_file_new.pkl')
ft=pd.read_csv('k60_gapfree_features_mid.csv'); cc=pd.read_csv('k60_capcheck_mid.csv'); tr=pd.read_csv('k60_transients_mid.csv')
man=pd.read_csv('../manual_checks_journal.csv')
s=s.rename(columns={'date_header':'date'})
s['cellkey']=s.date.astype(str)+'_s'+s.slice.astype(int).astype(str)+'_c'+s.cell.astype(int).astype(str)
g=ft.drop(columns=['session','date','slice','cell','holding_command_mV','duration_s','start_time']).merge(
  s[['new_name','cellkey','date','slice','cell','session','holding_command_mV','holding_listed_in_note','ljp_in_note','holding_flag','duration_s','start_time','T.Code','sex']],on='new_name')
g['t']=pd.to_datetime(g.date.astype(int).astype(str)+' '+g.start_time.astype(str),format='%Y%m%d %H:%M:%S')
# ---- cap checks nearest before / after
cc['t']=pd.to_datetime(cc.time); cc['cellkey']=cc.date.astype(str)+'_s'+cc.slice.astype(str)+'_c'+cc.cell.astype(str)
rows=[]
for r in g.itertuples():
    d=dict(new_name=r.new_name); x=cc[(cc.cellkey==r.cellkey)&cc.Ra_MOhm.notna()]
    for side in 'ba':
        y=x[x.t<r.t].sort_values('t').tail(1) if side=='b' else x[x.t>r.t].sort_values('t').head(1)
        if len(y):
            y=y.iloc[0]; d.update({f'Ra_{side}':y.Ra_MOhm,f'Rm_{side}':y.Rm_MOhm,f'Cm_{side}':y.Cm_pF,f'V_{side}':y.Vhold_capcheck_mV,
              f'gap_{side}_min':abs((r.t-y.t).total_seconds())/60})
    rows.append(d)
g=g.merge(pd.DataFrame(rows),on='new_name',how='left')
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
g['max_transient_pA']=g.new_name.map(ep.groupby('new_name').amp.max())
# ---- miniML both directions + detector 2
for nm,d in (('in','negative'),('out','positive')):
    x=mm[mm.direction==d].set_index('new_name')
    g['n_'+nm]=g.new_name.map(x.n_events); g['f_'+nm]=g.new_name.map(x.frequency_Hz); g['amp_'+nm]=g.new_name.map(x.amplitude_abs_median)
    g['valid_'+nm]=g.new_name.map(x.valid_time_s); g['noise_ml_'+nm]=g.new_name.map(x.noise_mad_pA)
a2=al.set_index('new_name')
for c in ['n_inward','n_outward','freq_inward_Hz','freq_outward_Hz']: g['d2_'+c]=g.new_name.map(a2[c])
g['abs_ihold']=g.Ihold_median_pA.abs(); g['abs_drift']=g.drift_total_pA.abs()
g.to_csv('mid_master_table.csv',index=False)
# ---- cells
cc2=cc[cc.Ra_MOhm.notna()].sort_values('t')
cr=[]
for k,d in cc2.groupby('cellkey'):
    span=(d.t.iloc[-1]-d.t.iloc[0]).total_seconds()/60
    cr.append(dict(cellkey=k,n_cc=len(d),Ra_first=d.Ra_MOhm.iloc[0],Ra_last=d.Ra_MOhm.iloc[-1],Ra_max=d.Ra_MOhm.max(),Rm_first=d.Rm_MOhm.iloc[0],Rm_last=d.Rm_MOhm.iloc[-1],Cm_first=d.Cm_pF.iloc[0],Cm_last=d.Cm_pF.iloc[-1],span_min=span))
cell_cc=pd.DataFrame(cr)
for c in ['Ra','Rm','Cm']:
    cell_cc['d'+c+'_pct']=np.where(cell_cc.span_min>=3,(cell_cc[c+'_last']-cell_cc[c+'_first'])/cell_cc[c+'_first']*100,np.nan)
man['cellkey']=man.date.astype(str)+'_s'+man.slice.astype(str)+'_c'+man.cell.astype(str)
mr=[]
for k,d in man.groupby('cellkey'):
    d=d.sort_values('check'); r=dict(cellkey=k)
    for c in ['Ra','Rm','Cm','Vm']:
        x=d[c].dropna(); r[c+'_m1']=x.iloc[0] if len(x) else np.nan; r[c+'_mL']=x.iloc[-1] if len(x) else np.nan; r['n_'+c]=len(x)
    mr.append(r)
cell_man=pd.DataFrame(mr)
for c in ['Ra','Rm','Cm']:
    cell_man['d'+c+'_man_pct']=np.where(cell_man['n_'+c]>=2,(cell_man[c+'_mL']-cell_man[c+'_m1'])/cell_man[c+'_m1']*100,np.nan)
H=[-50,-45,-40,-35,-30]
cells=g.groupby('cellkey').agg(date=('date','first'),slice=('slice','first'),cell=('cell','first'),n_rec=('new_name','size'),holdings=('holding_command_mV',lambda x:sorted(set(x)))).reset_index()
for h in H: cells[f'has_{abs(h)}']=cells.cellkey.map(g[g.holding_command_mV==h].groupby('cellkey').size()).fillna(0)>0
cells=cells.merge(cell_cc,on='cellkey',how='left').merge(cell_man,on='cellkey',how='left')
cells['holdings']=cells.holdings.astype(str)
cells.to_csv('mid_cells.csv',index=False)
print(len(g),len(cells)); print(g.groupby('holding_command_mV').agg(n=('new_name','size'),cells=('cellkey','nunique'),Ra_b=('Ra_b','median'),noise=('noise_lp3k_pA','median'),drift=('abs_drift','median'),ihold=('Ihold_median_pA','median'),nin=('n_in','median'),nout=('n_out','median')).round(1))
print(cells[['cellkey','n_rec','holdings','n_cc','Ra_first','Ra_last','Ra_max','Ra_m1','Ra_mL']].round(1).to_string())
