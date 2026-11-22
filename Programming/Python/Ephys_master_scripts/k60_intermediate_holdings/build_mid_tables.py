import pandas as pd, numpy as np, json
from scipy import stats
g=pd.read_csv('mid_master_table.csv'); cells=pd.read_csv('mid_cells.csv')
a60=pd.read_csv('../k60_analysis_set3.csv'); a60=a60[a60.grp=='-60 mV']; c60=pd.read_csv('../k60_cells.csv'); c60=c60[c60.has_m60]
H=[-50,-45,-40,-35,-30]; cols=[str(h) for h in H]+['All -30 to -50','-60 (reference)']
D={str(h):g[g.holding_command_mV==h] for h in H}; D['All -30 to -50']=g; D['-60 (reference)']=a60
a60=a60.copy(); a60['n_in']=a60.n_neg; a60['f_in']=a60.f_neg; a60['amp_in']=a60.amp_neg
D['-60 (reference)']=a60
def q(x,d=0):
    x=pd.Series(x).dropna()
    return 'n/a' if len(x)==0 else f"{x.median():.{d}f} [{x.quantile(.25):.{d}f}, {x.quantile(.75):.{d}f}]"
def mse(x,d=1):
    x=pd.Series(x).dropna()
    return 'n/a' if len(x)==0 else (f"{x.mean():.{d}f} ± {x.std(ddof=1)/np.sqrt(len(x)):.{d}f}" if len(x)>1 else f"{x.mean():.{d}f}")
def frac(m): m=pd.Series(m).dropna().astype(bool); return f"{int(m.sum())}/{len(m)} ({m.mean()*100:.0f}%)" if len(m) else 'n/a'
def ck(x): return x.cellkey
out={}
def R(name,fn,tab):
    tab.append([name]+[fn(D[c],c) for c in cols])
# table 1
t=[]
R('Recording days',lambda x,c:str(x.date.nunique()),t)
R('Slices',lambda x,c:str(x[['date','slice']].drop_duplicates().shape[0]),t)
R('Cells',lambda x,c:str(x.cellkey.nunique()),t)
R('Gap-free recordings',lambda x,c:str(len(x)),t)
R('Total recording time (min)',lambda x,c:f"{x.duration_s.sum()/60:.1f}",t)
R('Recordings per cell, median (range)',lambda x,c:f"{x.groupby('cellkey').size().median():.0f} ({x.groupby('cellkey').size().min()}-{x.groupby('cellkey').size().max()})",t)
R('Holding written in the journal',lambda x,c:str(len(x)) if c!='-60 (reference)' else str((x.label_status=='confirmed (journal/stated)').sum()),t)
R('Of these, from LJP-corrected sessions (7 and 8 Oct)',lambda x,c:str(int((x.ljp_in_note.fillna('').str.startswith('+10')).sum())) if c!='-60 (reference)' else 'n/a',t)
out['t1']=t
# table 2a
t=[]
R('Recordings with a cap check before them',lambda x,c:str(x.Ra_b.notna().sum()),t)
R('Ra before the recording (MOhm), mean ± SEM',lambda x,c:mse(x.Ra_b),t)
R('Rm before the recording (MOhm)',lambda x,c:mse(x.Rm_b),t)
R('Cm before the recording (pF)',lambda x,c:mse(x.Cm_b),t)
R('Ra < 35 MOhm before the recording',lambda x,c:frac(x.Ra_b.dropna()<35),t)
R('Minutes from that cap check to the recording, median',lambda x,c:f"{x.gap_b_min.median():.0f}" if 'gap_b_min' in x and x.gap_b_min.notna().any() else 'n/a',t)
R('Recordings with a cap check after them',lambda x,c:str(x.Ra_a.notna().sum()),t)
out['t2a']=t
# table 2b (cells that have that holding)
cg={str(h):cells[cells[f'has_{abs(h)}']] for h in H}; cg['All -30 to -50']=cells; cg['-60 (reference)']=c60
t=[]
def C(name,fn): t.append([name]+[fn(cg[c]) for c in cols])
C('Cells with journal manual checks',lambda c:str(c.Ra_m1.notna().sum()))
C('Ra, first to last check (MOhm)',lambda c:f"{c.Ra_m1.mean():.1f} to {c.Ra_mL.mean():.1f}")
C('Ra change, median [IQR] (%)',lambda c:q(c.dRa_man_pct))
C('Rm, first to last check (MOhm)',lambda c:f"{c.Rm_m1.mean():.0f} to {c.Rm_mL.mean():.0f}")
C('Cm, first to last check (pF)',lambda c:f"{c.Cm_m1.mean():.0f} to {c.Cm_mL.mean():.0f}")
C('Membrane potential at first check (mV)',lambda c:mse(c.Vm_m1))
C('Cells with Ra > 35 MOhm at any cap check',lambda c:frac(c.Ra_max.dropna()>35))
out['t2b']=t
# table 3
t=[]
R('Recordings',lambda x,c:str(len(x)),t)
R('Holding current (pA), median [IQR]',lambda x,c:q(x.Ihold_median_pA),t)
R('|Holding current| > 100 pA',lambda x,c:frac(x.abs_ihold>100),t)
R('Drift over 2 min, absolute (pA)',lambda x,c:q(x.abs_drift,1),t)
R('|Drift| > 25 pA',lambda x,c:frac(x.abs_drift>25),t)
R('Baseline wander, SD about trend (pA)',lambda x,c:q(x.wander_sd_pA,1),t)
R('Noise, MAD-SD (pA)',lambda x,c:q(x.noise_lp3k_pA,1),t)
R('Noise > 4 pA',lambda x,c:frac(x.noise_lp3k_pA>4),t)
R('Recordings with at least 1 spike-like transient',lambda x,c:frac(x.n_transients>0),t)
R('Spike-like transients per min',lambda x,c:f"{x.n_transients.sum()/(x.duration_s.sum()/60):.2f} (n = {int(x.n_transients.sum())})",t)
out['t3']=t
# table 4: events
def evtab(d,dd2,ref):
    t=[]
    n,f,am=f'n_{d}',f'f_{d}',f'amp_{d}'
    cf=lambda x:x.groupby('cellkey')[f].mean()
    def okref(x,c): return c!='-60 (reference)' or ref
    R('Recordings analysed',lambda x,c:str(len(x)) if okref(x,c) else 'not analysed',t)
    R('Total miniML events',lambda x,c:f"{int(x[n].sum()):,}" if okref(x,c) else 'n/a',t)
    R('Events per recording, median [IQR]',lambda x,c:q(x[n]) if okref(x,c) else 'n/a',t)
    R('Events per recording, range',lambda x,c:f"{int(x[n].min())}-{int(x[n].max())}" if okref(x,c) else 'n/a',t)
    R('Frequency per recording (Hz), median [IQR]',lambda x,c:q(x[f],2) if okref(x,c) else 'n/a',t)
    R('Frequency per cell (Hz), median [IQR]',lambda x,c:q(cf(x),2) if okref(x,c) else 'n/a',t)
    R('Pooled frequency (events / recorded time, Hz)',lambda x,c:f"{x[n].sum()/x.duration_s.sum():.2f}" if okref(x,c) else 'n/a',t)
    R('Recordings with < 5 events',lambda x,c:frac(x[n]<5) if okref(x,c) else 'n/a',t)
    R('Recordings with no event',lambda x,c:str(int((x[n]==0).sum())) if okref(x,c) else 'n/a',t)
    R('Event amplitude (pA), median of recording medians',lambda x,c:f"{x[am].median():.1f} (n = {x[am].notna().sum()})" if okref(x,c) else 'n/a',t)
    R('Amplitude / recording noise, median',lambda x,c:f"{(x[am]/x.noise_lp3k_pA).median():.1f}" if okref(x,c) else 'n/a',t)
    R('Second detector, same direction: events per recording',lambda x,c:q(x[dd2],0),t)
    return t
out['t4e']=evtab('in','d2_n_inward',True) ; out['t4i']=evtab('out','d2_n_outward',False)
# 60 ref second detector needs column names
for k in ('t4e','t4i'): pass
json.dump(out,open('mid_tables.json','w'),indent=1)
for k,v in out.items(): print('\n==',k); print(' | '.join(['']+cols)); print('\n'.join(' | '.join(r) for r in v))
