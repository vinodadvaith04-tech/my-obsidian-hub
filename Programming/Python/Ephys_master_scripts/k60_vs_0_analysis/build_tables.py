import pandas as pd, numpy as np, json
a=pd.read_csv('k60_analysis_set3.csv'); cells=pd.read_csv('k60_cells.csv'); cs=pd.read_csv('k60_gapfree_features_cs60.csv')
a['abs_drift']=a.drift_total_pA.abs(); a['abs_ihold']=a.Ihold_median_pA.abs()
G=['-60 mV','0 mV','+5 to +15 mV']
def q(x,d=0): 
    x=pd.Series(x).dropna(); return f"{x.median():.{d}f} [{x.quantile(.25):.{d}f}, {x.quantile(.75):.{d}f}]"
def mse(x,d=1):
    x=pd.Series(x).dropna(); return f"{x.mean():.{d}f} ± {x.std(ddof=1)/np.sqrt(len(x)):.{d}f}"
def frac(mask): return f"{int(mask.sum())}/{len(mask)} ({mask.mean()*100:.0f}%)"
out={}
# ---- table1
rows=[]
d={g:a[a.grp==g] for g in G}
def rowfn(name,fn): rows.append([name]+[fn(d[g],g) for g in G])
rowfn('Recording days',lambda x,g:str(x.date.nunique()))
rowfn('Slices',lambda x,g:str(x[['date','slice']].drop_duplicates().shape[0]))
rowfn('Cells',lambda x,g:str(x.cellkey.nunique()))
rowfn('Gap-free recordings',lambda x,g:str(len(x)))
rowfn('Total recording time (min)',lambda x,g:f"{x.duration_s.sum()/60:.1f}")
rowfn('Recordings per cell, median (range)',lambda x,g:f"{x.groupby('cellkey').size().median():.0f} ({x.groupby('cellkey').size().min()}–{x.groupby('cellkey').size().max()})")
rowfn('Holding written in journal or stated by you',lambda x,g:f"{(x.label_status=='confirmed (journal/stated)').sum()}")
rowfn('Holding assumed (default), cap checks at −60 mV',lambda x,g:f"{(x.label_status=='assumed, cap checks at -60').sum()}")
rowfn('Holding assumed (default), unverifiable',lambda x,g:f"{(x.label_status=='assumed, unverified').sum()}")
out['t1']=rows
# ---- table2
rows=[]
rowfn('Recordings with a cap check before them',lambda x,g:str(x.Ra_b.notna().sum()))
rowfn('Ra before the recording (MΩ), mean ± SEM',lambda x,g:mse(x.Ra_b))
rowfn('Rm before the recording (MΩ)',lambda x,g:mse(x.Rm_b))
rowfn('Cm before the recording (pF)',lambda x,g:mse(x.Cm_b))
rowfn('Ra < 35 MΩ before the recording',lambda x,g:frac(x.Ra_b.dropna()<35))
rowfn('Recordings with a cap check after them',lambda x,g:str(x.Ra_a.notna().sum()))
out['t2a']=rows
cg={'-60 mV':cells[cells.has_m60],'0 mV':cells[cells.has_0],'+5 to +15 mV':cells[cells.has_p]}
rows=[]
def crow(name,fn): rows.append([name]+[fn(cg[g]) for g in G])
crow('Cells with journal manual checks',lambda c:str(c.Ra_m1.notna().sum()))
crow('Ra, first → last check (MΩ)',lambda c:f"{c.Ra_m1.mean():.1f} → {c.Ra_mL.mean():.1f}")
crow('Ra change, median [IQR] (%)',lambda c:q(c.dRa_man_pct))
crow('Rm, first → last check (MΩ)',lambda c:f"{c.Rm_m1.mean():.0f} → {c.Rm_mL.mean():.0f}")
crow('Cm, first → last check (pF)',lambda c:f"{c.Cm_m1.mean():.0f} → {c.Cm_mL.mean():.0f}")
crow('Membrane potential at first check (mV)',lambda c:mse(c.Vm_m1))
crow('Cells with Ra > 35 MΩ at any cap check',lambda c:frac(c.Ra_max.dropna()>35))
out['t2b']=rows
# ---- table3
rows=[]
rowfn('Recordings',lambda x,g:str(len(x)))
rowfn('Holding current (pA), median [IQR]',lambda x,g:q(x.Ihold_median_pA))
rowfn('|Holding current| > 100 pA',lambda x,g:frac(x.abs_ihold>100))
rowfn('Drift over 2 min, absolute (pA)',lambda x,g:q(x.abs_drift,1))
rowfn('|Drift| > 25 pA',lambda x,g:frac(x.abs_drift>25))
rowfn('Baseline wander, SD about trend (pA)',lambda x,g:q(x.wander_sd_pA,1))
rowfn('Noise, MAD-SD (pA)',lambda x,g:q(x.noise_lp3k_pA,1))
rowfn('Noise > 4 pA',lambda x,g:frac(x.noise_lp3k_pA>4))
rowfn('Recordings with ≥1 spike-like transient',lambda x,g:frac(x.n_transients>0))
rowfn('Spike-like transients per min',lambda x,g:f"{x.n_transients.sum()/(x.duration_s.sum()/60):.2f} (n = {x.n_transients.sum()})")
# comparator
cs['abs_drift']=cs.drift_total_pA.abs(); cs['abs_ihold']=cs.Ihold_median_pA.abs()
tc=pd.read_csv('k60_transients_cs60.csv'); tc=tc[(tc.amp_pA>=200)&(tc.fwhm_ms<=3)]
comp=[str(len(cs)),q(cs.Ihold_median_pA),frac(cs.abs_ihold>100),q(cs.abs_drift,1),frac(cs.abs_drift>25),q(cs.wander_sd_pA,1),q(cs.noise_lp3k_pA,1),frac(cs.noise_lp3k_pA>4)]
out['t3']=rows; out['t3_comp']=comp
# ---- table4
rows=[]
rowfn('Direction analysed',lambda x,g:'Inward (EPSC)' if g=='-60 mV' else 'Outward (IPSC)')
rowfn('Recordings',lambda x,g:str(len(x)))
rowfn('Total miniML events',lambda x,g:f"{int(x.event_n.sum()):,}")
rowfn('Events per recording, median [IQR]',lambda x,g:q(x.event_n,0))
rowfn('Events per recording, range',lambda x,g:f"{int(x.event_n.min())}–{int(x.event_n.max())}")
rowfn('Frequency (Hz), median [IQR]',lambda x,g:q(x.event_f,2))
rowfn('Pooled frequency (events / recorded time, Hz)',lambda x,g:f"{x.event_n.sum()/x.duration_s.sum():.2f}")
rowfn('Recordings with < 5 events',lambda x,g:frac(x.event_n<5))
rowfn('Recordings with no event',lambda x,g:str(int((x.event_n==0).sum())))
rowfn('Event amplitude (pA), median of recording medians',lambda x,g:f"{x.event_amp.median():.1f} (n = {x.event_amp.notna().sum()})")
rowfn('Second detector, same direction: events per recording',lambda x,g:q(x.d2_n_inward if g=='-60 mV' else x.d2_n_outward,0))
out['t4']=rows
c60=a[a.grp=='-60 mV'].groupby('cellkey').agg(ev=('event_n','sum'),dur=('duration_s','sum')); c60['f']=c60.ev/c60.dur
out['c60']=dict(n=len(c60),ge01=int((c60.f>=0.1).sum()),ge03=int((c60.f>=0.3).sum()),lt005=int((c60.f<0.05).sum()),med=c60.f.median(),q1=c60.f.quantile(.25),q3=c60.f.quantile(.75))
for k,v in out.items(): print(k); print(v if not isinstance(v,list) else '\n'.join(' | '.join(r) for r in v))
json.dump(out,open('tables.json','w'),indent=1,default=str)
