"""Per-cell event measures for the lab-standards comparison (26 Nov 2026).
paper: digitised dots (Fig 3 D,E,F,H,I,J; Fig 6 F,G,H,J,K,L); Sam: his miniML per-cell values (amplitude = his per-cell median column; rise, charge, half-decay = his per-cell means, the only statistic exported); miniML groups: per-cell MEDIAN over all events of the cell's recordings at the holding."""
import pandas as pd, numpy as np
M='/mnt/user-data/uploads/Ephys pilot & practice/Ephys_master.xlsx'
L=pd.read_csv('/home/claude/work/lab/lab_cells_long.csv')            # same cells as the frequency figure
ev=pd.read_pickle('/home/claude/work/lab3/ev.pkl'); f=pd.read_excel(M,sheet_name='files')
ev['cellkey_rec']=ev.recording_name
rows=[]
# ---- paper
pdm=pd.read_csv('/home/claude/work/lab4/paper_dots_more.csv'); pdm=pdm[pdm.group=='CTL']
mm={'amp':'amplitude','rise':'risetime','decay':'decay_tau'}
for panel,g in pdm.groupby('panel'):
    div,kind,meas=panel.split('_'); 
    for i,v in enumerate(g.value.values):
        rows.append(dict(source=f'Paper {div} control',measure='E' if kind=='sEPSC' else 'I',variable=mm[meas],cell=f'{panel}_{i+1}',value=float(v)))
# ---- Sam
sam=pd.read_excel(M,sheet_name='sam_EPSC_control_cells')
for _,r in sam.iterrows():
    # Sam's batch notebook (miniML, GC_lstm_model, same settings as ours) saves decay_mean_ms = mean of miniML half-decay x 1000, so it is half-decay
    for var,col in (('amplitude','amplitude_abs_median'),('risetime','risetime_mean_ms'),('charge','charge_mean_pC'),('halfdecay','decay_mean_ms')):
        v=r[col]; rows.append(dict(source="Sam's processed EPSC, control cells (his pipeline)",measure='E',variable=var,cell=r.cell_id,value=abs(v) if var=='charge' else v))
# ---- miniML groups: cells exactly as in the frequency table
def stat(e,var,col):
    # like-for-like with published per-cell MEANS (paper rise; Sam rise, charge, half-decay): 99%-trimmed mean (top 1% of the cell's events dropped as artefact guard)
    # amplitude (Sam exports a per-cell median) and half-width (no published comparator): median
    x=e[col].abs().dropna()
    if var in ('amplitude','halfwidth'): return float(x.median())
    return float(x[x<=x.quantile(.99)].mean())
def cellev(recname,holding,direction):
    e=ev[(ev.recording_name==recname)&(ev.holding_command_mV==holding)&(ev.direction==direction)]
    return e
for src,meas,dirn in (('Helena KCl IPSC (miniML)','I','negative'),('Advaith KCl IPSC, Sam experiment (miniML)','I','negative'),('Advaith EPSC, Sam experiment (miniML)','E','negative')):
    for _,r in L[L.source==src].iterrows():
        e=cellev(r.cell,int(r.holding),dirn)
        if len(e)==0: print('no events',src,r.cell); continue
        for var,col in (('amplitude','amplitude_abs'),('risetime','risetime_ms'),('charge','charge'),('halfdecay','halfdecay_ms'),('halfwidth','halfwidth_ms')):
            rows.append(dict(source=src,measure=meas,variable=var,cell=r.cell,value=stat(e,var,col),n_events=len(e)))
# Cs: events of the cell's recordings at the holding (cellkey = date_slice_cell)
cs=pd.read_csv('/home/claude/work/cs/cs_analysis_set.csv')
for src,meas,h,dirn in (('Cs-gluconate, E at -80 mV (miniML)','E',-80,'negative'),('Cs-gluconate, I at 0 mV (miniML)','I',0,'positive')):
    for ck,g in cs[cs.holding_command_mV==h].groupby('cellkey'):
        if src.startswith('Cs-gluconate, E') and not (g.f_in.notna().any()): continue
        if src.startswith('Cs-gluconate, I') and not (g.f_out.notna().any()): continue
        e=ev[ev.new_name.isin(g.new_name)&(ev.direction==dirn)]
        if len(e)==0: print('no events',ck,h); continue
        for var,col in (('amplitude','amplitude_abs'),('risetime','risetime_ms'),('charge','charge'),('halfdecay','halfdecay_ms'),('halfwidth','halfwidth_ms')):
            rows.append(dict(source=src,measure=meas,variable=var,cell=ck,value=stat(e,var,col),n_events=len(e)))
D=pd.DataFrame(rows); D.to_csv('lab_measures_long.csv',index=False)
print(D.groupby(['variable','source','measure']).value.agg(['count','median','min','max']).round(2).to_string())
