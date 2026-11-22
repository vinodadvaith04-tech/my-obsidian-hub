"""Assemble per-cell E and I frequencies for: paper controls, current lab (Helena, Sam, Advaith in Sam's experiment), Cs-gluconate (Advaith)."""
import pandas as pd, numpy as np, re
MASTER='/mnt/user-data/uploads/Ephys pilot & practice/Ephys_master.xlsx'
rows=[]
# ---- 1. paper (digitised from Fig 3C,G; Fig 6E,I)
pd_=pd.read_csv('paper_dots.csv'); pd_=pd_[pd_.group=='CTL']
for panel,g in pd_.groupby('panel'):
    div,kind=panel.split('_'); 
    for i,v in enumerate(g.value.clip(lower=0)):
        rows.append(dict(era='Past lab (paper)',source=f'Paper {div} control',measure='E' if kind=='sEPSC' else 'I',cell=f'{panel}_{i+1}',freq=float(v),freq_d2=np.nan,n_rec=np.nan,valid_s=np.nan,holding='-65'))
# ---- 2. current lab: miniML per cell/holding
m=pd.read_excel(MASTER,sheet_name='miniML_per_cell_holding')
f=pd.read_excel(MASTER,sheet_name='files')
a2=pd.read_excel(MASTER,sheet_name='alt_det_per_file')
a2['recording_name']=a2.new_name.str.replace(r'_\d+\.abf$','',regex=True)
d2=a2.groupby('recording_name').agg(f_in=('freq_inward_Hz','mean'),f_out=('freq_outward_Hz','mean')).reset_index()
src={'Helena KCl IPSC':('Current lab','Helena KCl IPSC (miniML)','I'),'Sam IPSC':('Current lab','Advaith KCl IPSC, Sam experiment (miniML)','I'),'Sam EPSC (22 Oct)':('Current lab','Advaith EPSC, Sam experiment (miniML)','E')}
for s,(era,name,meas) in src.items():
    y=m[(m.session==s)&(m.direction=='negative')].merge(d2,on='recording_name',how='left')
    for _,r in y.iterrows():
        rows.append(dict(era=era,source=name,measure=meas,cell=r.recording_name,freq=r.frequency_Hz,freq_d2=r.f_in,n_rec=r.n_files,valid_s=r.valid_time_s,holding=str(int(r.holding_command_mV))))
sam=pd.read_excel(MASTER,sheet_name='sam_EPSC_control_cells')
for _,r in sam.iterrows():
    rows.append(dict(era='Current lab',source="Sam's processed EPSC, control cells (his pipeline)",measure='E',cell=r.cell_id,freq=r.frequency_Hz_corrected,freq_d2=np.nan,n_rec=r.n_recordings_in_sam_sheet1,valid_s=r.valid_time_s,holding='-60'))
# ---- 3. Cs-gluconate (Advaith): from the Cs holdings analysis
cs=pd.read_csv('/home/claude/work/cs/cs_cell_holding.csv')
for _,r in cs[cs.holding_command_mV==-80].dropna(subset=['f_in']).iterrows():
    rows.append(dict(era='Cs-gluconate (Advaith)',source='Cs-gluconate, E at -80 mV (miniML)',measure='E',cell=r.cellkey,freq=r.f_in,freq_d2=r.d2_freq_inward_Hz,n_rec=r.n_rec,valid_s=np.nan,holding='-80'))
for _,r in cs[cs.holding_command_mV==0].dropna(subset=['f_out']).iterrows():
    rows.append(dict(era='Cs-gluconate (Advaith)',source='Cs-gluconate, I at 0 mV (miniML)',measure='I',cell=r.cellkey,freq=r.f_out,freq_d2=r.d2_freq_outward_Hz,n_rec=r.n_rec,valid_s=np.nan,holding='0'))
D=pd.DataFrame(rows); D.to_csv('lab_cells_long.csv',index=False)
print(D.groupby(['era','source','measure']).freq.agg(['count','mean','median','min','max']).round(2).to_string())
