import pandas as pd, numpy as np
from scipy import stats
pd.set_option('display.width',250,'display.max_columns',80,'display.max_rows',400)
g=pd.read_csv('cs_master_table.csv')
a=g[(g.use_for_miniML=='Y')&g.holding_command_mV.notna()&(g.duration_s>=30)].copy()
print(len(a),'recordings used;',a.cellkey.nunique(),'cells')
V=['f_in','f_out','amp_in','amp_out','q_in','q_out','noise_lp3k_pA','abs_drift','Ihold_median_pA','d2_freq_inward_Hz','d2_freq_outward_Hz','Ra_b','n_transients']
ch=a.groupby(['cellkey','holding_command_mV'])[V].mean().reset_index()
ch['n_rec']=a.groupby(['cellkey','holding_command_mV']).size().values
ch.to_csv('cs_cell_holding.csv',index=False)
H=sorted(ch.holding_command_mV.unique())
def piv(v): return ch.pivot(index='cellkey',columns='holding_command_mV',values=v)
fi,fo=piv('f_in'),piv('f_out')
print('cells per holding',(fi.notna().sum()).to_dict())
print(pd.concat([fi.median().rename('f_in med'),fo.median().rename('f_out med'),fi.mean().rename('f_in mean'),fo.mean().rename('f_out mean')],axis=1).round(3))
print(fi.round(2).to_string()); print(fo.round(2).to_string())
# paired cells for ref combos
for e,i in [(-80,0),(-70,0),(-80,10),(-70,10)]:
    print(e,i,'cells with both',(fi[e].notna()&fo[i].notna()).sum())
