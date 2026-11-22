import pandas as pd, numpy as np, json
from scipy import stats
pd.set_option('display.width',260,'display.max_rows',400,'display.max_columns',60)
a=pd.read_csv('k60_analysis_set.csv'); g_all=pd.read_csv('k60_master_table.csv')
cc=pd.read_csv('k60_capcheck.csv'); man=pd.read_csv('manual_checks_journal.csv')
a['cellkey']=a.date.astype(str)+'_s'+a.slice.astype(str)+'_c'+a.cell.astype(str)
a['cellkey']=a.cellkey.str.replace('.0','',regex=False)
cc['cellkey']=cc.date.astype(str)+'_s'+cc.slice.astype(str)+'_c'+cc.cell.astype(str)
# ---------- cell-level seal test from cap checks
cc['t']=pd.to_datetime(cc.time)
rows=[]
for k,d in cc[cc.Ra_MOhm.notna()].sort_values('t').groupby('cellkey'):
    r=dict(cellkey=k,n_cc=len(d),Ra_first=d.Ra_MOhm.iloc[0],Ra_last=d.Ra_MOhm.iloc[-1],Ra_max=d.Ra_MOhm.max(),Rm_first=d.Rm_MOhm.iloc[0],Rm_last=d.Rm_MOhm.iloc[-1],Cm_first=d.Cm_pF.iloc[0],Cm_last=d.Cm_pF.iloc[-1],span_min=(d.t.iloc[-1]-d.t.iloc[0]).total_seconds()/60)
    rows.append(r)
cell_cc=pd.DataFrame(rows)
cell_cc['dRa_pct']=(cell_cc.Ra_last-cell_cc.Ra_first)/cell_cc.Ra_first*100
cell_cc['dRm_pct']=(cell_cc.Rm_last-cell_cc.Rm_first)/cell_cc.Rm_first*100
cell_cc['dCm_pct']=(cell_cc.Cm_last-cell_cc.Cm_first)/cell_cc.Cm_first*100
cell_cc.loc[cell_cc.span_min<3,['dRa_pct','dRm_pct','dCm_pct']]=np.nan   # need >=3 min between first and last check
# ---------- manual checks (journal). master cell for 6 Nov = slice3 cell2
man['slice_m']=man.slice; man['cell_m']=man.cell
man.loc[man.date==20261106,['slice_m','cell_m']]=[3,2]
man['cellkey']=man.date.astype(str)+'_s'+man.slice_m.astype(str)+'_c'+man.cell_m.astype(str)
rows=[]
for k,d in man.groupby('cellkey'):
    d=d.sort_values('check'); 
    def fl(c):
        x=d[c].dropna(); return (x.iloc[0],x.iloc[-1],len(x)) if len(x) else (np.nan,np.nan,0)
    r=dict(cellkey=k)
    for c in ['Ra','Rm','Cm','Vm']:
        r[c+'_m1'],r[c+'_mL'],r['n_'+c]=fl(c)
    rows.append(r)
cell_man=pd.DataFrame(rows)
for c in ['Ra','Rm','Cm']:
    cell_man['d'+c+'_man_pct']=np.where(cell_man['n_'+c]>=2,(cell_man[c+'_mL']-cell_man[c+'_m1'])/cell_man[c+'_m1']*100,np.nan)
# ---------- cell table
cells=a.groupby('cellkey').agg(date=('date','first'),slice=('slice','first'),cell=('cell','first'),
    has_m60=('holding_command_mV',lambda x:(x==-60).any()),has_0=('holding_command_mV',lambda x:(x==0).any()),has_p=('holding_command_mV',lambda x:x.isin([5,10,15]).any())).reset_index()
cells=cells.merge(cell_cc,on='cellkey',how='left').merge(cell_man,on='cellkey',how='left')
cells.to_csv('k60_cells.csv',index=False)
print(len(cells), cells.has_m60.sum(), cells.has_0.sum(), cells.has_p.sum(), (cells.has_m60&cells.has_0).sum(), (cells.has_m60&cells.has_p).sum())
def ms(x,d=1):
    x=pd.Series(x).dropna(); 
    if len(x)==0: return 'n/a'
    return f"{x.mean():.{d}f} ± {x.std(ddof=1)/np.sqrt(len(x)) if len(x)>1 else float('nan'):.{d}f}  | {x.median():.{d}f} [{x.quantile(.25):.{d}f}, {x.quantile(.75):.{d}f}] (n={len(x)})"
for nm,sel in [('cells at -60',cells.has_m60),('cells at 0',cells.has_0),('cells +5-15',cells.has_p)]:
    c=cells[sel]; print('=====',nm,len(c))
    for col in ['n_cc','Ra_first','Ra_last','Ra_max','dRa_pct','Rm_first','Rm_last','dRm_pct','Cm_first','Cm_last','dCm_pct','span_min','Ra_m1','Ra_mL','dRa_man_pct','Rm_m1','Rm_mL','dRm_man_pct','Cm_m1','Cm_mL','dCm_man_pct','Vm_m1']:
        print('  ',col,ms(c[col]))
    print('   Ra>35 at any check:',(c.Ra_max>35).sum(),'/',c.Ra_max.notna().sum(), ' dRa>25%:',(c.dRa_pct>25).sum(),'/',c.dRa_pct.notna().sum())
# recording-level Ra before
a['QC_Ra35']=a.Ra_b<35
print(a.groupby('grp').agg(n=('Ra_b','count'),QC=('QC_Ra35','sum'),Ra30=('Ra_b',lambda x:(x<=30).sum())))
a.to_csv('k60_analysis_set2.csv',index=False)
