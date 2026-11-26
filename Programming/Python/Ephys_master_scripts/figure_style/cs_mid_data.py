import sys; sys.path.insert(0,'/home/claude/work/style')
import numpy as np, pandas as pd
from scipy import stats
C='/home/claude/work/cs/'
a=pd.read_csv(C+'cs_analysis_set.csv')          # same per-cell set as Figure 3 and Options A-E (nothing left out)
cm=lambda c: a.groupby(['cellkey','holding_command_mV'])[c].mean().unstack()
fi,fo=cm('f_in'),cm('f_out')
MID=[-60,-50,-40,-30]
base=pd.DataFrame({'E0':fi[-80],'I0':fo[0]}).dropna(); base=base[(base.E0>0)&(base.I0>0)]
base['ratio_ref']=base.E0/base.I0; base['share_ref']=base.E0/(base.E0+base.I0)
rows=[]
for h in MID:
    t=pd.DataFrame({'fi':fi[h],'fo':fo[h]}).join(base,how='inner').dropna(); t=t[(t.fi>0)&(t.fo>0)]
    t['h']=h; t['ratio_h']=t.fi/t.fo; t['fold']=t.ratio_h/t.ratio_ref; t['share_h']=t.fi/(t.fi+t.fo); t['dshare']=(t.share_h-t.share_ref)*100
    t['abserr']=t.dshare.abs(); t['absl2']=np.abs(np.log2(t.fold)); rows.append(t.reset_index().rename(columns={'cellkey':'cell'}))
pc=pd.concat(rows,ignore_index=True)
# interpolated crossing holding per cell (as in the report: linear in log2 ratio, holdings -60..0)
cross={}
for c in base.index:
    xs,ys=[],[]
    for h in [-60,-50,-40,-30,0]:
        v1,v2=fi.at[c,h] if h in fi.columns else np.nan, fo.at[c,h] if h in fo.columns else np.nan
        if pd.notna(v1) and pd.notna(v2) and v1>0 and v2>0: xs.append(h); ys.append(np.log2(v1/v2)-np.log2(base.at[c,'E0']/base.at[c,'I0']))
    for k in range(len(xs)-1):
        if ys[k]>=0>=ys[k+1] and ys[k]!=ys[k+1]: cross[c]=xs[k]+(0-ys[k])/(ys[k+1]-ys[k])*(xs[k+1]-xs[k]); break
def short(c):
    d,s,cl=c.split('_'); mon={'11':'Nov'}[d[4:6]]; return f"{int(d[6:])} Nov {s} {cl}"
if __name__=='__main__':
    print(len(base),'cells with both refs'); print(pc.groupby('h').fold.agg(['size','median']).round(2))
    print('cross n',len(cross),'median',np.median(list(cross.values())).round(1),np.round(sorted(cross.values()),1))
    print(pc.groupby('h')[['dshare','abserr']].median().round(1))
