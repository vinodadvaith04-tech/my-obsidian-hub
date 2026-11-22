import pandas as pd, numpy as np, json
from scipy import stats
a=pd.read_csv('cs_analysis_set.csv'); P=json.load(open('cs_tables_part1.json')); Q=json.load(open('cs_tables_part2.json'))['primary']
H=[-80,-70,-60,-50,-40,-30,-20,0,10,15]
POT={-80:-96,-70:-86,-60:-76,-50:-66,-40:-56,-30:-46,-20:-36,0:-16,10:-6,15:-1}
def s(x,d=1): return f'{x:.{d}f}'.replace('-','−')
def mi(t,d=1):
    m,q1,q3=t
    return s(m,d) if (q1==q3) else f'{s(m,d)} [{s(q1,d)}, {s(q3,d)}]'
def ms(t,d=2): return s(t[0],d) if (t[1]!=t[1]) else f'{s(t[0],d)} ± {s(t[1],d)}'
def hl(h): return ('+' if h>0 else '−' if h<0 else '')+str(abs(h)) if h!=0 else '0'
def one(rowcells): return rowcells==1
tabs={}
# T1
L=['| Command holding (mV) | Approx. potential (mV) | Days | Slices | Cells | Recordings | Recording time (min) |','|---|---|---|---|---|---|---|']
for r in P['t1']:
    L.append(f"| {hl(r['h'])} | {s(POT[r['h']],0)} | {r['days']} | {r['slices']} | {r['cells']} | {r['recs']} | {s(r['minutes'],1)} |")
nsl=a[['date','slice']].drop_duplicates().shape[0]
L.append(f"| **All** | | {a.date.nunique()} | {nsl} | {a.cellkey.nunique()} | {len(a)} | {s(a.duration_s.sum()/60,1)} |")
tabs['t1']='\n'.join(L)
# T2
L=['| Command holding (mV) | Cells | Ra (MΩ) | Rm (MΩ) | Cm (pF) | Recordings with Ra ≥ 35 MΩ |','|---|---|---|---|---|---|']
for r in P['t2']:
    L.append(f"| {hl(r['h'])} | {r['cells']} | {mi(r['Ra'])} | {mi(r['Rm'],0)} | {mi(r['Cm'],0)} | {r['ra35_rec']} of {r['n_rec']} |")
tabs['t2']='\n'.join(L)
# T3
L=['| Command holding (mV) | Cells | Noise (pA) | Drift over 2 min (pA) | Holding current (pA) | Spike-like transients per recording |','|---|---|---|---|---|---|']
for r in P['t3']:
    L.append(f"| {hl(r['h'])} | {r['cells']} | {mi(r['noise'])} | {mi(r['drift'],0)} | {mi(r['ihold'],0)} | {mi(r['transients_per_rec'],1)} |")
tabs['t3']='\n'.join(L)
# T4
L=['| Command holding (mV) | Cells | Inward (EPSC) frequency, Hz | Inward, mean ± SEM | Outward (IPSC) frequency, Hz | Outward, mean ± SEM | Inward amplitude (pA) | Outward amplitude (pA) |','|---|---|---|---|---|---|---|---|']
for r in P['t4']:
    L.append(f"| {hl(r['h'])} | {r['cells']} | {mi(r['f_in'],2)} | {ms(r['f_in_mean'])} | {mi(r['f_out'],2)} | {ms(r['f_out_mean'])} | {s(r['amp_in'][0],1)} | {s(r['amp_out'][0],1)} |")
tabs['t4']='\n'.join(L)
# T5 purity
L=['| Command holding (mV) | Cells | Inward share of events, % | Outward share of events, % | Inward share, second detector, % |','|---|---|---|---|---|']
for r in Q['purity']:
    L.append(f"| {hl(int(r['h']))} | {r['n']} | {s(r['pctE'],0)} [{s(r['pctE_q'][0],0)}, {s(r['pctE_q'][1],0)}] | {s(r['pctI'],0)} [{s(r['pctI_q'][0],0)}, {s(r['pctI_q'][1],0)}] | {s(r['d2_pctE'],0)} |")
tabs['t5']='\n'.join(L)
# T6
e=Q['ei']['-80/0']; L=['| Command holding (mV) | Cells | EPSC frequency retained, % of that at −80 | IPSC frequency retained, % of that at 0 | E/I ratio at this holding ÷ reference ratio | Cells within ×1.5 | Wilcoxon p |','|---|---|---|---|---|---|---|']
for r in e['rows']:
    if r.get('rE') is None: continue
    L.append(f"| {hl(int(r['h']))} | {r['n']} | {s(r['rE'],0)} [{s(r['rE_q'][0],0)}, {s(r['rE_q'][1],0)}] | {s(r['rI'],0)} [{s(r['rI_q'][0],0)}, {s(r['rI_q'][1],0)}] | {s(r['fold'],2)} [{s(r['fold_q'][0],2)}, {s(r['fold_q'][1],2)}] | {r['within1p5']} of {r['n']} | {('<0.001' if r['p']<0.001 else s(r['p'],3))} |")
tabs['t6']='\n'.join(L)
# T7
L=['| Reference for E | Reference for I | Cells with a crossing | Matching command holding (mV), median [IQR] | Range (mV) |','|---|---|---|---|---|']
for k,(re_,ri_) in {'-80/0':('−80','0'),'-70/0':('−70','0'),'-80/10':('−80','+10'),'-70/10':('−70','+10')}.items():
    c=np.array(Q['ei'][k]['cross']); q=np.percentile(c,[25,50,75])
    L.append(f"| {re_} | {ri_} | {len(c)} of {Q['ei'][k]['n_ref']} | {s(q[1],1)} [{s(q[0],1)}, {s(q[2],1)}] | {s(c.min(),0)} to {s(c.max(),0)} |")
tabs['t7']='\n'.join(L)
for k,v in tabs.items(): print('\n##',k); print(v)
json.dump(tabs,open('cs_tables_md.json','w'),ensure_ascii=False)
# extras
cc=pd.read_csv('k60_capcheck_cs.csv'); cc['t']=pd.to_datetime(cc.time); cc['ck']=cc.date.astype(str)+'_s'+cc.slice.astype(str)+'_c'+cc.cell.astype(str)
c=cc[cc.Ra_MOhm.notna()].sort_values('t'); r=c.groupby('ck').agg(f=('Ra_MOhm','first'),l=('Ra_MOhm','last'),mx=('Ra_MOhm','max'),n=('Ra_MOhm','size'),rf=('Rm_MOhm','first'),rl=('Rm_MOhm','last'),sp=('t',lambda x:(x.iloc[-1]-x.iloc[0]).total_seconds()/60))
r=r[r.sp>=3]; ch=(r.l-r.f)/r.f*100
print('cells',len(c.ck.unique()),'with >=2 checks and span>=3 min',len(r),'dRa% med',ch.median().round(1),np.percentile(ch,[25,75]).round(1),'Ra increase >30%',(ch>30).sum(),'max Ra cells>=35',(r.mx>=35).sum(), 'all cells Ra>=35', c[c.Ra_MOhm>=35].ck.unique())
print('Ra median all checks',c.Ra_MOhm.median().round(1),'range',c.Ra_MOhm.min().round(1),c.Ra_MOhm.max().round(1),'n',len(c))
print('Rm<50:',c[c.Rm_MOhm<50].ck.unique(), 'Cm range',c.Cm_pF.min(),c.Cm_pF.max())
print('records with |Ihold|>1000',(a.Ihold_median_pA.abs()>1000).sum(), 'noise>8',(a.noise_lp3k_pA>8).sum(), a[a.noise_lp3k_pA>8].cellkey.unique())
print('recs with transients',(a.n_transients>0).sum(),len(a),'by cell top',a.groupby('cellkey').n_transients.sum().sort_values().tail(4).to_dict())
