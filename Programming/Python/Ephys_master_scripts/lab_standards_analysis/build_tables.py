import json, pandas as pd, numpy as np
S=json.load(open('lab_stats.json')); D=pd.read_csv('lab_cells_long.csv')
def f2(x): return f'{x:.2f}'
def row(d): return f"{d['median']:.2f} [{d['q1']:.2f}–{d['q3']:.2f}]",f"{d['mean']:.2f} ± {d['sem']:.2f}",f"{d['min']:.2f}–{d['max']:.2f}"
T={}
def freq_table(meas,items):
    L=['| Group | Cells | Median [IQR] (Hz) | Mean ± SEM (Hz) | Range (Hz) |','| --- | --- | --- | --- | --- |']
    for name,key in items:
        d=S['desc'][key]; a,b,c=row(d); L.append(f"| {name} | {d['n']} | {a} | {b} | {c} |")
    return '\n'.join(L)
T['E']=freq_table('E',[('Past lab, paper control DIV9','E|Paper DIV9'),('Past lab, paper control DIV21','E|Paper DIV21'),
 ("Current lab: Sam's processed data (K-Glu, gabazine, −60 mV)","E|src|Sam's processed EPSC, control cells (his pipeline)"),
 ("Current lab: Advaith in Sam's experiment (K-Glu, gabazine, −60 mV)",'E|src|Advaith EPSC, Sam experiment (miniML)'),
 ('Current lab, pooled','E|Current lab (pooled)'),('Advaith Cs-gluconate, −80 mV','E|Cs-gluconate (-80 mV)')])
T['I']=freq_table('I',[('Past lab, paper control DIV9','I|Paper DIV9'),('Past lab, paper control DIV21','I|Paper DIV21'),
 ('Current lab: Helena (KCl, DNQX+AP5, −60 mV)','I|src|Helena KCl IPSC (miniML)'),('Current lab: Advaith in Sam\'s experiment (KCl, DNQX+AP5, −60 mV)','I|src|Advaith KCl IPSC, Sam experiment (miniML)'),
 ('Current lab, pooled','I|Current lab (pooled)'),('Advaith Cs-gluconate, 0 mV','I|Cs-gluconate (0 mV)')])
def pf(p): return '<0.001' if p<0.001 else f'{p:.3f}'
for m in 'EI':
    L=['| Comparison | Fold (median ratio) | p (Mann–Whitney) | p (Holm) |','| --- | --- | --- | --- |']
    for r in S['tests'][m]: L.append(f"| {r['a']} vs {r['b']} | {r['fold']:.2f} | {pf(r['p'])} | {pf(r['p_holm'])} |")
    T['tests'+m]='\n'.join(L)
L=['| Group | Cells (E / I) | E/I ratio of means [95% CI] | E/I ratio of medians [95% CI] |','| --- | --- | --- | --- |']
for k,r in S['ratio'].items(): L.append(f"| {k} | {r['nE']} / {r['nI']} | {r['ratio_mean']:.2f} [{r['ci_mean'][0]:.2f}–{r['ci_mean'][1]:.2f}] | {r['ratio_med']:.2f} [{r['ci_med'][0]:.2f}–{r['ci_med'][1]:.2f}] |")
c=S['cs_within']; L.append(f"| Advaith Cs-gluconate, same cell (−80 mV / 0 mV) | {c['n']} cells | mean {c['mean']:.2f} | median {c['median']:.2f} [IQR {c['q1']:.2f}–{c['q3']:.2f}] |")
T['ratio']='\n'.join(L)
L=['| Command holding | Cells | E/I at the holding, median [IQR] | Separately measured E/I (same cells), median [IQR] | Fold, median [IQR] | Cells within ×1.5 | p (Wilcoxon) |','| --- | --- | --- | --- | --- | --- | --- |']
for h,r in S['cs_mid'].items():
    L.append(f"| {str(h).replace('-','−')} mV | {r['n']} | {r['ratio_h_med']:.2f} [{r['ratio_h_q1']:.2f}–{r['ratio_h_q3']:.2f}] | {r['ratio_ref_med']:.2f} [{r['ratio_ref_q1']:.2f}–{r['ratio_ref_q3']:.2f}] | {r['fold_med']:.2f} [{r['fold_q1']:.2f}–{r['fold_q3']:.2f}] | {r['n_within_1p5']} of {r['n']} | {pf(r['p'])} |")
T['mid']='\n'.join(L)
L=['| Group | Cells | Median [IQR], second detector (Hz) | Median, miniML (Hz) |','| --- | --- | --- | --- |']
for k,r in S['d2'].items():
    meas,src=k.split('|'); mm=S['desc'][f'{meas}|src|{src}']['median']
    L.append(f"| {meas}: {src.replace(' (miniML)','')} | {r['n']} | {r['median']:.2f} [{r['q1']:.2f}–{r['q3']:.2f}] | {mm:.2f} |")
T['d2']='\n'.join(L)
json.dump(T,open('lab_tables_md.json','w'),indent=1,ensure_ascii=False)
for k,v in T.items(): print('\n##',k); print(v)
