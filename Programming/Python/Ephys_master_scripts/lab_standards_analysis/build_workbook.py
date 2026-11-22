import pandas as pd, json
T=json.load(open('lab_tables_md.json')); S=json.load(open('lab_stats.json'))
def md2df(md):
    L=[l for l in md.split('\n') if l.startswith('|')]
    rows=[[c.strip() for c in l.strip('|').split('|')] for l in L if not set(l.replace('|','').strip())<=set('- ')]
    return pd.DataFrame(rows[1:],columns=rows[0])
D=pd.read_csv('lab_cells_long.csv')
readme=pd.DataFrame({'README':[
 'Lab standard values: EPSC/IPSC frequency and E/I ratio. Past lab (Peerboom et al. 2023 controls) vs current lab (Helena, Sam, Advaith in Sam experiment) vs Advaith Cs-gluconate.',
 'cells_long: one row per cell (frequency Hz = miniML events per valid second; paper = digitised from figure; Sam = his processed values). freq_d2 = second detector (where available).',
 'Table_1..Table_6 match the report. paper_digitised_dots: all dots read from Fig 3C,G and 6E,I (CTL and VU; only CTL used).',
 'E/I = EPSC Hz / IPSC Hz. Group-level ratio uses separate cells for E and I (bootstrap 95% CI).',
 'Paper Fig 6E y-axis: ticks 2.0, 1.5, 1.0 and baseline labelled 0 one step below 1.0; values read from ticks. If baseline is 0, DIV21 EPSC values are 0.5 Hz lower.']})
with pd.ExcelWriter('Lab_standards_analysis.xlsx') as w:
    readme.to_excel(w,'README',index=False)
    for i,k in enumerate(['E','I','testsE','testsI','ratio','mid','d2'],1):
        pass
    md2df(T['E']).to_excel(w,'Table_1_EPSC_freq',index=False); md2df(T['I']).to_excel(w,'Table_2_IPSC_freq',index=False)
    md2df(T['testsE']).to_excel(w,'Table_3a_tests_EPSC',index=False); md2df(T['testsI']).to_excel(w,'Table_3b_tests_IPSC',index=False)
    md2df(T['ratio']).to_excel(w,'Table_4_EI_individual',index=False); md2df(T['mid']).to_excel(w,'Table_5_EI_intermediate_Cs',index=False); md2df(T['d2']).to_excel(w,'Table_6_second_detector',index=False)
    D.to_excel(w,'cells_long',index=False)
    pd.read_csv('lab_box_stats.csv').to_excel(w,'box_stats',index=False)
    pd.read_csv('cs_E_I_same_cell.csv').to_excel(w,'Cs_E_I_same_cell',index=False)
    pd.read_csv('/home/claude/work/cs/cs_ei_per_cell_ref-80_0.csv').to_excel(w,'Cs_EI_per_cell_holding',index=False)
    pd.read_csv('paper_dots.csv').to_excel(w,'paper_digitised_dots',index=False)
