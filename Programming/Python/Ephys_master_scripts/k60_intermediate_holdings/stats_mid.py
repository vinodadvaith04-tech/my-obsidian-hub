import pandas as pd, numpy as np, json
from scipy import stats
g=pd.read_csv('mid_master_table.csv')
ch=g.groupby(['cellkey','holding_command_mV']).agg(f_in=('f_in','mean'),f_out=('f_out','mean'),noise=('noise_lp3k_pA','mean'),abs_drift=('abs_drift','mean'),ihold=('Ihold_median_pA','mean'),
   amp_in=('amp_in','mean'),amp_out=('amp_out','mean'),n_in=('n_in','sum'),n_out=('n_out','sum'),dur=('duration_s','sum'),Ra=('Ra_b','mean')).reset_index()
ch.to_csv('mid_cell_holding.csv',index=False)
print('cell-holding obs',len(ch))
res=[]
for m in ['f_in','f_out','noise','abs_drift','ihold','amp_in','amp_out','Ra']:
    x=ch.dropna(subset=[m]); r,p=stats.spearmanr(x.holding_command_mV,x[m]); res.append((m,len(x),round(r,2),round(p,4)))
print(pd.DataFrame(res,columns=['metric','n','rho(holding)','p']))
for a,b in [(-50,-40),(-50,-30),(-40,-30)]:
    w=ch[ch.holding_command_mV==a].set_index('cellkey').join(ch[ch.holding_command_mV==b].set_index('cellkey'),lsuffix='_a',rsuffix='_b',how='inner')
    print('\npaired',a,'vs',b,'cells',len(w))
    for m in ['f_in','f_out','noise','abs_drift','ihold','amp_out']:
        x=w[[m+'_a',m+'_b']].dropna()
        if len(x)>=4:
            try: p=stats.wilcoxon(x[m+'_a'],x[m+'_b']).pvalue
            except Exception: p=np.nan
            print(f'  {m}: n={len(x)} mean_a={x[m+"_a"].mean():.2f} mean_b={x[m+"_b"].mean():.2f} p={p:.3f}')
# cell dominance
for d in ['n_in','n_out']:
    t=ch.groupby('cellkey')[d].sum().sort_values(ascending=False); print(d,'total',t.sum(),'top3',t.head(3).to_dict(),'top1 share',round(t.iloc[0]/t.sum(),2))
for h in [-50,-45,-40,-35,-30]:
    x=ch[ch.holding_command_mV==h]
    print(h,'cells',len(x),'E mean',round(x.f_in.mean(),3),'med',round(x.f_in.median(),3),'| I mean',round(x.f_out.mean(),3),'med',round(x.f_out.median(),3),
      '| cells with both >=0.1Hz',int(((x.f_in>=.1)&(x.f_out>=.1)).sum()), '| both>=5 events',int(((x.n_in>=5)&(x.n_out>=5)).sum()))
# recording level both>=5
g['both5']=(g.n_in>=5)&(g.n_out>=5)
print(g.groupby('holding_command_mV').both5.agg(['sum','size']))
print(g.both5.sum(),len(g))
# E/I ratio per holding (pooled)
for h,d in g.groupby('holding_command_mV'): print(h,'pooled in/out freq',round(d.n_in.sum()/d.duration_s.sum(),3),round(d.n_out.sum()/d.duration_s.sum(),3),'out share',round(d.n_out.sum()/(d.n_in.sum()+d.n_out.sum()),2))
# ihold plausibility: flags
print(g[g.abs_ihold>200][['new_name','holding_command_mV','Ihold_median_pA','drift_total_pA','noise_lp3k_pA','n_in','n_out']].round(1).to_string())
print(g[(g.noise_lp3k_pA>6)|(g.n_transients>3)][['new_name','holding_command_mV','noise_lp3k_pA','n_transients','max_transient_pA']].round(1).to_string())
print(g[g.Ra_b>35][['new_name','holding_command_mV','Ra_b']].round(1).to_string())
print(g.groupby('date').agg(n=('new_name','size')).T)
