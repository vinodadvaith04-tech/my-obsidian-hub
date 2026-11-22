import pandas as pd, numpy as np, json, itertools
from scipy import stats
pd.set_option('display.width',250,'display.max_columns',80,'display.max_rows',400)
def run(csv,tag):
    a=pd.read_csv(csv)
    cm=lambda c: a.groupby(['cellkey','holding_command_mV'])[c].mean().unstack()
    fi,fo,qi,qo,ai,ao=cm('f_in'),cm('f_out'),cm('q_in'),cm('q_out'),cm('amp_in'),cm('amp_out')
    d2i,d2o=cm('d2_freq_inward_Hz'),cm('d2_freq_outward_Hz')
    res={}
    # ---- purity
    pur=[]
    for h in sorted(fi.columns):
        t=pd.DataFrame({'fi':fi[h],'fo':fo[h],'d2i':d2i[h],'d2o':d2o[h]}).dropna(subset=['fi','fo'])
        tot=t.fi+t.fo; pe=(t.fi/tot).dropna()*100; pi=(t.fo/tot).dropna()*100
        opp_in_pct=(t.fo/t.fi.replace(0,np.nan)*100).dropna(); opp_out_pct=(t.fi/t.fo.replace(0,np.nan)*100).dropna()
        pur.append(dict(h=h,n=len(t),pctE=np.median(pe),pctE_q=np.percentile(pe,[25,75]),pctI=np.median(pi),pctI_q=np.percentile(pi,[25,75]),
            fo_over_fi=np.median(opp_in_pct) if len(opp_in_pct) else np.nan, fi_over_fo=np.median(opp_out_pct) if len(opp_out_pct) else np.nan,
            fo_pooled_pct=t.fo.sum()/t.fi.sum()*100, fi_pooled_pct=t.fi.sum()/t.fo.sum()*100,
            d2_pctE=np.nanmedian((t.d2i/(t.d2i+t.d2o))*100), d2_pctI=np.nanmedian((t.d2o/(t.d2i+t.d2o))*100)))
    res['purity']=pur
    print(tag,'PURITY'); print(pd.DataFrame(pur).round(1).to_string())
    # paired comparisons of isolating candidates
    def w(x,y,lab):
        d=pd.DataFrame({'x':x,'y':y}).dropna()
        if len(d)<5: return (lab,len(d),np.nan)
        try: p=stats.wilcoxon(d.x,d.y).pvalue
        except Exception: p=np.nan
        return (lab,len(d),round(p,4),round(d.x.median(),3),round(d.y.median(),3))
    cmp=[w(fo[-80],fo[-70],'f_out -80 vs -70'),w(fi[-80],fi[-70],'f_in -80 vs -70'),w(fi[0],fi[10],'f_in 0 vs +10'),w(fo[0],fo[10],'f_out 0 vs +10'),
         w(fo[-80],fi[-80],'-80: f_out vs f_in'),w(fi[0],fo[0],'0: f_in vs f_out'),w(fi[10],fo[10],'+10: f_in vs f_out'),w(fo[-60],fo[-80],'f_out -60 vs -80'),w(fo[-50],fo[-80],'f_out -50 vs -80')]
    res['cmp']=cmp; print(cmp)
    # ---- E/I preservation
    mid=[-60,-50,-40,-30]
    out={}
    for (e,i) in [(-80,0),(-70,0),(-80,10),(-70,10)]:
        E0,I0=fi[e],fo[i]; QE,QI=(fi[e]*qi[e]),(fo[i]*qo[i])
        base=pd.DataFrame({'E0':E0,'I0':I0}).dropna(); base=base[(base.E0>0)&(base.I0>0)]
        rows=[]
        for h in mid:
            t=pd.DataFrame({'E0':E0,'I0':I0,'fi':fi[h],'fo':fo[h],'Q_E0':QE,'Q_I0':QI,'Qi':fi[h]*qi[h],'Qo':fo[h]*qo[h]}).dropna()
            t=t[(t.E0>0)&(t.I0>0)&(t.fi>0)&(t.fo>0)]
            if len(t)<3: rows.append(dict(h=h,n=len(t))); continue
            rE=t.fi/t.E0; rI=t.fo/t.I0; fold=rE/rI; lf=np.log2(fold)
            qf=(t.Qi/t.Qo)/(t.Q_E0/t.Q_I0)
            try: p=stats.wilcoxon(lf).pvalue
            except Exception: p=np.nan
            rows.append(dict(h=h,n=len(t),ratio_ref=np.median(t.E0/t.I0),ratio_h=np.median(t.fi/t.fo),rE=np.median(rE)*100,rE_q=np.percentile(rE,[25,75])*100,rI=np.median(rI)*100,rI_q=np.percentile(rI,[25,75])*100,
                fold=np.median(fold),fold_q=np.percentile(fold,[25,75]),log2fold_mean=lf.mean(),p=p,within1p5=int(((fold>1/1.5)&(fold<1.5)).sum()),
                qfold=np.median(qf),qfold_q=np.percentile(qf,[25,75])))
        # interpolation: holding where ratio_h equals reference, per cell (linear in log2 ratio vs holding, -60..0)
        hs=[-60,-50,-40,-30,0]; cross=[]
        for c in base.index:
            xs,ys=[],[]
            for h in hs:
                fih,foh=fi.at[c,h] if h in fi.columns else np.nan, fo.at[c,h] if h in fo.columns else np.nan
                if pd.notna(fih) and pd.notna(foh) and fih>0 and foh>0: xs.append(h); ys.append(np.log2(fih/foh)-np.log2(base.at[c,'E0']/base.at[c,'I0']))
            for k in range(len(xs)-1):
                if ys[k]>=0>=ys[k+1] and ys[k]!=ys[k+1]:
                    cross.append(xs[k]+(0-ys[k])/(ys[k+1]-ys[k])*(xs[k+1]-xs[k])); break
        out[f'{e}/{i}']=dict(rows=rows,n_ref=len(base),cross=cross)
        print(tag,'REF E',e,'I',i,'n_ref',len(base),'crossing',np.round(cross,1),'median',np.median(cross) if cross else None)
        print(pd.DataFrame(rows).round(2).to_string())
    res['ei']=out
    return res
r1=run('cs_analysis_set.csv','PRIMARY'); r2=run('cs_analysis_set_strict.csv','STRICT')
json.dump(dict(primary=r1,strict=r2),open('cs_tables_part2.json','w'),default=lambda o:float(o) if isinstance(o,(np.floating,np.integer)) else (o.tolist() if hasattr(o,'tolist') else str(o)))
