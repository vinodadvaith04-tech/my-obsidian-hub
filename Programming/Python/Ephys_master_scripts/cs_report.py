"""
cs_report.py - the Cs-gluconate pilot tables and figures, built from Ephys_master.xlsx.

Reads  : sheets miniML_per_file, alt_det_per_file, manual_qc of the master (nothing else)
Writes : <out>/Cs_pilot_tables.xlsx   ONE workbook, one sheet per table (README, recordings_long, per_recording_both_dir,
                                      per_cell_holding, pooled_by_holding, EI_ratio_sep_vs_int, vs_earlier_methods, qc_manual, floors)
         <out>/figs/fig1..fig7 .png
Usage  : python scripts/cs_report.py [--master Ephys_master.xlsx] [--out <Ephys folder>/Cs_pilot_package] [--params-hash HASH]
Direction semantics: negative = inward, positive = outward. All Cs numbers use COMMAND holdings (no LJP correction).
"""
import argparse, json, os, sys
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).resolve().parent))
import ephys_master as em
ap = argparse.ArgumentParser()
ap.add_argument("--master", type=Path, default=em.default_master())
ap.add_argument("--out", type=Path, default=None)
ap.add_argument("--params-hash", default=None)
A = ap.parse_args()
OUT = str(A.out or (A.master.parent / "Cs_pilot_package")); os.makedirs(OUT + "/figs", exist_ok=True)
T = {}   # sheet name -> DataFrame
plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,"axes.titlesize":10,"axes.labelsize":10,"figure.dpi":110,"savefig.dpi":170})
RED="#d6001c"; BLUE="#1f5fbf"; GREY="#6b6b6b"; LG="#c8c8c8"
_runs = em.read_sheet(A.master, "runs")
HASH = A.params_hash or _runs[_runs.tool == "miniML"].sort_values("last_processed").params_hash.iloc[-1]
SF = em.read_sheet(A.master, "miniML_per_file"); SF = SF[SF.params_hash == HASH].copy()
print("miniML run", HASH, "-", len(SF), "recording-directions")
lab={"20261104":"4 Nov","20261105":"5 Nov","20261106":"6 Nov","20261112":"12 Nov"}
cs=SF[SF.internal.astype(str).str.contains("Cs") & SF.new_name.str[:8].isin(lab)].drop_duplicates(["new_name","direction"]).copy()
cs["day"]=cs.new_name.str[:8]
cs["cellid"]=cs.day.map(lab)+" C"+cs.cell.astype(int).astype(str)
cs["hold"]=cs.holding_command_mV
T["cs_recordings_long"]=cs
# wide per recording
w=cs.pivot_table(index=["cellid","day","new_name","hold","duration_s"],columns="direction",values=["n_events","frequency_Hz","amplitude_abs_median","noise_mad_pA","hold_current_median_pA"],aggfunc="first")
w.columns=[f"{a}_{b}" for a,b in w.columns]; w=w.reset_index()
w["down"]=w.frequency_Hz_negative; w["up"]=w.frequency_Hz_positive
w["noise"]=w[["noise_mad_pA_negative","noise_mad_pA_positive"]].mean(axis=1); w["holdI"]=w[["hold_current_median_pA_negative","hold_current_median_pA_positive"]].mean(axis=1)
w["upfrac"]=w.up/(w.up+w.down)
order=["4 Nov C1","5 Nov C1","6 Nov C1"]+[f"12 Nov C{i}" for i in range(1,7)]
w["cellid"]=pd.Categorical(w.cellid,order,ordered=True); w=w.sort_values(["cellid","new_name"])
T["cs_per_recording_both_directions"]=w
# old detector
O=[em.read_sheet(A.master,"alt_det_per_file")]
O=pd.concat(O).drop_duplicates("new_name")[["new_name","freq_inward_Hz","freq_outward_Hz"]]
w=w.merge(O,on="new_name",how="left")
g=w.groupby(["cellid","hold"],observed=True).agg(n_rec=("new_name","count"),dur_s=("duration_s","sum"),down=("down","mean"),up=("up","mean"),old_down=("freq_inward_Hz","mean"),old_up=("freq_outward_Hz","mean"),
    amp_down=("amplitude_abs_median_negative","median"),amp_up=("amplitude_abs_median_positive","median"),noise=("noise","mean"),holdI=("holdI","median")).reset_index()
g["upfrac"]=g.up/(g.up+g.down); T["cs_per_cell_holding"]=g
# floors (wrong-direction rate)
floor_down=g[g.hold>=0].down.dropna(); floor_up=g[g.hold<=-70].up.dropna()
FL=dict(down_median=float(floor_down.median()),down_p90=float(floor_down.quantile(.9)),up_median=float(floor_up.median()),up_p90=float(floor_up.quantile(.9)),n_down=len(floor_down),n_up=len(floor_up))
print("floors",FL)
# ---------- Fig1 per cell
fig,axs=plt.subplots(3,3,figsize=(11,8.2),sharey=False)
for ax,c in zip(axs.ravel(),order):
    d=g[g.cellid==c].sort_values("hold")
    ax.axhline(FL["down_p90"],color=RED,lw=.6,ls=":"); ax.axhline(FL["up_p90"],color=BLUE,lw=.6,ls=":")
    ax.plot(d.hold,d.down,"o-",color=RED,ms=4,lw=1,label="inward (miniML)"); ax.plot(d.hold,d.up,"s-",color=BLUE,ms=4,lw=1,label="outward (miniML)")
    ax.plot(d.hold,d.old_down,"o",mfc="none",color=RED,ms=7,mew=.8); ax.plot(d.hold,d.old_up,"s",mfc="none",color=BLUE,ms=7,mew=.8)
    ax.set_title(c); ax.set_xlim(-85,20); ax.set_xlabel("command holding (mV)"); ax.set_ylabel("events per second")
    ax.set_ylim(bottom=0)
axs[0,0].legend(fontsize=7,frameon=False,loc="upper left",title="open = second detector",title_fontsize=7)
fig.suptitle("Inward and outward events per cell and holding (dotted lines: wrong-direction rate, 90th percentile)",fontsize=10,y=.995)
fig.tight_layout(); fig.savefig(OUT+"/figs/fig1_per_cell.png"); plt.close(fig)
# ---------- Fig2 pooled by holding
fig,axs=plt.subplots(1,3,figsize=(11.5,3.6))
hs=sorted(g.hold.unique())
for ax,(col,colr,ttl) in zip(axs[:2],[("down",RED,"Inward events"),("up",BLUE,"Outward events")]):
    for c in order:
        d=g[g.cellid==c].sort_values("hold"); ax.plot(d.hold,d[col],"-",color=LG,lw=.8)
    m=g.groupby("hold")[col].median(); n=g.groupby("hold")[col].count()
    ax.plot(m.index,m.values,"o-",color=colr,lw=2,ms=5); ax.set_title(ttl+" (median across cells)"); ax.set_xlabel("command holding (mV)"); ax.set_ylabel("events per second")
    for h_,v in m.items(): ax.annotate(str(n[h_]),(h_,v),textcoords="offset points",xytext=(0,6),ha="center",fontsize=7,color=GREY)
ax=axs[2]
for c in order:
    d=g[g.cellid==c].sort_values("hold"); ax.plot(d.hold,d.upfrac,"-",color=LG,lw=.8)
m=g.groupby("hold").upfrac.median(); ax.plot(m.index,m.values,"o-",color="k",lw=2,ms=5); ax.set_ylim(0,1.02); ax.axhline(.5,color=GREY,lw=.5,ls=":")
ax.set_title("Outward share of events"); ax.set_xlabel("command holding (mV)"); ax.set_ylabel("up / (up + down)")
fig.text(.01,.01,"Numbers above points: cells at that holding. Grey lines: individual cells.",fontsize=7,color=GREY)
fig.tight_layout(rect=(0,.03,1,1)); fig.savefig(OUT+"/figs/fig2_pooled.png"); plt.close(fig)
# ---------- QC
q=em.read_sheet(A.master,"manual_qc"); q=q[q.session.astype(str).str.contains("Cs")].copy(); T["qc_manual"]=q
qlab={"2026-11-04":"4 Nov","2026-11-05":"5 Nov","2026-11-06":"6 Nov","2026-11-12":"12 Nov"}
q["cellid"]=q.day.astype(str).map(qlab)+" C"+q.cell.astype(int).astype(str)
fig,axs=plt.subplots(1,4,figsize=(12.5,3.3))
for c,d in q.groupby("cellid"):
    for ax,col in zip(axs[:3],["Ra","Rm","Cm"]):
        dd=d.dropna(subset=[col]); 
        if len(dd)>1: ax.plot(dd.check,dd[col],"o-",lw=1,ms=3,label=c)
        elif len(dd)==1: ax.plot(dd.check,dd[col],"o",ms=4)
axs[0].axhline(35,color=RED,lw=.8,ls="--"); axs[0].set_title("Access resistance (MΩ)"); axs[1].set_title("Membrane resistance (MΩ)"); axs[2].set_title("Membrane capacitance (pF)")
for ax in axs[:3]: ax.set_xlabel("manual check number")
axs[0].text(.2,35.8,"exclusion limit 35 MΩ",color=RED,fontsize=7)
axs[0].legend(fontsize=6,frameon=False,ncol=2)
ax=axs[3]
for c in order:
    d=g[g.cellid==c].sort_values("hold"); ax.plot(d.hold,d.holdI,"o-",lw=.9,ms=3)
ax.set_title("Holding current (pA, median)"); ax.set_xlabel("command holding (mV)"); ax.axhline(0,color=GREY,lw=.5)
fig.tight_layout(); fig.savefig(OUT+"/figs/fig3_qc.png"); plt.close(fig)
# ---------- noise vs holding Cs vs KGlu
kg=SF[(SF.internal=="KGlu")&(SF.direction=="negative")&(SF.drug.isna())]
fig,ax=plt.subplots(figsize=(5.5,3.4))
for c in order:
    d=g[g.cellid==c].sort_values("hold"); ax.plot(d.hold,d.noise,"-",color=LG,lw=.8)
m=g.groupby("hold").noise.median(); ax.plot(m.index,m.values,"o-",color="k",label="Cs-gluconate median")
mk=kg.groupby("holding_command_mV").noise_mad_pA.median(); ax.plot(mk.index,mk.values,"s--",color=GREY,label="K-gluconate median (no blocker)")
ax.set_xlabel("command holding (mV)"); ax.set_ylabel("baseline noise, MAD (pA)"); ax.set_yscale("log"); ax.legend(fontsize=7,frameon=False)
fig.tight_layout(); fig.savefig(OUT+"/figs/fig4_noise.png"); plt.close(fig)
# ---------- comparison with earlier methods
ee=SF[(SF.session=="Sam EPSC (22 Oct)")&(SF.direction=="negative")]
ii=SF[(SF.session=="Sam IPSC")&(SF.direction=="negative")]
kg60=SF[(SF.internal=="KGlu")&(SF.drug.isna())&(SF.direction=="negative")&(SF.holding_command_mV==-60)]
def cellmean(d,col): return d.groupby("recording_name")[col].mean()
csE=g[g.hold.isin([-80,-70])].groupby("cellid",observed=True).agg(f=("down","mean"),a=("amp_down","median"))
csI=g[g.hold.isin([10,15])].groupby("cellid",observed=True).agg(f=("up","mean"),a=("amp_up","median"))
cmp={"EPSC":{"gabazine, K-gluc, −60 mV (n cells)":(cellmean(ee,"frequency_Hz"),cellmean(ee,"amplitude_abs_median")),"Cs-gluc, −80/−70 mV, no blocker":(csE.f,csE.a),"K-gluc, −60 mV, no blocker":(cellmean(kg60,"frequency_Hz"),cellmean(kg60,"amplitude_abs_median"))},
     "IPSC":{"KCl, DNQX+AP5, −60 mV":(cellmean(ii,"frequency_Hz"),cellmean(ii,"amplitude_abs_median")),"Cs-gluc, +10/+15 mV, no blocker":(csI.f,csI.a)}}
fig,axs=plt.subplots(2,2,figsize=(10,6.4))
for r,(k,dd) in enumerate(cmp.items()):
    for c,(idx,ttl) in enumerate([(0,"frequency (events per second)"),(1,"median amplitude (pA)")]):
        ax=axs[r,c]; names=list(dd); 
        for i,nm in enumerate(names):
            v=dd[nm][idx].dropna().values
            ax.scatter(np.random.default_rng(1).normal(i,.06,len(v)),v,s=18,color=RED if k=="EPSC" else BLUE,alpha=.7,edgecolor="none")
            if len(v): ax.hlines(np.median(v),i-.25,i+.25,color="k",lw=2)
            ax.text(i,ax.get_ylim()[1]*0.0,"",ha="center")
        ax.set_xticks(range(len(names))); ax.set_xticklabels([f"{n}\n(n={len(dd[n][idx].dropna())})" for n in names],fontsize=7)
        ax.set_title(f"{k}: {ttl}"); ax.set_ylim(bottom=0)
fig.tight_layout(); fig.savefig(OUT+"/figs/fig5_vs_earlier_methods.png"); plt.close(fig)
cmpt=[]
for k,dd in cmp.items():
    for nm,(f,a) in dd.items(): cmpt.append(dict(type=k,method=nm,n_cells=int(f.dropna().shape[0]),freq_median=float(np.nanmedian(f)),freq_iqr_lo=float(np.nanpercentile(f.dropna(),25)),freq_iqr_hi=float(np.nanpercentile(f.dropna(),75)),amp_median=float(np.nanmedian(a)),))
T["cs_vs_earlier_methods"]=pd.DataFrame(cmpt); print(pd.DataFrame(cmpt).round(2).to_string())
# ---------- upward fraction, Cs vs KGlu intermediates
kgi=SF[(SF.internal=="KGlu")&(SF.drug.isna())&(SF.holding_command_mV.isin([-60,-50,-40,-35]))]
kgw=kgi.pivot_table(index=["recording_name","holding_command_mV"],columns="direction",values="frequency_Hz",aggfunc="mean").reset_index()
kgw["upfrac"]=kgw.positive/(kgw.positive+kgw.negative)
fig,ax=plt.subplots(figsize=(5.8,3.6))
k=kgw.dropna(subset=["upfrac"]).groupby("holding_command_mV").upfrac.agg(["median","count"])
for c in order:
    d=g[g.cellid==c].sort_values("hold"); ax.plot(d.hold,d.upfrac,"-",color=LG,lw=.7)
mm=g[g.hold.isin([-60,-50,-40,-35,-30])].groupby("hold").upfrac.median()
ax.plot(mm.index,mm.values,"o-",color="k",label="Cs-gluconate (median)"); ax.plot(k.index,k["median"],"s--",color=GREY,label="K-gluconate (median)")
ax.set_xlim(-65,-25); ax.set_ylim(0,1); ax.set_xlabel("command holding (mV)"); ax.set_ylabel("outward share, up/(up+down)"); ax.legend(fontsize=7,frameon=False)
fig.tight_layout(); fig.savefig(OUT+"/figs/fig6_upfrac_cs_vs_kglu.png"); plt.close(fig)
print(k)
# ---------- E/I ratios within cell
sepE=lambda d: d[d.hold.isin([-80,-70])].down.mean(); sepI=lambda d: d[d.hold.isin([10,15])].up.mean()
rows=[]
for c in order:
    d=g[g.cellid==c]; e=sepE(d); i=sepI(d)
    for h in [-50,-40,-30,-20]:
        x=d[d.hold==h]
        if len(x): rows.append(dict(cellid=c,hold=h,E_sep=e,I_sep=i,ratio_sep=(e/i if i and i>0 else np.nan),down=float(x.down.iloc[0]),up=float(x.up.iloc[0]),ratio_int=float(x.down.iloc[0]/x.up.iloc[0]) if x.up.iloc[0]>0 else np.nan))
R=pd.DataFrame(rows); T["cs_EI_ratio_separate_vs_intermediate"]=R; print(R.round(2).to_string())
fig,ax=plt.subplots(figsize=(6.2,3.8))
for h,cl in zip([-50,-40,-30,-20],["#999","#444",BLUE,"#00a0a0"]):
    r=R[R.hold==h].dropna(subset=["ratio_sep","ratio_int"])
    ax.scatter(r.ratio_sep,r.ratio_int,s=35,color=cl,label=f"intermediate {h} mV",edgecolor="none")
lim=[0.01,100]; ax.plot(lim,lim,color=GREY,lw=.6,ls=":"); ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("E/I from separate holdings\n(inward at −80/−70, outward at +10/+15)"); ax.set_ylabel("inward/outward at the intermediate holding"); ax.legend(fontsize=7,frameon=False)
fig.tight_layout(); fig.savefig(OUT+"/figs/fig7_EI_ratio.png"); plt.close(fig)
# ---------- summary tables
pool=g.groupby("hold").agg(cells=("cellid","nunique"),recordings=("n_rec","sum"),down_med=("down","median"),up_med=("up","median"),upfrac_med=("upfrac","median"),amp_down=("amp_down","median"),amp_up=("amp_up","median"),noise_med=("noise","median")).reset_index()
T["cs_pooled_by_holding"]=pool; print(pool.round(2).to_string())
T["floors"]=pd.DataFrame([FL])

# ---------- one workbook
from openpyxl.styles import Font, PatternFill
from openpyxl.utils import get_column_letter
README = pd.DataFrame({"sheet": ["recordings_long","per_recording_both_dir","per_cell_holding","pooled_by_holding","EI_ratio_sep_vs_int","vs_earlier_methods","qc_manual","floors"],
    "what it is": ["one row per Cs recording x direction (miniML), straight from the master",
                   "one row per Cs recording, inward and outward side by side, plus the second detector frequencies",
                   "one row per cell x command holding (means over recordings): the main table",
                   "per holding, across cells (medians)",
                   "E/I ratio from separate holdings (inward at -80/-70, outward at +10/+15) vs inward/outward at one intermediate holding",
                   "Cs vs earlier methods (gabazine EPSC, KCl IPSC, K-gluconate)",
                   "manual patch checks (pipette, Cm, Rm, Ra) from the lab journal",
                   "wrong-direction rates (inward at >= 0 mV, outward at <= -70 mV): the false-detection floor"]})
sheets = {"README": README}
NAMES = {"cs_recordings_long": "recordings_long", "cs_per_recording_both_directions": "per_recording_both_dir",
         "cs_per_cell_holding": "per_cell_holding", "cs_pooled_by_holding": "pooled_by_holding",
         "cs_EI_ratio_separate_vs_intermediate": "EI_ratio_sep_vs_int", "cs_vs_earlier_methods": "vs_earlier_methods",
         "qc_manual": "qc_manual", "floors": "floors"}
sheets.update({NAMES[k]: v for k, v in T.items()})
xl = Path(OUT) / "Cs_pilot_tables.xlsx"
with pd.ExcelWriter(xl, engine="openpyxl") as xw:
    for name, df in sheets.items():
        name = name[:31]; df.to_excel(xw, sheet_name=name, index=False); ws = xw.sheets[name]
        for c in ws[1]: c.font = Font(bold=True); c.fill = PatternFill("solid", fgColor="D9D9D9")
        ws.freeze_panes = "A2"
        for j, col in enumerate(df.columns, 1):
            w = max([len(str(col))] + [len(f"{v:.4g}") if isinstance(v, float) else len(str(v)) for v in df[col].head(200).tolist()])
            ws.column_dimensions[get_column_letter(j)].width = min(max(w + 2, 8), 60)
print("written", xl)
