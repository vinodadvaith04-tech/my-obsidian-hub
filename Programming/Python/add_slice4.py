import pyabf, glob, re, shutil, os, numpy as np, pandas as pd, openpyxl
SRC="/mnt/user-data/uploads/Ephys pilot & practice/Ephys_master_table_DRAFT.xlsx"
OUT="/mnt/user-data/outputs/Ephys_master_table_DRAFT_plus_slice4.xlsx"
S4="/tmp/s4/Slice 4"; DATA="/tmp/ephys_data"
# holdings by cell and counter, from the 12-11 lab note (order of notes) matched to the file sequence (gap-free pairs)
H={1:{2:0,3:0,5:10},
   2:{7:-80,8:-80,10:-40,11:-40,13:-30,14:-30,16:0,17:0,18:10,19:10,22:-50,23:-50},
   3:{27:-80,28:-80,30:-50,31:-50,34:-40,35:-40,37:-70,38:-70,41:-30,42:-30,44:0,45:0,47:10,48:10},
   4:{51:-70,52:-70,53:-70},
   5:{59:-70,60:-70,62:-80},
   6:{65:-70,66:-70,69:-50,70:-50,73:-30,74:-30,76:10,77:10,81:-20,82:-20,85:-30,86:-30}}
FLAG={52:"short file (14.7 s), cell lost afterwards",56:"holding not in note (aborted attempt before the first -70 in the note), leaky: holding current about -1.1 nA",
      62:"-80 recording of cell 5 is only 39.5 s and holding current about -1.8 nA (leaky); mapping to the note's -80 by order",
      51:"cell lost afterwards (note)",53:"cell lost afterwards (note); 78.7 s"}
def dirs(h,internal="Cs-gluconate"):
    if h>=0: return "positive"
    if -50<=h<=-35: return "negative,positive"
    if h<=-60: return "negative"
    if h==-30: return "negative,positive"
    return "negative,positive"          # -20 (and any other unlisted value): both, nothing lost
rows=[]
for f in sorted(glob.glob(S4+"/*/*.abf")):
    cell=int(f.split("/")[-2].split()[-1]); ctr=int(re.search(r"_(\d+)\.abf",f).group(1)); old=os.path.basename(f)
    a=pyabf.ABF(f); proto=a.protocol or ""
    gap="GapFree" in proto
    ptype="sPSC gap-free" if gap else ("cap check" if "Cap" in proto else ("spike frequency" if "spike" in proto else "other"))
    new=f"20261112_Adv00_NT_unk_slice4_C{cell}_{ctr:03d}.abf"
    h=H.get(cell,{}).get(ctr,np.nan)
    use="Y" if (gap and a.dataLengthSec>=30 and not np.isnan(h)) else "N"
    flag=""
    if gap:
        flag=("from 12-11 note; matched by order to the file sequence" if not np.isnan(h) else "holding unknown")
        if ctr in FLAG: flag+="; "+FLAG[ctr]
    rows.append({"new_name":new,"old_name":old,"old_folder":"Sacha/Sacha recordings advaith/Caesium pilot/Slice 4/Cell %d"%cell,"new_folder":"Sacha/Sacha recordings advaith/Caesium pilot/Slice 4/Cell %d"%cell,
      "status":"renamed (name only in master; file not renamed on the Mac)","date_header":20261112,"start_time":a.abfDateTime.strftime("%H:%M:%S"),"session":"Cs internal pilot","S.No":"Adv00","T.Code":"NT","sex":"unk","sex_source":"not stated",
      "slice":4,"cell":float(cell),"counter":ctr,"protocol":proto,"protocol_type":ptype,"mode":"VC" if a.adcUnits[0]=="pA" else "CC","duration_s":round(a.dataLengthSec,1),"sweeps":a.sweepCount,
      "psc_content":"mixed E+I" if gap else np.nan,"drug":np.nan,"internal":"Cs-gluconate","internal_source":"stated","holding_header_mV":0,"holding_listed_in_note":h if gap else np.nan,"holding_command_mV":h if gap else np.nan,
      "holding_flag":flag if gap else np.nan,"use_for_miniML":use,"blockers":np.nan,"miniML_directions":dirs(h) if use=="Y" else np.nan})
    shutil.copy(f, os.path.join(DATA,new))
N=pd.DataFrame(rows); N.to_csv("/tmp/s4_rows.csv",index=False)
wb=openpyxl.load_workbook(SRC); ws=wb["files"]; hdr=[c.value for c in ws[1]]
for r in rows: ws.append([None if (isinstance(r.get(h),float) and np.isnan(r.get(h))) else r.get(h) for h in hdr])
wb.save(OUT)
print(N.groupby(["cell","use_for_miniML"]).size().unstack(fill_value=0)); print(N[N.use_for_miniML=="Y"].groupby("miniML_directions").size())
