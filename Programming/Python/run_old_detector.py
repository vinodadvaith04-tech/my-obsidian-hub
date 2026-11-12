"""Run the earlier matched-filter + fit detector (psc_detection.py) on all master recordings (use_for_miniML == Y).
Holdings, internal and cell labels come ONLY from the new master table (command holding); nothing from the old holdings.csv.
Resumable (skips files already done), stop file: <out>/STOP, newest recordings first."""
import sys, time, json, argparse, glob
from pathlib import Path
import numpy as np, pandas as pd, pyabf
ap=argparse.ArgumentParser(); ap.add_argument("--master",required=True); ap.add_argument("--data-root",required=True)
ap.add_argument("--psc-dir",required=True); ap.add_argument("--out",required=True); ap.add_argument("--max-hours",type=float)
a=ap.parse_args(); sys.path.insert(0,a.psc_dir); import psc_detection as P
out=Path(a.out); (out/"per_file").mkdir(parents=True,exist_ok=True)
m=pd.read_excel(a.master,sheet_name="files"); m=m[m.use_for_miniML=="Y"].copy()
m["_c"]=pd.to_numeric(m["counter"],errors="coerce"); m=m.sort_values(["date_header","_c"],ascending=False)
idx={p.name:p for p in Path(a.data_root).rglob("*.abf")}
t0=time.time()
for i,(_,r) in enumerate(m.iterrows(),1):
    if (out/"STOP").exists(): print("STOP file found"); break
    if a.max_hours and time.time()-t0>a.max_hours*3600: print("max-hours reached"); break
    f=out/"per_file"/(r.new_name[:-4]+".old.json")
    if f.exists(): continue
    ts=time.time(); ab=pyabf.ABF(str(idx[r.new_name])); ab.setSweep(0); I=np.asarray(ab.sweepY,float); fs=float(ab.dataRate)
    ok,rej,di=P.detect(I,fs)
    ev=[dict(new_name=r.new_name,polarity=e["polarity"],peak_s=e["peak_s"],amp_pA=e["amp_pA"],rise1090_ms=e["rise1090_ms"],tau_decay_ms=e["td"],halfwidth_ms=e["halfwidth_ms"],snr=e["snr"],isolated=bool(e["isolated"])) for e in ok]
    dur=di["duration_s"]-di["masked_s"]
    summ=dict(new_name=r.new_name,holding_command_mV=r.holding_command_mV,internal=r.internal,duration_s=di["duration_s"],analysed_s=dur,noise_lp_pA=di["noise_lp_pA"],
        limit_inward_pA=di["limit_inward_pA"],limit_outward_pA=di["limit_outward_pA"],n_inward=sum(e["polarity"]=="inward" for e in ev),n_outward=sum(e["polarity"]=="outward" for e in ev))
    summ["freq_inward_Hz"]=summ["n_inward"]/dur; summ["freq_outward_Hz"]=summ["n_outward"]/dur
    pd.DataFrame(ev).to_csv(out/"per_file"/(r.new_name[:-4]+".old.events.csv"),index=False); f.write_text(json.dumps(summ,default=float))
    print(f"[{i}/{len(m)}] {r.new_name} {time.time()-ts:.0f}s",flush=True)
S=[json.load(open(x)) for x in glob.glob(str(out/"per_file"/"*.old.json"))]
pd.DataFrame(S).to_csv(out/"oldDetector_file_summary.csv",index=False)
ev=[pd.read_csv(x) for x in glob.glob(str(out/"per_file"/"*.old.events.csv")) if Path(x).stat().st_size>5]
if ev: pd.concat(ev).to_csv(out/"oldDetector_events_all.csv",index=False)
print("done",len(S),"recordings")
