import re, glob, os, pandas as pd, numpy as np
base="/mnt/user-data/uploads/Lab Journal/Ephys"
files=sorted(glob.glob(base+"/1. Practice/*.md")+glob.glob(base+"/Longer recordings/*.md"))
def num(s):
    s=s.strip()
    m=re.search(r'-?\d+(?:\.\d+)?',s)
    return float(m.group()) if m else np.nan
rows=[]
for fp in files:
    nm=os.path.basename(fp)
    date=re.match(r'(\d{4})-(\d\d)-(\d\d)',nm); date=int(''.join(date.groups()))
    sl=1; cell=None; chk=0
    for ln in open(fp,encoding='utf8'):
        h=re.match(r'#+\s*Slice\s*(\d+)',ln,re.I)
        if h: sl=int(h.group(1)); continue
        h=re.match(r'#+\s*Cell\s*(\d+)',ln,re.I)
        if h: cell=int(h.group(1)); chk=0; continue
        if ln.startswith('|') and cell is not None:
            c=[x.strip() for x in ln.strip().strip('|').split('|')]
            if len(c)<6 or c[0].lower().startswith('protocol') or (all(re.fullmatch(r'[-: ]*',x) for x in c) and any('-' in x for x in c)): continue
            if not c[0].lower().startswith('manual') : 
                if not any(c[1:]): continue
            chk+=1
            raw=c[1:6]
            if not any(raw): continue
            rows.append(dict(date=date,slice=sl,cell=cell,check=chk,pip=num(raw[0]),Cm=num(raw[1]),Rm=num(raw[2]),Ra=num(raw[3]),Vm=num(raw[4]),Ra_raw=raw[3],Rm_raw=raw[2],Vm_raw=raw[4],Cm_raw=raw[1]))
d=pd.DataFrame(rows)
d['Vm']=-d.Vm.abs()
d.to_csv('manual_checks_journal.csv',index=False)
print(len(d)); print(d.groupby(['date','cell']).size().unstack(fill_value=0))
print(d[(d.Ra_raw.str.contains(r'[^\d\.]'))|(d.Rm_raw.str.contains(r'[^\d\.]'))|(d.Cm_raw.str.contains(r'[^\d\.]'))][['date','cell','check','Cm_raw','Rm_raw','Ra_raw']])
