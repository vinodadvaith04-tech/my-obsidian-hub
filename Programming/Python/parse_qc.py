import re, pandas as pd, glob, os
D="/mnt/user-data/uploads/Lab Journal/Ephys/3. Cs internal pilot/"
rows=[]
for f in sorted(glob.glob(D+"2026-11-*.md")):
    day=os.path.basename(f)[:10]; cell=None; sl=None
    for line in open(f,encoding="utf-8"):
        m=re.match(r"#+\s*Slice\s*(\d+)",line)
        if m: sl=int(m.group(1))
        m=re.match(r"#+\s*Cell\s*(\d+)",line)
        if m: cell=int(m.group(1)); continue
        m=re.match(r"\|\s*Manual check (\d+)\s*\|(.*)\|\s*$",line)
        if m and cell:
            cols=[c.strip() for c in m.group(2).split("|")]
            num=lambda s: float(s) if re.fullmatch(r"[\d.]+",s or "") else None
            rows.append(dict(day=day,cell=cell,check=int(m.group(1)),pip=num(cols[0]),Cm=num(cols[1]),Rm=num(cols[2]),Ra=num(cols[3])))
q=pd.DataFrame(rows); q=q.dropna(subset=["Cm","Rm","Ra"],how="all"); q.to_csv("cs_qc_manual.csv",index=False)
print(q.groupby(["day","cell"]).agg(n=("check","count"),Ra_first=("Ra","first"),Ra_last=("Ra","last"),Ra_max=("Ra","max"),Rm_first=("Rm","first"),Rm_last=("Rm","last"),Cm_first=("Cm","first"),Cm_last=("Cm","last")))
