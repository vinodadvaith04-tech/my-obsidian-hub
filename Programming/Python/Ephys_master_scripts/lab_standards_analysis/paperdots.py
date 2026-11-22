import pymupdf as fitz, json, numpy as np
d=fitz.open('/mnt/user-data/uploads/3c7df21c-5483.full.pdf')
S=72/150
panels={ # name: (page idx, bbox px x0,y0,x1,y1, split_x px between CTL and VU)
 'DIV9_sEPSC':(7,(95,340,350,505),None),
 'DIV9_sIPSC':(7,(95,555,350,720),None),
 'DIV21_sEPSC':(10,(95,455,350,610),None),
 'DIV21_sIPSC':(10,(95,685,350,840),None),
}
out={}
for name,(pn,bb,_) in panels.items():
    pg=d[pn]; x0,y0,x1,y1=[v*S for v in bb]
    R=fitz.Rect(x0,y0,x1,y1)
    dr=pg.get_drawings()
    circ=[]
    for x in dr:
        r=x['rect']
        if R.contains(r) and x['type']=='fs' and len(x['items'])==4 and 2<r.width<6 and abs(r.width-r.height)<0.5:
            circ.append(((r.x0+r.x1)/2,(r.y0+r.y1)/2,x['fill']))
    words=[w for w in pg.get_text('words') if R.contains(fitz.Rect(w[:4]))]
    print(name,len(circ)); print([ (w[4],round((w[1]+w[3])/2,1),round(w[0],1)) for w in words])
    out[name]={'circ':circ,'words':[(w[4],(w[1]+w[3])/2,w[0]) for w in words]}
json.dump(out,open('paper_raw.json','w'))
