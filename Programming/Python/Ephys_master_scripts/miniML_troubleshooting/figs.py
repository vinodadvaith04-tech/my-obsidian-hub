import sys; sys.path.insert(0,'.')
import pandas as pd, numpy as np, pickle, pyabf
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':13,'axes.spines.top':False,'axes.spines.right':False,'axes.titlesize':14,'axes.labelsize':13,'legend.frameon':False,'font.family':'DejaVu Sans'})
M_C='#b5451b'; D_C='#2f5d8a'; GR='#8a8a8a'
t=pd.read_pickle('tru.pkl'); key=pd.read_pickle('key.pkl')
# ---------- Fig 1: recall by group, random windows
r=pd.read_csv('summary_random.csv'); a=pd.read_csv('summary_all.csv')
names={'Cs-gluconate up':'Cs-gluconate\noutward (IPSC-like)','Cs-gluconate down':'Cs-gluconate\ninward (EPSC-like)','KGlu down':'K-gluconate\ninward','KGlu up':'K-gluconate\noutward','KCl':'KCl\ninward (IPSC)'}
fig,ax=plt.subplots(figsize=(11,5.2))
x=np.arange(len(r)); w=0.36
b1=ax.bar(x-w/2,r.rec_M,w,color=M_C,label='miniML'); b2=ax.bar(x+w/2,r.rec_D2,w,color=D_C,label='Second detector (matched filter + fit)')
for i,row in r.iterrows():
    ax.text(x[i]-w/2,row.rec_M+0.02,f'{row.rec_M*100:.0f}%',ha='center',fontsize=12,color=M_C)
    ax.text(x[i]+w/2,row.rec_D2+0.02,f'{row.rec_D2*100:.0f}%',ha='center',fontsize=12,color=D_C)
ax.set_xticks(x); ax.set_xticklabels([names[g]+f'\nn={n} events' for g,n in zip(r.group,r.truth)])
ax.set_ylim(0,1.12); ax.set_ylabel('Recall (share of audited events found)'); ax.legend(loc='upper right',ncol=2,bbox_to_anchor=(1,1.12))
ax.set_title('Recall by recording type, random audit windows',loc='left',pad=26)
fig.tight_layout(); fig.savefig('Fig1_recall_by_type.png',dpi=170); plt.close(fig)
# ---------- Fig 2: recall vs rise time and vs SNR
df=pd.read_pickle('/home/claude/work/audit2/truth_feat2.pkl')
df['grp']=np.where(df.internal=='KCl','KCl',df.internal+' '+df.pol)
df['gg']=np.select([df.grp=='Cs-gluconate up',df.grp.isin(['Cs-gluconate down','KGlu down']),df.grp=='KCl'],['Cs outward (n={})','Cs and K-glu inward (n={})','KCl inward (n={})'],'other')
df=df[df.gg!='other']
df['rb']=pd.cut(df.rise_ms,[0,1,2,3,5,100],labels=['<1','1–2','2–3','3–5','>5'])
df['sb']=pd.cut(df.snr,[0,4,5,6,8,100],labels=['<4','4–5','5–6','6–8','>8'])
order=['Cs-gluconate up','inw','KCl']
df['panel']=np.select([df.grp=='Cs-gluconate up',df.grp=='KCl'],[0,2],1)
titles=['Cs-gluconate outward','Cs- and K-gluconate inward','KCl inward']
fig,axs=plt.subplots(2,3,figsize=(14,7.4),sharey=True)
tab=[]
for j in range(3):
    d=df[df.panel==j]
    for i,(c,lab) in enumerate((('rb','Rise time of the event (ms, 10–90%)'),('sb','Event SNR (amplitude / local noise)'))):
        ax=axs[i,j]; g=d.groupby(c,observed=True).agg(n=('found_M','size'),M=('found_M','mean'),D=('found_O','mean'))
        xs=np.arange(len(g)); ax.plot(xs,g.M,'-o',color=M_C,lw=2.5,ms=8,label='miniML'); ax.plot(xs,g.D,'-s',color=D_C,lw=2.5,ms=8,label='Second detector')
        ax.set_xticks(xs); ax.set_xticklabels([f'{k}\nn={n}' for k,n in zip(g.index,g.n)],fontsize=11); ax.set_ylim(0,1.05); ax.set_xlabel(lab)
        if j==0: ax.set_ylabel('Recall')
        if i==0: ax.set_title(titles[j],loc='left')
        for k_,row in g.iterrows(): tab.append(dict(panel=titles[j],feature=c,bin=k_,n=row.n,miniML=row.M,second=row.D))
axs[0,0].legend(loc='lower left')
fig.tight_layout(); fig.savefig('Fig2_recall_vs_rise_and_snr.png',dpi=170); plt.close(fig)
pd.DataFrame(tab).to_csv('fig2_data.csv',index=False)
# ---------- Fig 3: parameter sweep (Cs outward)
sw=pd.DataFrame([
 ('default (window 600, threshold 0.5)',.485,.992),('threshold 0.7',.389,1.0),('threshold 0.3',.515,.918),('threshold 0.2',.542,.904),('threshold 0.1',.584,.838),
 ('window 150',.019,.5),('window 300',.153,.771),('window 1200',.576,.903),('window 2400',.527,.908),
 ('window 1200 + threshold 0.3',.649,.885),('window 1200 + threshold 0.2',.672,.859),('peak width 10',.366,1.0),('2 kHz low-pass first',.466,.976)],columns=['cfg','rec','prec'])
sw['n']=np.arange(1,len(sw)+1)
fig=plt.figure(figsize=(14,6)); ax=fig.add_axes([0.06,0.13,0.5,0.78]); lg=fig.add_axes([0.6,0.13,0.39,0.78]); lg.axis('off')
ax.axvspan(0.8,1.0,color='#e9eef5',zorder=0); ax.text(0.9,0.52,'KCl recall\nabout 0.89',ha='center',fontsize=11,color=D_C)
nud={1:(0,0),2:(-0.012,0.012),12:(0.0,-0.0),13:(0,0),3:(0,0),9:(0,0),4:(0,0),8:(0,0)}
for _,r_ in sw.iterrows():
    c=M_C if r_.n==1 else '#555'
    ax.scatter([r_.rec],[r_.prec],s=250,color=c,zorder=3)
    ax.text(r_.rec,r_.prec,str(r_.n),ha='center',va='center',color='white',fontsize=10.5,zorder=4,fontweight='bold')
for i,(_,r_) in enumerate(sw.iterrows()):
    lg.text(0,0.97-i*0.07,f'{r_.n}',fontsize=12,fontweight='bold',color=M_C if r_.n==1 else '#555',va='center')
    lg.text(0.07,0.97-i*0.07,f'{r_.cfg}:  recall {r_.rec+1e-9:.2f}, precision {r_.prec+1e-9:.2f}',fontsize=11,va='center')
ax.set_xlim(0,1); ax.set_ylim(0.45,1.05); ax.set_xlabel('Recall (Cs outward events)'); ax.set_ylabel('Precision')
fig.suptitle('Every miniML setting tried trades precision for recall (Cs outward events, 262 audited events)',x=0.07,ha='left',fontsize=14)
fig.savefig('Fig3_parameter_sweep_Cs_outward.png',dpi=170); plt.close(fig)
sw.to_csv('fig3_data.csv',index=False)
# ---------- Fig 4: network probability of missed events
x=pd.read_csv('pmax_bands.csv',index_col=0)
order=['Cs-gluconate up','Cs-gluconate down','KGlu down','KGlu up','KCl']; x=x.loc[order]
p=x.div(x.sum(axis=1),axis=0)
cols=['#7a1f0c','#c8663f','#e8b08f','#9db4cc']; labs=['< 0.1 (network sees nothing)','0.1–0.3','0.3–0.5 (a lower threshold could catch)','≥ 0.5 (network sees it, lost afterwards)']
fig,ax=plt.subplots(figsize=(11,5.2)); left=np.zeros(len(p))
for c,col,lb in zip(x.columns,cols,labs):
    ax.barh(np.arange(len(p)),p[c],left=left,color=col,label=lb,edgecolor='white',linewidth=1.5)
    for i,(v,l) in enumerate(zip(p[c],left)):
        if v>0.06: ax.text(l+v/2,i,f'{v*100:.0f}%',ha='center',va='center',color='white' if col in ('#7a1f0c','#c8663f') else '#222',fontsize=12)
    left+=p[c].values
ax.set_yticks(np.arange(len(p))); ax.set_yticklabels([names[g].replace('\n',' ')+f'  (n={int(n)} missed)' for g,n in zip(p.index,x.sum(axis=1))]); ax.invert_yaxis()
ax.set_xlim(0,1); ax.set_xlabel('Share of missed events (by highest classifier output within ±50 ms)')
ax.legend(loc='upper center',bbox_to_anchor=(0.45,-0.2),ncol=2,fontsize=11)
ax.set_title("What miniML's classifier outputs at the events it misses",loc='left')
fig.tight_layout(); fig.savefig('Fig4_network_output_at_misses.png',dpi=170); plt.close(fig)
# ---------- Fig 5: example window
win=pd.read_excel('/mnt/user-data/uploads/Ephys pilot & practice/audit/audit_labels.xlsx','audit_windows').set_index('window')
fl=pd.read_pickle('/home/claude/work/audit2/master_sheets.pkl')['files'].set_index('new_name')
wid='W203'; rw=win.loc[wid]; nm=rw.file.replace('_Adv6_1B12_','_Adv3_1B12_')
ROOT='/mnt/user-data/uploads/Ephys pilot & practice/'
abf=pyabf.ABF(ROOT+fl.loc[nm,'path_on_mac']); abf.setSweep(0); y=abf.sweepY; fs=abf.dataRate; tt=np.arange(len(y))/fs
pr=pickle.load(open(f'/home/claude/work/audit2/preds/{wid}.pkl','rb'))['res']['positive']
a_,b_=rw.start_s,rw.end_s
m=(tt>=a_)&(tt<=b_); ptime=pr['off']+np.arange(len(pr['pred']))*pr['dt']; pm=(ptime>=a_)&(ptime<=b_)
T=t[t.window==wid]; Tu=T[T.polarity=='up']; K=key[key.window==wid]
Mk=K[(K.detector=='miniML')&(K.polarity=='up')&(~K.false)]; Dk=K[(K.detector!='miniML')&(K.polarity=='up')&(~K.false)]
import scipy.signal as ss
yl=ss.sosfiltfilt(ss.butter(4,1000,'low',fs=fs,output='sos'),y)
fig,(ax1,ax2)=plt.subplots(2,1,figsize=(13,6.4),sharex=True,gridspec_kw={'height_ratios':[2.3,1]})
ax1.plot(tt[m],yl[m],color='#444',lw=0.9)
top=np.nanmax(yl[m])+8
ax1.plot(Tu.time_s,[top+14]*len(Tu),'|',color='k',ms=16,mew=1.8,label=f'audited events (n={len(Tu)})')
ax1.plot(Mk.time_s,[top+7]*len(Mk),'v',color=M_C,ms=8,label=f'miniML markers (n={len(Mk)})')
ax1.plot(Dk.time_s,[top]*len(Dk),'v',color=D_C,ms=8,label=f'second detector markers (n={len(Dk)})')
ax1.set_ylabel('Current (pA, 1 kHz low-pass)'); ax1.legend(loc='upper center',ncol=3,bbox_to_anchor=(0.5,1.14),fontsize=11)
ax2.plot(ptime[pm],pr['pred'][pm],color=M_C,lw=1.4); ax2.axhline(0.5,color='#333',ls='--',lw=1); ax2.text(a_+0.02,0.53,'default threshold 0.5',fontsize=10.5)
ax2.set_ylim(0,1.05); ax2.set_ylabel('miniML output'); ax2.set_xlabel('Time (s)'); ax2.set_xlim(a_,b_)
fig.suptitle(f'Cs-gluconate, {int(rw.holding_mV)} mV: classifier output stays flat at many real outward IPSCs ({wid})',x=0.01,ha='left',fontsize=14,y=0.995)
fig.tight_layout(rect=(0,0,1,0.96)); fig.savefig('Fig5_example_window.png',dpi=170); plt.close(fig)
print('ok', len(Tu), len(Mk), len(Dk), nm)
