import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
m=pd.read_csv('k60_master_table.csv')
m=m[m.label_status!='contradicted by cap check'].copy()
m['f']=np.where(m.grp=='-60 mV',m.f_neg,m.f_pos)
m['cellkey']=m.date.astype(str)+'_s'+m.slice.astype(str)+'_c'+m.cell.astype(str)
c=m.groupby(['grp','cellkey']).agg(f=('f','mean'),nrec=('f','size')).reset_index()
c.to_csv('k60_event_freq_per_cell.csv',index=False)
grp=['-60 mV','0 mV','+5 to +15 mV']; cols={'-60 mV':'#1f4e79','0 mV':'#b5451b','+5 to +15 mV':'#c9892b'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(1,2,figsize=(9.4,4.6)); rows=[]
for pi,ax in enumerate(axs):
  rng=np.random.default_rng(3)
  for k,g in enumerate(grp):
      v=c[c.grp==g].f.values; q1,med,q3=np.percentile(v,[25,50,75]); mn,mx=v.min(),v.max(); mu=v.mean()
      rows.append((g,len(v),mu,med,q1,q3,mn,mx)) if pi==0 else None; col=cols[g]; w=.28
      ax.add_patch(plt.Rectangle((k-w,q1),2*w,q3-q1,fc=col,alpha=.18,ec=col,lw=1.2,zorder=1))
      ax.plot([k-w,k+w],[med,med],color=col,lw=2.4,zorder=3)
      ax.plot([k,k],[mn,q1],color=col,lw=1.2,zorder=1); ax.plot([k,k],[q3,mx],color=col,lw=1.2,zorder=1)
      for y in (mn,mx): ax.plot([k-w/2,k+w/2],[y,y],color=col,lw=1.2,zorder=1)
      ax.scatter([k],[mu],marker='D',s=42,fc='white',ec=col,lw=1.6,zorder=5)
      ax.scatter(k+rng.uniform(-.14,.14,len(v)),v,s=22,color=col,alpha=.8,lw=.4,ec='white',zorder=4)
      ax.text(k,-0.035,f'n = {len(v)} cells',ha='center',va='top',fontsize=8,transform=ax.get_xaxis_transform(),color='#333')
  ax.set_xticks(range(3)); ax.set_xticklabels(['-60','0','+5 to +15']); ax.tick_params(axis='x',pad=14)
  ax.set_xlabel('Command holding (mV)',labelpad=4); ax.set_xlim(-.6,2.6)
  if pi==0:
      ax.set_ylim(0,8); ax.set_ylabel('miniML event frequency (Hz)'); ax.set_title('A  Full range',loc='left',fontsize=9)
      ax.annotate('6 Nov C2\n7.5 Hz',xy=(0.05,7.51),xytext=(0.55,7.0),fontsize=8,arrowprops=dict(arrowstyle='-',color='#777',lw=.8),va='center')
  else:
      ax.set_ylim(0,2.0); ax.set_ylabel('miniML event frequency (Hz)'); ax.set_title('B  Zoom to 0-2 Hz (same data)',loc='left',fontsize=9)
      ax.legend(handles=[Line2D([0],[0],color='#555',lw=2.4,label='Median'),Line2D([0],[0],marker='D',color='w',mfc='white',mec='#555',mew=1.6,ms=6,label='Mean'),
        Line2D([0],[0],marker='s',color='w',mfc='#ccc',mec='#555',ms=8,label='Box = IQR'),Line2D([0],[0],color='#555',lw=1.2,label='Whiskers = min, max'),
        Line2D([0],[0],marker='o',color='w',mfc='#555',ms=5,label='One cell')],frameon=False,fontsize=8,loc='upper right')
fig.tight_layout(); fig.savefig('fig2.png',dpi=220)
pd.DataFrame(rows,columns=['grp','n_cells','mean','median','q1','q3','min','max']).round(3).to_csv('k60_event_freq_box_stats.csv',index=False)
print(pd.DataFrame(rows,columns=['grp','n','mean','median','q1','q3','min','max']).round(2)); print(c.round(2).to_string())
