"""A compact biochemical-marker panel; selection is not based on significance."""
from pathlib import Path
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
r=Path(__file__).resolve().parent/'atlas';src=r/'results';(r/'supplement').mkdir(exist_ok=True)
a=pd.read_csv(src/'relative_representation.tsv.gz',sep='\t',index_col=0);meta=pd.read_csv(src/'sample_analysis_manifest.tsv',sep='\t',index_col=0);order=pd.read_csv(src/'heatmap_sample_key.tsv',sep='\t',index_col=0).index
panels=[('Carbon',['K00855','K01601','K15230','K00198','K00448']),('Methane / C1',['K00399','K10944','K16157','K14028']),('Sulfur',['K11180','K00394','K17218','K17222','K17486']),('Phosphorus',['K02040','K01077','K06163','K18916','K00937']),('Iron',['K04758','K04759']),('Hydrogen',['K00440','K03617']),('Respiration / light',['K00404','K00425','K02689','K08928'])]
inv=pd.read_csv(src/'KO_inventory.tsv',sep='\t').set_index('KO');ks=[k for c,l in panels for k in l];data=a.loc[ks,order].to_numpy();cmap=plt.get_cmap('plasma');import json
vmax=json.loads((src/'viewer_manifest.json').read_text())['colour_max_log1p'];norm=Normalize(0,vmax);rgba=cmap(norm(np.log10(1+np.nan_to_num(data))));rgba[data==0]=[1,1,1,1];rgba[np.isnan(data)]=[.776,.8,.824,1]
fig,ax=plt.subplots(figsize=(10,8.6));ax.imshow(rgba,aspect='auto',interpolation='nearest');names={k:inv.loc[k,'gene'].split(',')[0] for k in ks};names.update({'K10944':'Amo/Pmo A','K00198':'CooS/AcsA','K18916':'PtxD'})
ax.set_yticks(range(len(ks)),[names[k]+'  '+k for k in ks],fontsize=8.5);ticks=np.linspace(0,len(order)-1,min(12,len(order)),dtype=int);ax.set_xticks(ticks,[f'{meta.loc[order[i],"lon"]:.0f}°' for i in ticks],fontsize=9);ax.set_xlabel(f'{len(order)} metagenomes · decreasing signed longitude (east → west)',fontsize=10)
start=0
for cycle,k in panels:
    end=start+len(k);ax.axhline(end-.5,color='#a8a2af',lw=.6);ax.text(1.01,(start+end-1)/2,cycle,transform=ax.get_yaxis_transform(),fontsize=9,va='center');start=end
fig.subplots_adjust(left=.22,right=.83,top=.94,bottom=.15);fig.suptitle('Biogeographic representation of selected functional markers',fontsize=13,y=.985)
bar=fig.add_axes([.27,.055,.45,.018]);cb=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=bar,orientation='horizontal',ticks=np.arange(vmax+1));cb.ax.set_xticklabels([f'{10**i-1:,.0f}' for i in np.arange(vmax+1)]);cb.set_label('KO assignments per million total KO assignments',fontsize=9)
for ext in ['png','pdf']:fig.savefig(r/'supplement'/('Selected_marker_geography.'+ext),dpi=300)
plt.close(fig);pd.DataFrame([{'cycle':c,'KO':k} for c,kk in panels for k in kk]).to_csv(r/'supplement/Selected_marker_gene_key.tsv',sep='\t',index=False)

# Environmental design panel (the method used for reference Supplementary Fig. S8).
design=r.parent/'analysis_nitrogen/design_review/environmental_information_by_study.tsv'
try:
    d=pd.read_csv(design,sep='\t')
except pd.errors.EmptyDataError:
    d=pd.DataFrame()
if not d.empty:
    d=d[d.context=='joint_complete_cases'];wide=d.pivot(index='study_id',columns='predictor',values='raw_within_variance_share').fillna(0)
    wide=wide.loc[wide.max(axis=1)>1e-12]
    fig,axes=plt.subplots(1,3,figsize=(12,max(4,len(wide)*.42)),sharey=True)
    for ax,column,colour in zip(axes,['Temperature','Organic carbon','pH'],['#6a00a8','#cb4679','#fca636']):
        ax.barh(wide.index,wide[column]*100,color=colour);ax.set_xlim(0,100);ax.set_title(column);ax.set_xlabel('Within-study variation (%)')
    axes[0].invert_yaxis();fig.tight_layout()
    for ext in ['png','pdf']:fig.savefig(r/'supplement'/('Environmental_design_support.'+ext),dpi=300)
    plt.close(fig)
