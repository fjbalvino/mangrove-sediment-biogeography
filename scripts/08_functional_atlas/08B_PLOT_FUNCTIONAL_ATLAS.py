#!/usr/bin/env python3
"""Draw Figure 4 and the functional-count sensitivity figure from saved results."""
from pathlib import Path
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT=Path(__file__).resolve().parents[2]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--results',type=Path,default=ROOT/'results/functional')
parser.add_argument('--output',type=Path,default=ROOT/'figures/functional')
args=parser.parse_args()
R=args.results
OUT=args.output;OUT.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelsize':12,
                     'axes.titlesize':12,'axes.linewidth':.8,'xtick.labelsize':10,'ytick.labelsize':10,
                     'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
realm_codes={'Tropical Atlantic':'TA','Central Indo-Pacific':'CIP','Western Indo-Pacific':'WIP','Tropical Eastern Pacific':'TEP'}
realm_colors={'TA':'#E78AC3','CIP':'#52B698','WIP':'#E6BB26','TEP':'#A6D854'}
class_order=['Absolute core','Near-global','Realm-core','Periphery']
class_colors={'Absolute core':'#7F0000','Near-global':'#C62828','Realm-core':'#E57373','Periphery':'#B7B7B7'}
classes=pd.read_csv(R/'functional_class_summary.tsv',sep='\t').set_index('class').loc[class_order]
scores=pd.read_csv(R/'PCA_sample_scores.tsv',sep='\t')
stats=json.loads((R/'functional_results_summary.json').read_text())
pca=stats['PCA'];coupling=stats['coupling'];decay=stats['decay']['metagenome']
pairs=pd.read_csv(R/'coupling_pairs.tsv.gz',sep='\t')
decay_pairs=pd.read_csv(R/'metagenome_distance_decay_pairs.tsv.gz',sep='\t')
perf=pd.read_csv(R/'functional_PERMANOVA.tsv',sep='\t')
pair=perf.loc[(perf.comparison=='Tropical Atlantic vs Central Indo-Pacific')&(perf.geometry=='Aitchison')].iloc[0]

def clean(ax,grid=None):
    ax.spines[['top','right']].set_visible(False)
    if grid:ax.grid(axis=grid,color='#E8E8E8',lw=.55);ax.set_axisbelow(True)

def label(ax,text,x=-.115):
    ax.text(x,1.04,text,transform=ax.transAxes,fontsize=18,fontweight='bold',va='bottom')

fig=plt.figure(figsize=(13.12,8.64))
gs=fig.add_gridspec(2,2,left=.13,right=.98,bottom=.085,top=.94,wspace=.31,hspace=.43)
a=fig.add_subplot(gs[0,0]); b=fig.add_subplot(gs[0,1]);c=fig.add_subplot(gs[1,0]);d=fig.add_subplot(gs[1,1])

# Both percentages are visible even when the richness bar is very short.
y=np.arange(4)
for i,cls in enumerate(class_order):
    row=classes.loc[cls];color=class_colors[cls]
    a.barh(i-.15,row.KO_richness_percent,height=.25,color=color,alpha=.34)
    a.barh(i+.15,row.KO_assignment_percent,height=.25,color=color)
    a.text(row.KO_richness_percent+1.1,i-.15,f'{row.KO_richness_percent:.2f}%' if row.KO_richness_percent<1 else f'{row.KO_richness_percent:.1f}%',va='center',fontsize=10)
    a.text(row.KO_assignment_percent+1.1,i+.15,f'{row.KO_assignment_percent:.1f}%',va='center',fontsize=10)
a.set_yticks(y,[f'{cls}\n{int(classes.loc[cls,"KOs"]):,} KOs ({classes.loc[cls,"KO_assignment_percent"]:.1f}%)' for cls in class_order])
a.set_ylim(3.55,-.9);a.set_xlim(0,102);a.set_xticks([0,25,50,75,100]);a.set_xlabel('Share of richness or KO annotation counts (%)')
a.tick_params(axis='y',length=0,labelsize=10.7)
a.legend(handles=[Patch(facecolor='#999999',alpha=.34,label='KO richness'),Patch(facecolor='#999999',label='KO annotation counts')],
         frameon=False,ncol=1,loc='upper right',fontsize=9.6,handlelength=1.2,borderaxespad=0)
clean(a,'x');label(a,'a',-.34)

# PCA maintains the original CLR(count + 1) definition; its count dependence is
# made explicit by the accompanying diagnostic figure and manuscript caveat.
for realm in ['Tropical Atlantic','Tropical Eastern Pacific','Western Indo-Pacific','Central Indo-Pacific']:
    sub=scores[scores.REALM==realm];code=realm_codes[realm]
    b.scatter(sub.PC1,sub.PC2,s=35,c=realm_colors[code],edgecolors='#444444',lw=.5,alpha=.84,label=f'{code} (n = {len(sub)})')
b.set_xlabel(f'PC1 ({pca["PC1_variance_percent"]:.1f}%)');b.set_ylabel(f'PC2 ({pca["PC2_variance_percent"]:.1f}%)')
b.legend(frameon=False,loc='upper left',fontsize=9,handletextpad=.25,borderaxespad=.2)
b.margins(x=.08,y=.3)
b.text(.98,.025,f'TA–CIP: site-level PERMANOVA\nR² = {pair.R2:.3f}; Holm P = {pair.Holm_P_pairwise_within_geometry:.4f}',transform=b.transAxes,ha='right',va='bottom',fontsize=9.4,
       bbox={'facecolor':'white','edgecolor':'none','alpha':.87,'pad':2})
inset=b.inset_axes([.52,.84,.46,.10]);left=0
for cls in class_order:
    val=pca['top100_absolute_PC1_loading_percent_by_class'][cls]
    inset.barh(0,val,left=left,color=class_colors[cls],height=.8)
    if val>=15:inset.text(left+val/2,0,f'{val:.1f}%',ha='center',va='center',fontsize=8.5,color='white' if cls=='Near-global' else '#222222')
    left+=val
inset.set_xlim(0,100);inset.set_ylim(-.55,.55);inset.axis('off')
inset.set_title('Top 100 absolute PC1 loadings',fontsize=8.8,pad=2)
clean(b);label(b,'b')

c.scatter(pairs.taxonomic_Aitchison,pairs.functional_Aitchison,s=4,c='#6E8299',alpha=.08,rasterized=True,edgecolors='none')
fit=np.polyfit(pairs.taxonomic_Aitchison,pairs.functional_Aitchison,1)
xx=np.linspace(pairs.taxonomic_Aitchison.min(),pairs.taxonomic_Aitchison.max(),200)
c.plot(xx,np.polyval(fit,xx),color='#202020',lw=1.6)
c.set_xlabel('Taxonomic Aitchison dissimilarity');c.set_ylabel('Functional Aitchison dissimilarity')
c.text(.03,.98,f'Mantel r = {coupling["metagenome_level"]["r"]:.3f}; P = {coupling["metagenome_level"]["P"]:.4f}\nSite-centroid r = {coupling["site_centroid_level"]["r"]:.3f}; P = {coupling["site_centroid_level"]["P"]:.4f}',transform=c.transAxes,va='top',fontsize=10.3,
       bbox={'facecolor':'white','edgecolor':'none','alpha':.88,'pad':2})
c.text(.97,.025,f'{stats["functional_profiles"]} matched metagenomes',transform=c.transAxes,ha='right',fontsize=9.5,color='#555555')
c.margins(.03);clean(c);label(c,'c',-.34)

colors={'taxonomic':'#0072B2','functional':'#D55E00'}
xx=np.linspace(decay_pairs.log10_distance_km.min(),decay_pairs.log10_distance_km.max(),250)
for name in ['taxonomic','functional']:
    d.scatter(decay_pairs.log10_distance_km,decay_pairs[name+'_similarity'],s=3,c=colors[name],alpha=.045,edgecolors='none',rasterized=True)
    d.plot(xx,decay['intercept_'+name]+decay['slope_'+name]*xx,color=colors[name],lw=2,label=name.capitalize())
d.set_xticks([-1,0,1,2,3,4],['0.1','1','10','100','1,000','10,000']);d.set_xlabel('Geographic distance (km, log scale)')
d.set_ylabel('Bray–Curtis similarity');d.set_ylim(0,1.22);d.set_yticks([0,.2,.4,.6,.8,1])
d.set_xlim(min(xx)-.06,max(xx)+.06)
d.text(.03,.98,f'Taxonomic slope = {decay["slope_taxonomic"]:.4f}; P = {decay["P_taxonomic_slope_two_sided"]:.4f}',transform=d.transAxes,va='top',color=colors['taxonomic'],fontsize=10)
d.text(.03,.90,f'Functional slope = {decay["slope_functional"]:.4f}; P = {decay["P_functional_slope_two_sided"]:.4f}',transform=d.transAxes,va='top',color=colors['functional'],fontsize=10)
d.text(.03,.82,f'Δ slope = {decay["difference_tax_minus_function"]:.4f}; P = {decay["P_difference_two_sided"]:.4f}',transform=d.transAxes,va='top',fontsize=10)
d.text(.98,.025,f'{stats["georeferenced_profiles"]} metagenomes · {stats["coordinate_sites"]} sites',transform=d.transAxes,ha='right',fontsize=9.5,color='#555555')
clean(d);label(d,'d')
fig.savefig(OUT/'Figure_4.png',dpi=500)
fig.savefig(OUT/'Figure_4.pdf',dpi=300)
plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(12,4.7),gridspec_kw={'wspace':.33})
fig.subplots_adjust(left=.085,right=.98,bottom=.15,top=.9)
for realm in ['Tropical Atlantic','Tropical Eastern Pacific','Western Indo-Pacific','Central Indo-Pacific']:
    sub=scores[scores.REALM==realm];code=realm_codes[realm]
    axs[0].scatter(np.log10(sub.KO_assignment_total),sub.PC1,s=31,c=realm_colors[code],edgecolors='#444444',lw=.45,label=code)
    axs[1].scatter(sub.Hellinger_PC1,sub.Hellinger_PC2,s=31,c=realm_colors[code],edgecolors='#444444',lw=.45,label=code)
axs[0].set_xlabel('Total KO annotation counts (log₁₀)');axs[0].set_ylabel('PC1 of CLR-transformed KO counts')
axs[0].text(.03,.97,f'Spearman ρ = {pca["PC1_spearman_rho_KO_assignment_total"]:.3f}\nn = {stats["functional_profiles"]} metagenomes',transform=axs[0].transAxes,va='top',fontsize=10.5)
axs[0].legend(frameon=False,loc='lower right',fontsize=10)
axs[1].set_xlabel(f'Hellinger PC1 ({pca["Hellinger_PC1_variance_percent"]:.1f}%)');axs[1].set_ylabel(f'Hellinger PC2 ({pca["Hellinger_PC2_variance_percent"]:.1f}%)')
axs[1].text(.03,.97,'Sensitivity using square-root relative abundances',transform=axs[1].transAxes,va='top',fontsize=10,
            bbox={'facecolor':'white','edgecolor':'none','alpha':.85})
axs[1].margins(.15,.19)
for ax,letter in zip(axs,['a','b']):clean(ax);label(ax,letter)
fig.savefig(OUT/'Supplementary_Figure_S4.png',dpi=400)
fig.savefig(OUT/'Supplementary_Figure_S4.pdf',dpi=300)
plt.close(fig)
print(OUT/'Figure_4.png')
