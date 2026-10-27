#!/usr/bin/env python3
"""Assemble Figure 5 and two supplementary diagnostics from validated tables."""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

COLORS = {'Prevalent':'#0072B2', 'Shell':'#D55E00', 'Cloud':'#777777'}


def tidy(ax):
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(length=3, labelsize=11)


def title(ax, letter, text):
    ax.set_title(text, loc='left', fontsize=12, fontweight='bold', pad=14)
    ax.text(-.12,1.08,letter, transform=ax.transAxes, fontsize=17,fontweight='bold')


def save(fig, out, name, formats=('png','pdf','tiff')):
    for ext in formats:
        kwargs = {'pil_kwargs':{'compression':'tiff_lzw'}} if ext=='tiff' else {}
        target=out/(name+'.'+ext)
        buffer=io.BytesIO()
        fig.savefig(buffer,format=ext,dpi=600,facecolor='white',**kwargs)
        data=buffer.getvalue()
        assert len(data)>100
        if ext=='pdf':
            assert data.startswith(b'%PDF-') and b'%%EOF' in data[-64:]
        temporary=target.with_suffix(target.suffix+'.partial')
        with temporary.open('wb') as handle:
            handle.write(data);handle.flush();os.fsync(handle.fileno())
        temporary.replace(target)
        assert hashlib.sha256(target.read_bytes()).digest()==hashlib.sha256(data).digest()
    plt.close(fig)


def panel_a(ax,cross):
    persistence=['Absolute core','Near-global','Realm-core','Periphery']
    prevalence=['Prevalent','Shell','Cloud']
    z=cross.pivot(index='persistence_class',columns='prevalence_class',values='n_KOs').reindex(index=persistence,columns=prevalence).fillna(0).values
    cmap=LinearSegmentedColormap.from_list('blue_cells',['#FAFAFA','#D3E8F2','#74AEC8'])
    ax.imshow(np.log10(z+1),cmap=cmap,aspect='auto',vmin=0,vmax=np.log10(z.max()+1))
    for i in range(4):
        for j in range(3):
            ax.text(j,i,f'{int(z[i,j]):,}\n({100*z[i,j]/z.sum():.1f}%)' if z[i,j] else '0',ha='center',va='center',fontsize=11)
    ax.set_xticks(range(3),[f'{g}\n{int(z[:,j].sum()):,} KOs' for j,g in enumerate(prevalence)],fontsize=11)
    ax.set_yticks(range(4),persistence,fontsize=11)
    ax.set_xticks(np.arange(-.5,3,1),minor=True);ax.set_yticks(np.arange(-.5,4,1),minor=True)
    ax.grid(which='minor',color='white',linewidth=2)
    ax.tick_params(which='both',length=0)
    for s in ax.spines.values():s.set_visible(False)
    title(ax,'a','Persistence and prevalence')


def panel_b(ax,df):
    order=['Prevalent','Realm-core Shell','Peripheral Shell','Peripheral Cloud']
    rows=[]
    for g in order:
        z=df[df.group==g]
        if len(z):rows.append((g,z.R2.median(),*z.R2.quantile([.05,.95]),int(z.n_KOs.iloc[0])))
    for i,(g,m,lo,hi,n) in enumerate(rows):
        color=COLORS['Shell' if 'Shell' in g else ('Cloud' if 'Cloud' in g else 'Prevalent')]
        ax.plot([lo,hi],[i,i],color=color,lw=2)
        ax.scatter(m,i,s=50,color=color,zorder=3)
    ax.set_yticks(range(len(rows)),[g.replace(' ','\n',1)+f'\n({n:,} KOs)' if ' ' in g else g+f'\n({n:,} KOs)' for g,_,_,_,n in rows],fontsize=10)
    ax.invert_yaxis();ax.set_ylim(len(rows)-.5,-.5)
    ax.set_xlim(0,max(r[3] for r in rows)*1.1)
    ax.set_xlabel('Aitchison PERMANOVA $R^2$',fontsize=11)
    tidy(ax);title(ax,'b','Regional differentiation')


def panel_prevalence(ax,df):
    groups=['Prevalent','Shell','Cloud']
    pos=[];data=[];colors=[]
    for j,design in enumerate(['3 x 6','3 x 4']):
        for i,g in enumerate(groups):
            pos.append(j*4+i);data.append(df[(df.design==design)&(df.group==g)].R2.values);colors.append(COLORS[g])
    violin=ax.violinplot(data,positions=pos,widths=.8,showextrema=False)
    for body,color in zip(violin['bodies'],colors):body.set_facecolor(color);body.set_alpha(.55);body.set_edgecolor(color)
    ax.boxplot(data,positions=pos,widths=.19,showfliers=False,patch_artist=True,
               boxprops={'facecolor':'white','edgecolor':'#333333'},medianprops={'color':'black'},
               whiskerprops={'color':'#444444'},capprops={'color':'#444444'})
    ax.set_xticks(pos,groups*2,rotation=20,ha='right',fontsize=10)
    ax.set_ylabel('Aitchison PERMANOVA $R^2$',fontsize=11)
    hi=max(float(np.max(z)) for z in data)
    ax.set_ylim(0,hi*1.18)
    for x,t in [(1,'3 realms × 6 localities'),(5,'3 realms × 4 localities')]:
        ax.text(x,hi*1.09,t,ha='center',fontsize=10)
    ax.axvline(3,color='#dddddd',lw=1)
    tidy(ax);title(ax,'a','Balanced sample selections')


def panel_contrast(ax,raw):
    groups=[('Primary design','Primary: 3 realms × 6 localities'),
            ('Alternative design','Alternative: 3 realms × 4 localities'),
            ('Adjusted for KO annotation total','Annotation-total adjustment'),
            ('Equal-KO resampling','Equal KO counts'),('Detected in >= 10 samples','Occurrence ≥10 profiles'),
            ('Prevalence 90/10','Prevalence 90% / 10%'),('Prevalence 97.5/20','Prevalence 97.5% / 20%'),
            ('Leave-one-country-out','Country exclusions'),('Leave-one-study-out','Study / BioProject exclusions')]
    rows=[]
    for g,label in groups:
        z=raw.loc[raw.analysis==g,'delta_R2'].dropna()
        if not len(z):continue
        deletion=g.startswith('Leave')
        lo,hi=(z.min(),z.max()) if deletion else z.quantile([.05,.95])
        rows.append((label,z.median(),lo,hi,deletion,len(z)))
    for i,(label,m,lo,hi,deletion,n) in enumerate(rows):
        ax.plot([lo,hi],[i,i],color='#555555',lw=1.6)
        ax.scatter(m,i,s=45,marker='D' if deletion else 'o',facecolor='white' if deletion else '#D55E00',edgecolor='#333333',zorder=3)
    ax.axvline(0,color='#666666',ls='--',lw=1)
    ax.set_yticks(range(len(rows)),[x[0] for x in rows],fontsize=10)
    ax.set_ylim(len(rows)-.5,-.5)
    ax.set_xlabel('Shell $R^2$ − max(Prevalent, Cloud) $R^2$',fontsize=10)
    tidy(ax);title(ax,'b','Paired Shell contrasts')
    return rows


def panel_process(fig, heat_ax, effect_ax, process, heat):
    order=['sulfur_respiration','methanogenic_hdr_transfer','nitrate_sulfite_respiration','thiocyanate_degradation','aromatic_phenolic_transformation','environmental_sensing_osmoregulation']
    labels=['Anaerobic sulfur\nrespiration','Methanogenic/Hdr\nelectron transfer','Anaerobic nitrate and\nsulfite respiration','Thiocyanate degradation','Aromatic and phenolic\ntransformation','Environmental sensing\nand osmoregulation']
    realms=['Central Indo-Pacific','Tropical Atlantic','Western Indo-Pacific']
    effects=process[process.design=='3 x 6'].set_index('process_id').loc[order]
    values=heat.pivot(index='process_id',columns='realm',values='balanced_mean_process_z').loc[order,realms].to_numpy()
    cmap=LinearSegmentedColormap.from_list('process_diverging',['#2166AC','#F7F7F7','#B2182B'])
    im=heat_ax.imshow(values,aspect='auto',vmin=-.8,vmax=.8,cmap=cmap)
    for i in range(6):
        for j in range(3):
            heat_ax.text(j,i,f'{values[i,j]:+.2f}',ha='center',va='center',fontsize=10,color='white' if abs(values[i,j])>.6 else '#222222')
    heat_ax.set_yticks(range(6),[t+f' (n = {int(n)})' for t,n in zip(labels,effects.n_KOs)],fontsize=11)
    heat_ax.set_xticks(range(3),['Central\nIndo-Pacific','Tropical\nAtlantic','Western\nIndo-Pacific'],fontsize=11)
    heat_ax.xaxis.tick_top();heat_ax.tick_params(which='both',length=0,pad=7)
    heat_ax.set_xticks(np.arange(-.5,3,1),minor=True);heat_ax.set_yticks(np.arange(-.5,6,1),minor=True)
    heat_ax.grid(which='minor',color='white',lw=1.5)
    for s in heat_ax.spines.values():s.set_visible(False)
    for i,r in enumerate(effects.itertuples()):
        effect_ax.plot([r.R2_q05,r.R2_q95],[i,i],color='#888888',lw=2)
        effect_ax.scatter(r.R2_median,i,color='#D55E00',s=44,zorder=3)
    effect_ax.set_ylim(5.5,-.5);effect_ax.set_xlim(0,.6)
    effect_ax.set_yticks([]);effect_ax.set_xticks([0,.2,.4,.6])
    effect_ax.set_xlabel('Regional $R^2$\nof process score',fontsize=11)
    tidy(effect_ax)
    pos=heat_ax.get_position();fig.text(.025,pos.y1+.04,'c',fontsize=18,fontweight='bold')
    cb=fig.add_axes([pos.x0+.025,pos.y0-.035,pos.width-.05,.010 if fig.get_size_inches()[1]>8 else .025])
    cbar=fig.colorbar(im,cax=cb,orientation='horizontal',ticks=[-.8,0,.8])
    cbar.ax.tick_params(labelsize=9,length=2);cbar.set_label('Balanced mean process score (z)',fontsize=10,labelpad=2)


def draw_world(ax,world):
    from matplotlib.patches import Polygon
    for feature in world['features']:
        geom=feature['geometry']
        polygons=geom['coordinates'] if geom['type']=='MultiPolygon' else [geom['coordinates']]
        for poly in polygons:
            ax.add_patch(Polygon(poly[0],closed=True,facecolor='#F2F2F2',edgecolor='#000000',linewidth=.65,zorder=0))


def panel_map(fig,ax,data,world):
    from matplotlib.colors import Normalize
    z=data[data.mappable].sort_values('n_profiles',ascending=False)
    norm=Normalize(.25,1.15)
    draw_world(ax,world)
    ax.scatter(z.lon,z.lat,s=18+7*np.minimum(z.n_profiles,20),c=z.turnover_point_estimate,cmap='viridis',norm=norm,edgecolors='#222222',linewidths=.7,zorder=3)
    ax.set(xlim=(-120,140),ylim=(-38,50));ax.set_aspect('equal',adjustable='box');ax.axis('off')
    ax.text(-.038,1.04,'d',transform=ax.transAxes,fontsize=18,fontweight='bold')
    ax.text(0,1.06,'Shell compositional turnover across named localities',transform=ax.transAxes,fontsize=12,fontweight='bold')
    for bounds,xlim,ylim,label in [([.025,.025,.19,.36],(-39.3,-38.1),(-13.3,-12.35),'Brazil'),([.76,.025,.22,.36],(108,115),(20.4,23.1),'South China')]:
        inset=ax.inset_axes(bounds)
        draw_world(inset,world)
        q=z[z.lon.between(*xlim)&z.lat.between(*ylim)]
        inset.scatter(q.lon,q.lat,s=.55*(18+7*np.minimum(q.n_profiles,20)),c=q.turnover_point_estimate,cmap='viridis',norm=norm,edgecolors='#222222',linewidths=.5,zorder=3)
        inset.set(xlim=xlim,ylim=ylim,xticks=[],yticks=[])
        inset.set_title(label,fontsize=9,pad=3,bbox={'facecolor':'white','edgecolor':'none','pad':.8})
        for sp in inset.spines.values():sp.set_color('#BBBBBB');sp.set_linewidth(.6)
    pos=ax.get_position()
    cb=fig.add_axes([.15,pos.y0-.035,.29,.010 if fig.get_size_inches()[1]>8 else .025])
    cbar=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap='viridis'),cax=cb,orientation='horizontal',ticks=[.3,.5,.7,.9,1.1])
    cbar.ax.tick_params(labelsize=9,length=2);cbar.set_label('RMS Aitchison distance from the global centroid',fontsize=9,labelpad=2)
    handles=[ax.scatter([],[],s=18+7*n,facecolor='white',edgecolor='#222222',linewidth=.7,label=str(n)) for n in [1,5,10,20]]
    ax.legend(handles=handles,title='Metagenomes (capped at 20)',ncol=4,frameon=False,loc='upper center',bbox_to_anchor=(.77,-.02),fontsize=9,title_fontsize=9,handletextpad=.35,columnspacing=1)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--tables',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True)
    p.add_argument('--map-only',action='store_true',help='Regenerate only panel d and the assembled Figure 5 as PNG files.')
    p.add_argument('--world-map',type=Path,required=True,help='External Natural Earth GeoJSON; see config/maps/SOURCE.json.')
    a=p.parse_args()
    if not a.world_map.is_file():p.error(f'Missing external basemap: {a.world_map}')
    world=json.loads(a.world_map.read_text())
    a.outdir.mkdir(parents=True,exist_ok=True)
    def read(n):return pd.read_csv(a.tables/n,sep='\t')
    cross=read('180_persistence_prevalence_crosswalk.tsv')
    inter=read('180_panel_b_intersection_balanced_R2.tsv')
    prev=read('180_panel_c_prevalence_balanced_R2.tsv')
    raw=read('180_panel_d_robustness_raw.tsv')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42,'ps.fonttype':42,'axes.labelcolor':'black'})
    process=read('181_process_regional_summary.tsv')
    heat=read('181_process_heatmap.tsv')
    locality=read('182_shell_locality_turnover_summary.tsv')
    fig=plt.figure(figsize=(10.5,13.2))
    axa=fig.add_axes([.19,.718,.30,.235]);axb=fig.add_axes([.71,.718,.26,.235])
    panel_a(axa,cross);panel_b(axb,inter)
    ah=fig.add_axes([.34,.388,.33,.25]);ae=fig.add_axes([.74,.388,.23,.25])
    panel_process(fig,ah,ae,process,heat)
    axm=fig.add_axes([.06,.060,.91,.225])
    panel_map(fig,axm,locality,world)
    save(fig,a.outdir,'Figure_5_functional_regional_identity_combined',formats=('png',) if a.map_only else ('png','pdf','tiff'))
    if a.map_only:
        fig=plt.figure(figsize=(10.5,4.8));axm=fig.add_axes([.06,.24,.91,.65])
        panel_map(fig,axm,locality,world)
        save(fig,a.outdir,'Figure_5d',formats=('png',))
        print('Regenerated panel d and assembled Figure 5 PNGs only.',flush=True)
        return
    for letter,fn,data in [('a',panel_a,cross),('b',panel_b,inter)]:
        fig,ax=plt.subplots(figsize=(7,4.6));fig.subplots_adjust(left=.36 if letter=='b' else .25,right=.97,top=.86,bottom=.18)
        fn(ax,data);save(fig,a.outdir,'Figure_5'+letter)
    fig=plt.figure(figsize=(10.5,5))
    ah=fig.add_axes([.34,.23,.33,.64]);ae=fig.add_axes([.74,.23,.23,.64])
    panel_process(fig,ah,ae,process,heat)
    save(fig,a.outdir,'Figure_5c')
    fig=plt.figure(figsize=(10.5,4.8));axm=fig.add_axes([.06,.24,.91,.65])
    panel_map(fig,axm,locality,world)
    save(fig,a.outdir,'Figure_5d')
    fig,axes=plt.subplots(2,1,figsize=(9.5,9.5));fig.subplots_adjust(left=.40,right=.97,bottom=.08,top=.94,hspace=.42)
    panel_prevalence(axes[0],prev);summary=panel_contrast(axes[1],raw)
    save(fig,a.outdir,'Supplementary_functional_regional_sensitivity')
    fig,ax=plt.subplots(figsize=(9.5,8.2));fig.subplots_adjust(left=.39,right=.97,bottom=.09,top=.94)
    colors={'Central Indo-Pacific':'#0072B2','Tropical Atlantic':'#D55E00','Western Indo-Pacific':'#009E73','Tropical Eastern Pacific':'#CC79A7'}
    zz=locality.sort_values(['realm','turnover_point_estimate']).reset_index(drop=True)
    for i,r in zz.iterrows():
        color=colors[r.realm]
        ax.plot([r.bootstrap_q05,r.bootstrap_q95],[i,i],color=color,lw=1.2)
        ax.scatter(r.turnover_point_estimate,i,s=36,marker='o' if r.mappable else 'D',facecolor=color if r.mappable else 'white',edgecolor=color,zorder=3)
    ax.set_yticks(range(len(zz)),[r.locality+f' (n = {r.n_profiles})'+(' *' if not r.mappable else '') for r in zz.itertuples()],fontsize=10)
    ax.set_ylim(len(zz)-.5,-.5);ax.set_xlim(.2,1.2)
    ax.set_xlabel('RMS Aitchison distance from the global centroid')
    ax.set_title('Shell turnover across all 26 named localities',loc='left',fontsize=13,pad=14)
    tidy(ax);ax.tick_params(axis='y',labelsize=10)
    save(fig,a.outdir,'Supplementary_shell_locality_turnover')
    # Supplementary diagnostic: paired effects of annotation-total adjustment.
    primary=prev[prev.design=='3 x 6']
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.2));fig.subplots_adjust(left=.09,right=.98,bottom=.20,top=.85,wspace=.3)
    for i,g in enumerate(COLORS):
        z=primary[primary.group==g]
        axes[0].scatter(z.R2,z.adjusted_unique_R2,s=9,alpha=.25,color=COLORS[g],label=g)
        axes[1].boxplot([z.permdisp_P],positions=[i],widths=.55,showfliers=False,patch_artist=True,
                       boxprops={'facecolor':COLORS[g],'alpha':.5},medianprops={'color':'black'})
    limit=max(primary.R2.max(),primary.adjusted_unique_R2.max())*1.06
    axes[0].plot([0,limit],[0,limit],ls='--',color='#777777',lw=1)
    axes[0].set(xlim=(0,limit),ylim=(0,limit),xlabel='Unadjusted regional $R^2$',ylabel='Adjusted unique regional $R^2$')
    axes[0].legend(frameon=False,fontsize=9,markerscale=2)
    axes[1].axhline(.05,ls='--',color='#777777',lw=1)
    axes[1].set_xticks(range(3),list(COLORS));axes[1].set(ylim=(0,1.02),ylabel='PERMDISP permutation P')
    title(axes[0],'a','Annotation-total sensitivity');title(axes[1],'b','Dispersion diagnostics')
    for ax in axes:tidy(ax)
    save(fig,a.outdir,'Supplementary_functional_regional_diagnostics')
    exports={p.name:{'size_bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in a.outdir.iterdir() if p.suffix in ['.png','.pdf','.tiff']}
    (a.outdir/'figure_validation.json').write_text(json.dumps({'status':'PASS','figure5_panels':4,'dpi':600,'source':'current Figure 5 tables','deletion_intervals':'range across exclusion medians','sampling_intervals':'5th–95th percentiles','supplementary_sensitivity_rows':summary,'panel_c':'six selected process scores','panel_d':'named-locality RMS Aitchison turnover','exports':exports},indent=2)+'\n')
    print('Figure 5 and supplementary figures exported.',flush=True)


if __name__=='__main__':main()
