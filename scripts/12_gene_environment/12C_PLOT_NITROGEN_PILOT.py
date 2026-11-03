#!/usr/bin/env python3
"""Static and standalone interactive nitrogen atlas from 12B outputs."""
import argparse
import json
import html
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, to_hex
from matplotlib.patches import Rectangle, FancyArrowPatch, Circle, Polygon
from matplotlib.cm import ScalarMappable

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'svg.fonttype': 'none',
                     'axes.spines.top': False, 'axes.spines.right': False})
INK = '#173f47'
PALE = '#edf5f4'
CMAP = plt.get_cmap('viridis')


def country_polygons(world):
    for f in world['features']:
        g = f['geometry']
        if g['type'] == 'Polygon':
            yield g['coordinates'][0]
        elif g['type'] == 'MultiPolygon':
            for p in g['coordinates']:
                yield p[0]


def paint(v, norm):
    if pd.isna(v):
        return '#c6ccd2'
    if v == 0:
        return '#ffffff'
    return to_hex(CMAP(norm(np.log10(1 + v))))


def gene_tiles(ax, genes, values, norm, x, y, width=1.08):
    for i, (ko, name) in enumerate(genes):
        xx = x + i * (width + .06)
        v = values.get(ko, np.nan)
        patch = Rectangle((xx, y), width, .6, facecolor=paint(v,norm), edgecolor='#5a7178', lw=.6,
                          hatch='///' if pd.isna(v) else None)
        patch.set_gid('box_' + ko)
        ax.add_patch(patch)
        dark = pd.notna(v) and v > 0 and norm(np.log10(1+v)) < .62
        t = ax.text(xx+width/2, y+.40, name, ha='center', va='center', fontsize=9,
                    color='white' if dark else INK)
        t.set_gid('label_' + ko)
        t = ax.text(xx+width/2, y+.15, 'ND' if pd.isna(v) else f'{v:.1f}', ha='center', va='center',
                    fontsize=8, color='white' if dark else INK)
        t.set_gid('value_' + ko)


def cycle_figure(values, norm, title, path):
    fig, ax = plt.subplots(figsize=(15, 9))
    fig.patch.set_facecolor('white'); ax.set_facecolor('white')
    ax.set_xlim(0,15); ax.set_ylim(0,9); ax.axis('off')
    ax.text(.15,8.7,'Ciclo del nitrógeno · potencial genético',fontsize=22,weight='bold',color=INK)
    ax.text(.15,8.25,title,fontsize=11,color=INK)
    nodes={'N2':(1.8,6.8,'N₂\nNitrógeno'),'NH4':(1.8,3.1,'NH₄⁺\nAmonio'),
           'NO2':(7.2,3.1,'NO₂⁻\nNitrito'),'NO3':(7.2,5.7,'NO₃⁻\nNitrato'),
           'NO':(12.7,3.1,'NO\nÓxido nítrico'),'N2O':(12.7,6.8,'N₂O\nÓxido nitroso')}
    def arrow(a,b,rad=0,uncertain=False):
        axy=np.array(nodes[a][:2]);bxy=np.array(nodes[b][:2]);vec=bxy-axy;vec=vec/np.linalg.norm(vec)
        ax.add_patch(FancyArrowPatch(axy+.57*vec,bxy-.57*vec,
            connectionstyle=f'arc3,rad={rad}',arrowstyle='-' if uncertain else '-|>',
            mutation_scale=17,lw=1.7,color='#a47d34' if uncertain else '#456873',
            linestyle='--' if uncertain else '-'))
    arrow('N2','NH4');arrow('NH4','NO2',uncertain=True);arrow('NO2','NH4',-.42)
    arrow('NO3','NO2',.16);arrow('NO2','NO3',.16,True)
    arrow('NO2','NO');arrow('NO','N2O');arrow('N2O','N2')
    for x,y,label in nodes.values():
        ax.add_patch(Circle((x,y),.57,color=PALE,ec='#9db9bc',lw=1))
        ax.text(x,y,label,ha='center',va='center',fontsize=10,color=INK,weight='bold')
    ax.text(.2,5.65,'Fijación de N₂',color=INK)
    gene_tiles(ax,[('K02588','nifH'),('K02586','nifD'),('K02591','nifK')],values,norm,.2,4.9,.95)
    ax.text(6.0,7.7,'Reducción de N₂O',color=INK)
    gene_tiles(ax,[('K00376','nosZ')],values,norm,6.3,6.95)
    ax.text(11.2,5.4,'NO → N₂O',color=INK)
    gene_tiles(ax,[('K04561','norB'),('K02305','norC')],values,norm,11.2,4.6)
    ax.text(9.0,4.0,'NO₂⁻ → NO · alternativas',color=INK)
    gene_tiles(ax,[('K00368','nirK'),('K15864','nirS')],values,norm,9.0,3.35)
    ax.text(8.0,5.4,'NO₃⁻ → NO₂⁻',color=INK)
    gene_tiles(ax,[('K02567','napA'),('K02568','napB')],values,norm,8.0,4.6)
    ax.text(4.1,5.65,'Nar / Nxr: dirección\nsin resolver',color='#815f23',fontsize=10)
    gene_tiles(ax,[('K00370','narG/nxrA'),('K00371','narH/nxrB')],values,norm,4.1,4.6,1.2)
    ax.text(3.0,3.65,'Amo/Pmo: sustrato sin resolver',color='#815f23',fontsize=10)
    gene_tiles(ax,[('K28504','amoA actual'),('K10944','Amo/Pmo A')],values,norm,3.0,2.65,1.3)
    ax.text(3.6,1.65,'NO₂⁻ → NH₄⁺',color=INK)
    gene_tiles(ax,[('K03385','nrfA'),('K15876','nrfH')],values,norm,3.6,.85)
    ax.text(9.3,1.9,'Anammox: información no disponible',fontsize=10,color=INK,va='top')
    ax.text(9.3,1.5,'hzsABC y hdh carecen de fila en la matriz.\nUrea, cianato y asimilación: panel completo\nen el atlas y el visor.',fontsize=9,color=INK,linespacing=1.5,va='top')
    fig.subplots_adjust(left=.035,right=.985,top=.97,bottom=.13)
    cb=fig.colorbar(ScalarMappable(norm=norm,cmap=CMAP),ax=ax,orientation='horizontal',fraction=.035,pad=.04,shrink=.48)
    ticks=[0,1,2,3,4];ticks=[v for v in ticks if v<=norm.vmax]
    cb.set_ticks(ticks);cb.set_ticklabels([str(10**v-1) for v in ticks])
    cb.set_label('Asignaciones por millón de asignaciones KO · escala log₁₀(1 + valor)',fontsize=9)
    fig.text(.04,.04,'Blanco: cero observado. Gris rayado: KO sin fila (ND). Líneas discontinuas: asignación no resuelta.',fontsize=9,color=INK)
    fig.text(.04,.018,'Las flechas describen transformaciones de referencia; el ADN no mide actividad ni flujo. No se suman las subunidades.',fontsize=9,color=INK)
    for ext in ['png','pdf','svg']:
        fig.savefig(path.with_suffix('.'+ext),dpi=220,facecolor='white')
    plt.close(fig)


def clean_json(value):
    if isinstance(value,dict): return {str(k):clean_json(v) for k,v in value.items()}
    if isinstance(value,list): return [clean_json(v) for v in value]
    if isinstance(value,np.generic): value=value.item()
    if isinstance(value,float) and not np.isfinite(value): return None
    return value


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['results','reactions','world','outdir']:p.add_argument('--'+k,required=True)
    a=p.parse_args();src=Path(a.results);out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
    reactions=json.loads(Path(a.reactions).read_text());world=json.loads(Path(a.world).read_text())
    rel=pd.read_csv(src/'nitrogen_per_million_KO_assignments.tsv',sep='\t',index_col=0)
    groups=pd.read_csv(src/'study_coordinate_groups.tsv',sep='\t',index_col=0)
    samples=pd.read_csv(src/'sample_analysis_manifest.tsv',sep='\t',index_col=0)
    grouped=pd.read_csv(src/'group_relative_representation.tsv',sep='\t',index_col=0)
    inv=pd.read_csv(src/'nitrogen_KO_inventory.tsv',sep='\t').set_index('KO')
    vmax=float(np.ceil(np.log10(1+np.nanmax(rel.to_numpy())))*1.0);norm=Normalize(0,vmax)
    cycle_figure(grouped.median(),norm,f'Mediana de {len(groups)} resúmenes estudio × coordenada · {len(samples)} perfiles',out/'01_Ciclo_nitrogeno')
    order=groups.sort_values(['REALM','Country','lon','lat','study_id']).index
    fig,ax=plt.subplots(figsize=(18,14))
    cmap=CMAP.copy();cmap.set_bad('#c6ccd2')
    values=grouped.loc[order,rel.index].T.to_numpy()
    im=ax.imshow(np.ma.masked_invalid(np.log10(1+values)),aspect='auto',cmap=cmap,norm=norm,interpolation='none')
    # Zero observations remain white and distinct from missing KO rows.
    for row,col in zip(*np.where(values==0)):
        ax.add_patch(Rectangle((col-.5,row-.5),1,1,facecolor='white',edgecolor='none'))
    ax.set_yticks(np.arange(len(rel)));ax.set_yticklabels([f'{inv.loc[k,"gene"]} · {k}' for k in rel.index],fontsize=9)
    ax.set_xticks(np.arange(len(order)));ax.set_xticklabels([f'G{i+1:02d}' for i in range(len(order))],rotation=90,fontsize=6)
    ax.set_title(f'Nitrógeno: representación de los {len(rel)} KOs del panel\n{len(groups)} grupos estudio × coordenada; escala común, sin estandarización por fila',loc='left',fontsize=16,pad=43,color=INK)
    old=None
    for i,g in enumerate(order):
        realm=groups.loc[g,'REALM']
        if realm!=old:
            short={'Central Indo-Pacific':'CIP','Western Indo-Pacific':'WIP','Tropical Atlantic':'TA','Tropical Eastern Pacific':'TEP'}.get(realm,realm)
            ax.axvline(i-.5,color=INK,lw=1.4); ax.text(i,-1.2,short,fontsize=8,weight='bold',ha='right' if short=='TEP' else 'left');old=realm
    fig.colorbar(im,ax=ax,shrink=.45,pad=.02,label='log₁₀(1 + asignaciones por millón)')
    fig.text(.02,.02,'Blanco = cero observado. Gris = KO sin fila. Grupos cartográficos descriptivos; no equivalen a réplicas biológicas verificadas.',fontsize=10)
    fig.tight_layout(rect=(0,.04,1,1))
    for ext in ['png','pdf']:fig.savefig(out/('02_Heatmap_nitrogeno.'+ext),dpi=220)
    plt.close(fig)
    groups.loc[order].assign(heatmap_id=[f'G{i+1:02d}' for i in range(len(order))]).to_csv(src/'heatmap_group_key.tsv',sep='\t')
    # Geographic layers: one point per coordinate, median of study-level summaries.
    coord=groups.groupby('coordinate_group')[['lat','lon']].first()
    cg=grouped.join(groups.coordinate_group).groupby('coordinate_group').median()
    mapkos=['K02588','K02567','K00368','K15864','K00376','K03385']
    fig,axs=plt.subplots(3,2,figsize=(16,10),sharex=True,sharey=True)
    for ax,ko in zip(axs.flat,mapkos):
        ax.set_facecolor('#f0f6fa')
        for poly in country_polygons(world):
            ax.add_patch(Polygon(poly,closed=True,fc='#dfe5e3',ec='white',lw=.25))
        ax.scatter(coord.lon,coord.lat,c=[paint(v,norm) for v in cg.reindex(coord.index)[ko]],s=36,edgecolors=INK,lw=.45,zorder=3)
        ax.set_xlim(-125,155);ax.set_ylim(-40,40);ax.set_aspect('equal');ax.grid(alpha=.15,lw=.4)
        ax.set_title(f'{inv.loc[ko,"gene"]} · {ko}',loc='left',weight='bold',color=INK)
    fig.suptitle('Distribución geográfica de seis marcadores del nitrógeno',fontsize=19,color=INK,x=.06,ha='left')
    fig.subplots_adjust(left=.055,right=.95,top=.92,bottom=.14,hspace=.2,wspace=.12)
    cax=fig.add_axes([.3,.095,.4,.02]);fig.colorbar(ScalarMappable(norm=norm,cmap=CMAP),cax=cax,orientation='horizontal',label='log₁₀(1 + asignaciones por millón) · misma escala en los seis mapas')
    fig.text(.055,.018,f'{len(coord)} coordenadas. Si coinciden estudios, el punto muestra la mediana de sus resúmenes. Cero = blanco. Base: Natural Earth 1:110 millones.',fontsize=9)
    for ext in ['png','pdf']:fig.savefig(out/('03_Mapas_genes_nitrogeno.'+ext),dpi=220)
    plt.close(fig)
    diagnostics=pd.read_csv(src/'environment_association_diagnostics.tsv',sep='\t')
    diagnostics=diagnostics[diagnostics.outcome=='log2_positive_group_mean']
    fig,axs=plt.subplots(1,3,figsize=(14,6),sharey=True)
    selected=['K02588','K02567','K00368','K15864','K00376','K03385']
    for ax,var,title in zip(axs,['bio1','soc_mean','phh2o_mean'],['Temperatura climática','Carbono orgánico cartográfico','pH cartográfico']):
        for i,k in enumerate(selected):
            row=diagnostics[(diagnostics.KO==k)&(diagnostics.predictor==var)]
            if len(row)!=1:continue
            row=row.iloc[0];ax.plot([row.loo_min,row.loo_max],[i,i],color='#78979d',lw=2)
            ax.scatter(row.beta_per_sd,i,s=46,color=INK,zorder=4)
        ax.axvline(0,color='#b37433',ls='--',lw=1)
        ax.set_title(title,fontsize=12,color=INK);ax.set_xlabel('Cambio log₂ por 1 DE del predictor',fontsize=10)
        ax.grid(axis='x',alpha=.2);ax.set_yticks(range(len(selected)));ax.set_yticklabels([inv.loc[k,'gene'] for k in selected]);ax.invert_yaxis()
    fig.suptitle('Diagnóstico de sensibilidad: retirar un estudio puede cambiar el signo',fontsize=17,color=INK)
    fig.subplots_adjust(left=.09,right=.97,top=.82,bottom=.25,wspace=.27)
    fig.text(.06,.10,'Punto: pendiente ajustada por estudio y los otros dos predictores. Línea: mínimo y máximo al retirar estudios con contraste ambiental.',fontsize=10)
    fig.text(.06,.055,'Las líneas NO son intervalos de confianza. Diagnóstico exploratorio de resúmenes cartográficos positivos; sin valores p ni q.',fontsize=10,color='#815f23')
    for ext in ['png','pdf']:fig.savefig(out/('04_Sensibilidad_ambiental.'+ext),dpi=220)
    plt.close(fig)
    # Stable-coordinate map embedded as SVG; no network requests on opening the viewer.
    def xy(lon,lat):return ((lon+180)/360*1000,(85-lat)/145*400)
    paths=[]
    for poly in country_polygons(world):
        pts=[xy(*pt[:2]) for pt in poly]
        paths.append('<path d="M'+' L'.join(f'{x:.2f},{y:.2f}' for x,y in pts)+' Z" fill="#dfe7e5" stroke="white" stroke-width="0.5"/>')
    for g,row in coord.iterrows():
        x,y=xy(row.lon,row.lat)
        paths.append(f'<circle class="point" data-coord="{html.escape(g)}" cx="{x:.3f}" cy="{y:.3f}" r="4.5" stroke="#173f47" stroke-width="1"><title>{html.escape(g)}</title></circle>')
    mapsvg='<svg viewBox="0 0 1000 400" aria-label="Mapa mundial de coordenadas" style="background:#eff6fa">'+''.join(paths)+'</svg>'
    cyclesvg=(out/'01_Ciclo_nitrogeno.svg').read_text();cyclesvg=cyclesvg[cyclesvg.index('<svg'):]
    # Main static title is replaced in HTML by the current selection's title above the plot.
    cyclesvg=re.sub(r'Mediana de .*?perfiles','Selección actual · valores por KO',cyclesvg)
    D=clean_json(dict(reactions=reactions,groups=groups.to_dict('index'),samples=samples[['analysis_group','coordinate_group','study_id','lat','lon','Country','Locality','bio1','soc_mean','phh2o_mean']].to_dict('index'),
        group_values=grouped.to_dict('index'),sample_values=rel.T.to_dict('index'),coordinate_values=cg.to_dict('index'),
        global_values=grouped.median().to_dict(),palette=[to_hex(CMAP(i/255)) for i in range(256)],vmax=vmax))
    template=Path(__file__).with_name('nitrogen_viewer_template.html').read_text()
    for key,value in {'__N_PROFILES__':len(samples),'__N_KOS__':len(rel),'__N_COORDS__':len(coord),'__N_GROUPS__':len(groups),'__N_AVAILABLE__':int(inv.row_available.sum())}.items():
        template=template.replace(key,str(value))
    document=template.replace('__MAP_SVG__',mapsvg).replace('__CYCLE_SVG__',cyclesvg).replace('__DATA__',json.dumps(D,ensure_ascii=False,allow_nan=False).replace('</','<\\/'))
    (out/'Visor_ciclo_nitrogeno_191.html').write_text(document)
    print(json.dumps({'figures':4,'html':True,'n_map_points':len(coord),'n_group_options':len(groups),'n_run_options':len(samples),'color_max_log1p':vmax}))


if __name__=='__main__':main()
