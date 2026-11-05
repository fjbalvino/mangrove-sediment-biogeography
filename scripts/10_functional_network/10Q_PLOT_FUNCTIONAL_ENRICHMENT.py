#!/usr/bin/env python3
"""Plot observed functions and conditional support without hiding null results."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import textwrap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import FancyBboxPatch, Patch
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('annotation',Path(__file__).with_name('10O_EXTEND_FUNCTIONAL_ANNOTATION.py'))
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)

DISPLAY=[
 ('PATH:ko00970','tRNA aminoacylation'),
 ('BRITE:ko01007','Amino-acid-related enzymes'),
 ('DETAIL:ko03016:9d6a59d62ff7','Prokaryotic aminoacyl-tRNA synthetases'),
 ('BRITE:ko03400','DNA repair and recombination'),
 ('CATEGORY:09105','Amino-acid metabolism'),
 ('CATEGORY:09104','Nucleotide metabolism'),
 ('CATEGORY:09108','Cofactor and vitamin metabolism'),
 ('CATEGORY:09101','Carbohydrate metabolism'),
 ('CATEGORY:09103','Lipid metabolism'),
 ('PATH:ko00670','One-carbon metabolism via folate'),
 ('PATH:ko00190','Oxidative phosphorylation'),
 ('MODULE:M00144','NADH–quinone oxidoreductase'),
 ('BRITE:ko02000','Transporters'),
 ('PATH:ko02010','ABC transporters'),
 ('DETAIL:ko02000:01153ab05944','Branched-chain amino-acid transporter'),
 ('DETAIL:ko02000:d427f5eb7ba9','Peptide/nickel transporter'),
 ('PATH:ko02024','Quorum-sensing map (M06: transport genes)'),
 ('BRITE:ko03110','Chaperones and protein folding'),
 ('CATEGORY:09132','Signal transduction'),
]


def save_figure(fig,path,dpi=300):
    with tempfile.NamedTemporaryFile(dir=path.parent,suffix=path.suffix,delete=False) as tmp:
        name=Path(tmp.name)
    fig.savefig(name,dpi=dpi,facecolor='white',bbox_inches='tight')
    name.replace(path)


def manuscript_figure():
    """Assemble Figure 5 from frozen graph coordinates and 10Q plotted values.

    This is a presentation-only operation: no statistical test, graph layout,
    module assignment, abundance or detection count is recomputed.
    """
    p=argparse.ArgumentParser(description='Replace historical Figure 5 b/c with the expanded functional panels')
    p.add_argument('--figure5-from',type=Path,required=True,help='Directory of the existing 10Q plotted tables')
    p.add_argument('--network-nodes',type=Path,required=True)
    p.add_argument('--network-edges',type=Path,required=True)
    p.add_argument('--outdir',type=Path,required=True)
    args=p.parse_args();out=args.outdir;out.mkdir(parents=True,exist_ok=True)
    base=args.figure5_from
    composition=pd.read_csv(base/'plotted_functional_composition.tsv',sep='\t')
    components=pd.read_csv(base/'plotted_component_context.tsv',sep='\t').fillna('')
    codetection=json.loads((base/'component_codetection.json').read_text())
    profile_count=json.loads((base/'figure_provenance.json').read_text())['profiles']
    nodes=pd.read_csv(args.network_nodes,sep='\t')
    edges=pd.read_csv(args.network_edges,sep='\t')
    if nodes.KO.duplicated().any() or edges.edge_key.duplicated().any():
        raise ValueError('Duplicate graph identifiers')
    members=set(nodes.KO)
    if not (set(edges['from'])|set(edges['to']))<=members:
        raise ValueError('An edge has an endpoint outside the retained network panel')
    modules=sorted(nodes.module.unique());sizes=nodes.groupby('module').size()
    if modules!=sorted(composition.module.unique()):raise ValueError('Module identities differ between panels')
    if composition.duplicated(['term_id','module']).any():raise ValueError('Duplicate heatmap cells')
    for row in composition.itertuples():
        if not np.isclose(row.percent_module,100*row.n_KOs/sizes[row.module]):
            raise ValueError('Module sizes differ between graph and heatmap')
    palette=dict(zip(modules,['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9']))
    if len(modules)>len(palette):raise ValueError('More modules than the declared manuscript palette')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10.5,'axes.titleweight':'bold','pdf.fonttype':42})
    fig=plt.figure(figsize=(12.8,11.8));dark='#172238';muted='#536174'
    fig.text(.035,.974,'a',fontsize=17,weight='bold',color=dark)
    fig.text(.061,.974,'Co-abundance network',fontsize=13,weight='bold',color=dark)
    fig.text(.363,.974,'b',fontsize=17,weight='bold',color=dark)
    fig.text(.389,.974,'Functional composition and enrichment',fontsize=13,weight='bold',color=dark)
    ax=fig.add_axes([.026,.435,.285,.510])
    ax.add_collection(LineCollection(edges[['x','y','xend','yend']].to_numpy().reshape(-1,2,2),
                                    colors='#7A7A7A',linewidths=.22,alpha=.28))
    shapes={'Absolute global core':'s','Near-global':'o','Realm-core':'^','Peripheral':'D'}
    for tier,marker in shapes.items():
        group=nodes[nodes.tier.eq(tier)]
        ax.scatter(group.x,group.y,c=[palette[m] for m in group.module],marker=marker,
                   s=12,edgecolors='#222222',linewidths=.22,alpha=.96,zorder=3)
    dx=[.05,.10,-.30,.38,.25,-.08];dy=[.30,-.35,.17,.38,-.22,-.23]
    xr=nodes.x.max()-nodes.x.min();yr=nodes.y.max()-nodes.y.min()
    for i,module in enumerate(modules):
        g=nodes[nodes.module.eq(module)];x,y=g[['x','y']].median()
        lx,ly=x+dx[i]*xr,y+dy[i]*yr
        ax.plot([x,lx],[y,ly],color=palette[module],linewidth=.55,zorder=2)
        ax.text(lx,ly,f'{module}\n{len(g)} KOs',ha='center',va='center',fontsize=10.5,weight='bold',
                bbox=dict(facecolor='white',edgecolor='none',pad=2),zorder=4,color=dark)
    ax.set_aspect('equal');ax.autoscale_view();ax.margins(x=.20,y=.12);ax.axis('off')
    handles=[Line2D([],[],marker=mark,color='none',markerfacecolor='white',markeredgecolor='#5b6371',
                    markeredgewidth=.7,markersize=5,label=tier) for tier,mark in shapes.items() if tier in set(nodes.tier)]
    ax.legend(handles=handles,loc='upper left',bbox_to_anchor=(.05,-.015),frameon=False,
              title='Persistence class',fontsize=9,title_fontsize=9.5,handletextpad=.5,borderaxespad=0)
    terms=composition[['term_id','label']].drop_duplicates()
    counts=composition.pivot(index='term_id',columns='module',values='n_KOs').reindex(index=terms.term_id,columns=modules)
    percentages=composition.pivot(index='term_id',columns='module',values='percent_module').reindex(index=terms.term_id,columns=modules)
    stars=composition.pivot(index='term_id',columns='module',values='supported').reindex(index=terms.term_id,columns=modules)
    if counts.isna().any().any():raise ValueError('Incomplete plotted function-by-module grid')
    bx=fig.add_axes([.655,.445,.31,.484])
    cmap=plt.get_cmap('plasma').copy();cmap.set_bad('#f0f2f6')
    # Vector cells preserve sharp artwork at publication size and in PDF editors.
    bx.set_facecolor('#f0f2f6')
    im=bx.pcolormesh(np.arange(len(modules)+1)-.5,np.arange(len(terms)+1)-.5,
        np.ma.masked_where(counts.to_numpy()==0,percentages),cmap=cmap,vmin=0,vmax=50,rasterized=False)
    bx.set_xlim(-.5,len(modules)-.5);bx.set_ylim(len(terms)-.5,-.5)
    labels=['\n'.join(textwrap.wrap(label,37,break_long_words=False,break_on_hyphens=False)) for label in terms.label]
    bx.set_yticks(range(len(terms)),labels,fontsize=10)
    for tick in bx.get_yticklabels():tick.set_linespacing(.95)
    bx.set_xticks(range(len(modules)),[f'{m}\n{sizes[m]} KOs' for m in modules],fontsize=10)
    bx.xaxis.tick_top();bx.tick_params(axis='both',length=0,pad=6)
    bx.set_xticks(np.arange(-.5,len(modules)),minor=True);bx.set_yticks(np.arange(-.5,len(terms)),minor=True)
    bx.grid(which='minor',color='white',linewidth=1.2);bx.tick_params(which='minor',length=0)
    for spine in bx.spines.values():spine.set_visible(False)
    for i in range(len(terms)):
        for j in range(len(modules)):
            n=int(counts.iloc[i,j]);pct=percentages.iloc[i,j]
            bx.text(j,i,str(n)+('*' if stars.iloc[i,j] else ''),ha='center',va='center',fontsize=10.5,
                weight='bold' if stars.iloc[i,j] else 'normal',color='#8390a2' if n==0 else (dark if pct>=32 else 'white'))
    cb=fig.colorbar(im,cax=fig.add_axes([.67,.407,.28,.011]),orientation='horizontal')
    cb.solids.set_rasterized(False)
    cb.solids.set_edgecolor('face');cb.solids.set_linewidth(.25)
    cb.set_ticks([0,10,20,30,40,50]);cb.ax.tick_params(labelsize=9,length=2)
    cb.set_label('KOs assigned to the function (% of module)',fontsize=9.5,labelpad=4)
    fig.text(.035,.337,'c',fontsize=17,weight='bold',color=dark)
    fig.text(.061,.337,'Respiratory machinery',fontsize=13,weight='bold',color=dark)
    fig.text(.56,.337,'d',fontsize=17,weight='bold',color=dark)
    fig.text(.586,.337,'Nutrient transport',fontsize=13,weight='bold',color=dark)
    cx=fig.add_axes([.041,.075,.45,.24]);cx.set_xlim(0,4);cx.set_ylim(0,4.9);cx.axis('off')
    cx.text(0,4.7,'NADH-quinone oxidoreductase | M03',fontsize=10.5,color=dark)
    respiratory=components[components.system.eq('NADH_quinone_oxidoreductase')]
    for i,r in enumerate(respiratory.itertuples()):
        rr,col=divmod(i,4);y=3.50-rr*.91;highlight=r.highlighted_module=='M03'
        cx.add_patch(FancyBboxPatch((col+.03,y),.91,.79,boxstyle='round,pad=.02,rounding_size=.05',
                    facecolor=palette['M03'] if highlight else '#e7eaf0',edgecolor='white'))
        color='white' if highlight else dark
        for yy,txt,size,weight in [(.60,r.label,10.5,'bold'),(.33,r.KO,9,'normal'),(.08,f'{r.profiles_detected}/{profile_count}',8.8,'normal')]:
            cx.text(col+.485,y+yy,txt,ha='center',fontsize=size,weight=weight,color=color)
    fig.text(.041,.061,'B/C/D alternatives: B + C + D, B + CD, or BCD.',fontsize=9,color=muted)
    fig.text(.041,.043,f"Full KEGG rule met in {codetection['1']['NADH_rule_met']}/{profile_count} community profiles.",fontsize=9,color=dark)
    tx=fig.add_axes([.566,.075,.40,.24]);tx.set_xlim(0,3);tx.set_ylim(0,4.9);tx.axis('off')
    tx.text(0,4.7,'Branched-chain amino acids | M06',fontsize=10.5,color=dark)
    indexed=components.set_index('KO')
    for col,(title,kos) in enumerate([('Binding',['K01999']),('Membrane',['K01997','K01998']),('ATP binding',['K01996','K01995'])]):
        tx.text(col+.46,3.9,title,ha='center',fontsize=10,weight='bold',color=dark)
        for rr,ko in enumerate(kos):
            y=2.40-rr*1.22;r=indexed.loc[ko]
            tx.add_patch(FancyBboxPatch((col+.035,y),.87,1.02,boxstyle='round,pad=.02,rounding_size=.05',
                        facecolor=palette['M06'] if r.highlighted_module=='M06' else '#e7eaf0',edgecolor='white'))
            for yy,txt,size,weight in [(.73,r.label,11,'bold'),(.43,ko,9.5,'normal'),(.15,f'{r.profiles_detected}/{profile_count}',9,'normal')]:
                tx.text(col+.47,y+yy,txt,ha='center',fontsize=size,weight=weight,color=dark)
    fig.text(.566,.061,f"Five components co-detected in {codetection['1']['five_Liv_components_codetected']}/{profile_count} profiles.",fontsize=9,color=dark)
    fig.text(.566,.043,'Box colors match network modules; gray: outside that module.',fontsize=9,color=muted)
    # Detailed statistical and biological interpretation belongs in the caption.
    fig.text(.035,.013,'* Primary conditional enrichment supported. Box fractions: profiles with ≥1 annotation record. Community co-detection does not establish same-organism completion.',fontsize=8.4,color=muted)
    for ext in ['png','pdf']:save_figure(fig,out/('Figure_5_functional_coabundance.'+ext),dpi=600)
    plt.close(fig)
    inputs=[args.network_nodes,args.network_edges,*[base/f for f in ['plotted_functional_composition.tsv',
            'plotted_component_context.tsv','component_codetection.json','figure_provenance.json']]]
    a.write_json(out/'Figure_5_provenance.json',dict(
        inputs={str(path):a.digest(path) for path in inputs},source={str(Path(__file__).relative_to(ROOT)):a.digest(__file__)},
        profiles=profile_count,nodes=len(nodes),edges=len(edges),module_sizes={str(k):int(v) for k,v in sizes.items()},
        heatmap_cells=int(counts.size),component_boxes=len(components),
        panel_mapping={'a':'Previous network a; graph geometry and persistence unchanged; module IDs and sizes retained',
            'b':'New enrichment panel a','c':'New enrichment panel b','d':'New enrichment panel c'},
        presentation_changes=['Remove historical function labels from network','Match component highlights to network module colors'],
        statistical_analysis_rerun=False))
    print(json.dumps(dict(files=['Figure_5_functional_coabundance.png','Figure_5_functional_coabundance.pdf'],
                         nodes=len(nodes),edges=len(edges),heatmap_cells=int(counts.size))))


def main():
    if '--figure5-from' in sys.argv:
        manuscript_figure();return
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['enrichment','annotation','matrix','samples','outdir']:
        p.add_argument('--'+key,type=Path,required=True)
    args=p.parse_args();out=args.outdir;out.mkdir(parents=True,exist_ok=True)
    results=pd.read_csv(args.enrichment/'enrichment_with_sensitivity.tsv',sep='\t')
    inventory=pd.read_csv(args.annotation/'KO_annotation_inventory.tsv.gz',sep='\t').fillna('').set_index('KO')
    catalogue=pd.read_csv(args.enrichment/'tested_term_catalogue.tsv.gz',sep='\t').set_index('term_id')
    modules=sorted(results.statistical_module.unique())
    # The displayed panel is a declared biological summary, not the testing
    # family. All 900 tests remain in the result table.
    display=[(key,label) for key,label in DISPLAY if key in set(results.term_id)]
    represented={key for key,_ in display}
    # Include every supported gene set; identical membership can share a row.
    equivalences=[]
    for row in results[results.supported].itertuples():
        if row.term_id in represented:continue
        target=set(str(catalogue.loc[row.term_id,'KOs']).split(';')) & set(inventory.index[inventory.statistical_module!=''])
        alias=next((key for key,_ in display if
            (set(str(catalogue.loc[key,'KOs']).split(';')) & set(inventory.index[inventory.statistical_module!='']))==target),None)
        if alias:
            equivalences.append(dict(term_id=row.term_id,displayed_as=alias,reason='Identical KO membership in primary universe'))
        else:
            display.append((row.term_id,row.term));represented.add(row.term_id)
    matrix=a.read_matrix(args.matrix,args.samples)
    expected=json.loads((args.annotation/'provenance.json').read_text())['inputs']['matrix']['sha256']
    if a.digest(args.matrix)!=expected:raise ValueError('Component matrix differs from annotation source')
    if not matrix.gt(0).sum(axis=1).reindex(inventory.index).eq(inventory.detected_profiles).all():
        raise ValueError('Sample selection differs from the annotated inventory')
    definitions=pd.read_csv(args.annotation/'KEGG_module_definitions.tsv',sep='\t').set_index('module_id')
    respiratory=str(definitions.loc['M00144','definition'])
    liv={'K01999':'LivK','K01997':'LivH','K01998':'LivM','K01996':'LivF','K01995':'LivG'}
    co={}
    for threshold in [1,2]:
        values={k:matrix.loc[k].ge(threshold).to_numpy(np.int8) for k in matrix.index}
        co[str(threshold)]={'NADH_rule_met':int((a.logic.evaluate(a.logic.parse(respiratory),values,len(matrix.columns))==1).sum()),
            'five_Liv_components_codetected':int((a.logic.evaluate(a.logic.parse('+'.join(liv)),values,len(matrix.columns))==1).sum())}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titleweight':'bold','pdf.fonttype':42})
    fig=plt.figure(figsize=(15.6,11.3))
    fig.text(.035,.976,'Functional organization of the co-abundance network',fontsize=20,weight='bold',color='#172238')
    fig.text(.035,.947,'191 metagenomes  ·  six major modules  ·  observed functions and conditional enrichment',fontsize=11,color='#536174')
    ax=fig.add_axes([.27,.17,.35,.71])
    counts=np.zeros((len(display),len(modules)));percent=counts.copy();stars=np.zeros_like(counts,dtype=bool)
    plotted=[]
    for i,(term,label) in enumerate(display):
        for j,module in enumerate(modules):
            r=results[(results.term_id==term)&(results.statistical_module==module)].iloc[0]
            counts[i,j]=r.observed_KOs;percent[i,j]=r.percent_module;stars[i,j]=r.supported
            plotted.append(dict(term_id=term,label=label,module=module,n_KOs=int(r.observed_KOs),
                percent_module=float(r.percent_module),supported=bool(r.supported),q_BY=float(r.q_BY),
                p_maxT=float(r.p_maxT_joint_null),supported_sensitivity=bool(r.supported_sensitivity)))
    norm=Normalize(0,50);cmap=plt.get_cmap('plasma').copy();cmap.set_bad('#f0f2f6')
    im=ax.imshow(np.ma.masked_where(counts==0,percent),aspect='auto',cmap=cmap,norm=norm)
    ax.set_yticks(range(len(display)),[label for _,label in display],fontsize=10)
    sizes=results.groupby('statistical_module').module_KOs.first()
    ax.set_xticks(range(len(modules)),[f'{m}\n{int(sizes[m])} KOs' for m in modules],fontsize=10)
    ax.xaxis.tick_top();ax.tick_params(axis='both',length=0,pad=8)
    ax.set_xticks(np.arange(-.5,len(modules),1),minor=True)
    ax.set_yticks(np.arange(-.5,len(display),1),minor=True)
    ax.grid(which='minor',color='white',linewidth=1.3);ax.tick_params(which='minor',length=0)
    for spine in ax.spines.values():spine.set_visible(False)
    for i in range(len(display)):
        for j in range(len(modules)):
            value=int(counts[i,j]);color='#8390a2' if value==0 else ('#18223b' if percent[i,j]>=32 else 'white')
            ax.text(j,i,str(value)+('*' if stars[i,j] else ''),ha='center',va='center',fontsize=10,
                    weight='bold' if stars[i,j] else 'normal',color=color)
    fig.text(.035,.905,'a',weight='bold',fontsize=16)
    fig.text(.057,.905,'Functional composition',weight='bold',fontsize=12)
    cbax=fig.add_axes([.30,.112,.28,.014]);cb=fig.colorbar(im,cax=cbax,orientation='horizontal')
    cb.set_label('KOs assigned to the function (% of module)',fontsize=10)
    cb.set_ticks([0,10,20,30,40,50])
    fig.text(.035,.055,'Numbers: distinct KOs in each module.  * BY-adjusted P ≤ 0.05 and joint-null maxT P ≤ 0.05.',fontsize=9.5,color='#354359')
    fig.text(.035,.032,'Terms overlap. All tests use the complete 150-term × 6-module family; unstarred functions remain descriptive.',fontsize=9.5,color='#354359')
    purple='#92278f';orange='#ef9b28';gray='#e7eaf0';dark='#213047'
    bx=fig.add_axes([.68,.48,.29,.40]);bx.set_xlim(0,4);bx.set_ylim(0,5.25);bx.axis('off')
    bx.text(-.08,5.2,'b',weight='bold',fontsize=16)
    bx.text(.15,5.2,'Respiratory machinery',weight='bold',fontsize=12)
    bx.text(0,4.8,'NADH–quinone oxidoreductase · M03',fontsize=10,color=dark)
    genes=[('K%05d'%i,'Nuo'+chr(65+i-330)) for i in range(330,344)]+[('K13378','NuoCD'),('K13380','NuoBCD')]
    components=[]
    for i,(ko,label) in enumerate(genes):
        row,col=divmod(i,4);x=col;y=3.72-row*.79
        highlighted=inventory.loc[ko,'statistical_module']=='M03'
        bx.add_patch(FancyBboxPatch((x+.035,y),.91,.68,boxstyle='round,pad=.02,rounding_size=.05',
            facecolor=purple if highlighted else gray,edgecolor='white',linewidth=1))
        color='white' if highlighted else dark
        bx.text(x+.49,y+.47,label,ha='center',fontsize=10,weight='bold',color=color)
        bx.text(x+.49,y+.29,ko,ha='center',fontsize=8.5,color=color)
        count=int(inventory.loc[ko,'detected_profiles'])
        bx.text(x+.49,y+.11,f'{count}/191',ha='center',fontsize=8,color=color)
        components.append(dict(system='NADH_quinone_oxidoreductase',KO=ko,label=label,
            highlighted_module='M03' if highlighted else '',profiles_detected=count))
    bx.text(0,.37,'B/C/D alternatives: B + C + D, B + CD, or BCD.',fontsize=8.7,color=dark)
    bx.text(0,.05,f"Full KEGG rule met in {co['1']['NADH_rule_met']}/191 community profiles.",fontsize=9,color=dark)
    cx=fig.add_axes([.68,.17,.29,.26]);cx.set_xlim(0,3);cx.set_ylim(0,3.6);cx.axis('off')
    cx.text(-.06,3.5,'c',weight='bold',fontsize=16)
    cx.text(.13,3.5,'Nutrient transport',weight='bold',fontsize=12)
    cx.text(0,3.12,'Branched-chain amino acids · M06',fontsize=10,color=dark)
    for col,(title,ks) in enumerate([('Binding',['K01999']),('Membrane',['K01997','K01998']),('ATP binding',['K01996','K01995'])]):
        cx.text(col+.46,2.68,title,ha='center',fontsize=9,weight='bold',color=dark)
        for row,ko in enumerate(ks):
            y=1.72-row*.8
            highlight=inventory.loc[ko,'statistical_module']=='M06'
            cx.add_patch(FancyBboxPatch((col+.025,y),.88,.66,boxstyle='round,pad=.02,rounding_size=.05',
                facecolor=orange if highlight else gray,edgecolor='white'))
            count=int(inventory.loc[ko,'detected_profiles'])
            cx.text(col+.46,y+.45,liv[ko],ha='center',fontsize=10,weight='bold',color=dark)
            cx.text(col+.46,y+.26,ko,ha='center',fontsize=8.5,color=dark)
            cx.text(col+.46,y+.08,f'{count}/191',ha='center',fontsize=8,color=dark)
            components.append(dict(system='branched_chain_amino_acid_transport',KO=ko,label=liv[ko],
                highlighted_module='M06' if highlight else '',profiles_detected=count))
    cx.text(0,.44,f"Five components co-detected in {co['1']['five_Liv_components_codetected']}/191 profiles.",fontsize=9,color=dark)
    cx.text(0,.05,'Gene-set context; same-organism completion is not established.',fontsize=8.5,color='#536174')
    fig.text(.68,.121,'Colored: in the highlighted network module.',fontsize=9,color=dark)
    fig.text(.68,.101,'Gray: detected in the matrix, outside that module.',fontsize=9,color=dark)
    fig.text(.68,.079,'Fractions in boxes: profiles with ≥1 annotation record.',fontsize=9,color=dark)
    for suffix in ['png','pdf']:save_figure(fig,out/('06_Functional_module_enrichment.'+suffix))
    plt.close(fig)
    a.write_tsv(pd.DataFrame(plotted),out/'plotted_functional_composition.tsv')
    a.write_tsv(pd.DataFrame(components),out/'plotted_component_context.tsv')
    a.write_json(out/'component_codetection.json',co)
    a.write_json(out/'display_aliases.json',equivalences)
    sources=[Path(__file__),Path(__file__).with_name('10O_EXTEND_FUNCTIONAL_ANNOTATION.py'),
        Path(__file__).parent/'lib/biogeochemical_context.py',a.logic.__file__]
    a.write_json(out/'figure_provenance.json',dict(
        inputs={str(path):a.digest(path) for path in [args.matrix,args.samples,
            args.enrichment/'enrichment_with_sensitivity.tsv',args.annotation/'KEGG_module_definitions.tsv',
            args.annotation/'KO_annotation_inventory.tsv.gz',args.enrichment/'tested_term_catalogue.tsv.gz']},
        source={str(Path(path).relative_to(ROOT)):a.digest(path) for path in sources},
        versions=dict(matplotlib=matplotlib.__version__,numpy=np.__version__,pandas=pd.__version__),
        profiles=len(matrix.columns),support_symbols='Primary comparison only; see results for background sensitivity'))
    print(json.dumps(dict(displayed_terms=len(display),components=co,files=['06_Functional_module_enrichment.png','06_Functional_module_enrichment.pdf'])))


if __name__=='__main__':main()
