"""Descriptive figures, no inferential geographical comparisons."""
from pathlib import Path
import json,textwrap,argparse
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch
ROOT=Path(__file__).resolve().parents[2]
ap=argparse.ArgumentParser();ap.add_argument('--results',required=True);ap.add_argument('--figures',required=True);args=ap.parse_args()
O=Path(args.results);F=Path(args.figures);F.mkdir(exist_ok=True,parents=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42})
summary=pd.read_csv(O/'resumen_reglas.tsv',sep='\t',keep_default_na=False)
scope=pd.read_csv(O/'alcance_65_procesos.tsv',sep='\t',keep_default_na=False).set_index('process_id')
meta=pd.read_csv(O/'muestras_diseno_y_deteccion.tsv',sep='\t',keep_default_na=False).set_index('run_accession')
diag=json.loads((O/'diagnostico_y_validacion.json').read_text())
n=len(meta)
colors={'Central Indo-Pacific':'#2B8CBE','Western Indo-Pacific':'#C79A28','Tropical Atlantic':'#3B8E70','Tropical Eastern Pacific':'#9B6C9E'}
short={'Central Indo-Pacific':'CIP','Western Indo-Pacific':'WIP','Tropical Atlantic':'TA','Tropical Eastern Pacific':'TEP'}
for realm in sorted(set(meta.REALM)-set(colors)):
    colors[realm]=plt.get_cmap('tab20')(len(colors)%20)
    short[realm]=realm
selected=[('C01_M00165','Calvin: componentes del módulo'),('M01_MCR','McrABG: dirección sin resolver'),('N01_nifHDK','Nitrogenasa NifHDK'),('N01_nifHDKENB','NifHDKENB: criterio adicional'),('N06_nir','Nitrito → NO: NirK o NirS'),('N08_nosZ','N₂O reductasa: NosZ'),('N09_nrfAH','Amonificación de nitrito: NrfAH'),('S01_sat_apr_dsr','Sat–AprAB–DsrAB: dirección sin resolver'),('S04_Sox_core','SoxAXBYZ: núcleo sin SoxCD'),('S05_sqr','Oxidación de sulfuro: Sqr'),('S10_dmdA','DMSP: desmetilación inicial DmdA'),('P01_pstSCAB','Captación de fosfato: PstSCAB'),('P03_phnJ','C–P liasa: marcador PhnJ'),('P03_phnGHIJKLM','C–P liasa: conjunto PhnGHIJKLM')]
fig,(ax,bx)=plt.subplots(1,2,figsize=(16,9),gridspec_kw={'width_ratios':[1.5,1]})
fig.subplots_adjust(left=.285,right=.98,top=.84,bottom=.19,wspace=.27)
for i,(rid,label) in enumerate(selected):
    r=summary[summary.rule_id==rid].set_index('umbral_registros')
    a,b=int(r.loc[1,'n_cumple']),int(r.loc[2,'n_cumple'])
    ax.plot([b,a],[i,i],color='#A7B8C2',lw=3,zorder=1)
    ax.scatter(a,i,s=55,c='#126B81',zorder=3)
    ax.scatter(b,i,s=40,facecolors='white',edgecolors='#C6652D',linewidth=1.6,zorder=4)
    ax.text(n*1.025,i,f'{a} / {b}',va='center',fontsize=9)
ax.set_yticks(range(len(selected)),[x[1] for x in selected]);ax.invert_yaxis();ax.set_xlim(0,n*1.20)
ax.set_xticks(sorted(set([0,int(n*.25),int(n*.5),int(n*.75),n])));ax.grid(axis='x',color='#E4E9EB',zorder=0)
ax.set_xlabel(f'Metagenomas que satisfacen la regla (de {n})',labelpad=12)
ax.set_title('A  Componentes observados',loc='left',fontweight='bold',pad=17)
ax.text(n*1.025,-1.0,'≥1 / ≥2',fontsize=9,fontweight='bold')
ax.scatter([],[],s=55,c='#126B81',label='≥1 registro por KO')
ax.scatter([],[],s=40,facecolors='white',edgecolors='#C6652D',label='≥2 registros por KO')
ha,la=ax.get_legend_handles_labels()
fig.legend(ha,la,frameon=False,loc='lower left',bbox_to_anchor=(.282,.072),fontsize=9,ncol=2)
for realm,g in meta.groupby('REALM'):
    g=g[g.total_asignaciones_KO>0]
    bx.scatter(g.total_asignaciones_KO,g.n_reglas_cumplidas_ge1,s=34,c=colors[realm],alpha=.75,edgecolors='white',linewidth=.3,label=f'{short[realm]} (n={len(g)})')
bx.set_xscale('log');bx.set_xlabel('Total de asignaciones KO\n(escala logarítmica)',labelpad=10)
bx.set_ylabel(f'Reglas satisfechas de las {diag["n_reglas"]} evaluadas')
bx.set_title('B  Dependencia del rendimiento de anotación',loc='left',fontweight='bold',pad=17,fontsize=11)
bx.grid(color='#E4E9EB',zorder=0);bx.set_axisbelow(True)
bx.text(.03,.95,('Spearman ρ = no estimable' if diag['rho_rendimiento_anotacion_numero_reglas'] is None else f'Spearman ρ = {diag["rho_rendimiento_anotacion_numero_reglas"]:.3f}'),transform=bx.transAxes,va='top',bbox=dict(facecolor='white',edgecolor='none',alpha=.8))
hb,lb=bx.get_legend_handles_labels()
fig.legend(hb,lb,frameon=False,fontsize=9,loc='lower left',bbox_to_anchor=(.73,.060),ncol=2)
fig.suptitle(f'Evidencia funcional en {n} metagenomas de sedimentos de manglar',x=.04,ha='left',y=.97,fontsize=20,fontweight='bold')
fig.text(.04,.91,'Reglas de genes y subunidades aplicadas a la matriz actual. Comparación descriptiva de dos umbrales.',fontsize=12,color='#47545D')
fig.text(.04,.018,'Codetección comunitaria: no demuestra ruta completa en un organismo, actividad ni flujo. ≥2 es una sensibilidad exploratoria.\nEl número de reglas incluye conjuntos solapados: sirve como diagnóstico de detección, no como índice de diversidad funcional.',fontsize=10,color='#47545D')
for ext in ['png','pdf','svg']:fig.savefig(F/f'01_resumen_evidencia_{n}.{ext}',dpi=300,facecolor='white')
plt.close(fig)

# Three atlas pages preserve readable labels and all 65 categories.
status=pd.read_csv(O/'matriz_65_reglas_representativas_ge1.tsv',sep='\t',index_col=0,keep_default_na=False)
state_codes={'requiere_otra_evidencia':0,'sin_filas':1,'sin_deteccion_con_vacios':2,'sin_deteccion':3,'componentes_parciales':4,'regla_cumplida':5}
state_colors=['#A9ABB2','#755A94','#DAD0E5','#F3F4F5','#E4AF63','#177E89']
state_labels=['Requiere otra evidencia','Sin filas de los componentes','Sin detección; hay KOs sin fila','Sin detección en las filas evaluadas','Componentes parciales','Regla representativa satisfecha']
order=meta.reset_index().assign(realm_order=lambda x:x.REALM.map({'Central Indo-Pacific':0,'Western Indo-Pacific':1,'Tropical Atlantic':2,'Tropical Eastern Pacific':3})).sort_values(['realm_order','REALM','study_id','sampling_site_id','run_accession']).run_accession.tolist()
pd.DataFrame({'posicion':range(1,len(order)+1),'run_accession':order,'REALM':meta.loc[order,'REALM'].values,'study_id':meta.loc[order,'study_id'].values}).to_csv(O/'orden_muestras_atlas.tsv',sep='\t',index=False)
pages=[('02_atlas_carbono_metano',['C','CH4_C1'],'Carbono y metano'),('03_atlas_nitrogeno_azufre',['N','S'],'Nitrógeno y azufre'),('04_atlas_fosforo_hierro_energia',['P','Fe','H2_electrones','O2_luz'],'Fósforo, hierro y contexto energético')]
for fname,cycles,title in pages:
    pids=scope[scope.ciclo.isin(cycles)].index.tolist()
    arr=status.loc[pids,order].map(state_codes.__getitem__).values.astype(int)
    fig,ax=plt.subplots(figsize=(16,9));fig.subplots_adjust(left=.40,right=.985,bottom=.16,top=.85)
    im=ax.imshow(arr,aspect='auto',interpolation='nearest',cmap=ListedColormap(state_colors),norm=BoundaryNorm(np.arange(-.5,6.5),6))
    labels=[]
    for pid in pids:
        r=scope.loc[pid];lab=r.alcance_evaluado or r.categoria
        tag={'ambiguo':' [A]','contexto':' [C]','no_resuelto':' [NE]'}.get(r.tipo,'')
        labels.append(pid+'  '+lab+tag)
    ax.set_yticks(range(len(pids)),labels,fontsize=9)
    ax.tick_params(axis='y',length=0);ax.set_xticks([])
    ordered=meta.loc[order];start=0
    for realm,g in ordered.groupby('REALM',sort=False):
        end=start+len(g)
        ax.axvline(start-.5,color='white',lw=2)
        ax.text((start+end-1)/2,-1.1,short[realm]+f' ({len(g)})',ha='center',va='bottom',fontsize=10,color=colors[realm],fontweight='bold')
        start=end
    ax.set_xlabel(f'{n} metagenomas ordenados por región, estudio y sitio. Identificadores en orden_muestras_atlas.tsv.',fontsize=9,labelpad=11)
    fig.suptitle(title+' · reglas representativas',x=.035,ha='left',y=.975,fontsize=21,fontweight='bold')
    fig.text(.035,.92,'Umbral ≥1 registro por KO. Cada fila evalúa el alcance indicado; no resume todas las vías de la categoría.',fontsize=11)
    handles=[Patch(facecolor=c,label=l) for c,l in zip(state_colors,state_labels)]
    fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.03,.055),ncol=3,frameon=False,fontsize=9)
    fig.text(.035,.025,'[A] La identidad o dirección del proceso sigue ambigua. [C] Contexto o componente parcial. [NE] Requiere otra evidencia. No se infiere actividad.',fontsize=9,color='#47545D')
    for ext in ['png','pdf']:fig.savefig(F/f'{fname}.{ext}',dpi=300,facecolor='white')
    plt.close(fig)
print('BIOGEOCHEMICAL_PLOTS_PASS: four local figures; PNG/PDF and summary SVG.')
