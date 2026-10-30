#!/usr/bin/env python3
from pathlib import Path
import json,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.ticker import FuncFormatter,MaxNLocator
from PIL import Image
import argparse
parser=argparse.ArgumentParser(description='Draw checked figure-derived sampling positions with a provisional 1 km radius.')
parser.add_argument('--evidence-dir',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
D=args.evidence_dir;r=json.loads((D/'georeferencing_evidence.json').read_text())
def x(lon):return math.radians(lon)*6378137
def y(lat):return math.log(math.tan(math.pi/4+math.radians(lat)/2))*6378137
fig,axes=plt.subplots(2,1,figsize=(12,11));fig.subplots_adjust(top=.90,bottom=.12,hspace=.28,left=.10,right=.97)
fig.suptitle('Ubicación aproximada de las seis muestras',x=.10,y=.97,ha='left',fontsize=21,fontweight='bold')
fig.text(.10,.935,'Cinco marcadores recuperados del mapa de Jin et al. (2019)',fontsize=13,color='#334155')
for ax,panel,name,bbox in zip(axes,['B','C'],['beihai_wide','sanya_reference'],[(109.073,109.167,21.551,21.62),(109.479,109.725,18.204,18.29)]):
 e=r['panels'][panel]['reference_export']['extent'];im=Image.open(D/(name+'.png'));ax.imshow(im,extent=[e['xmin'],e['xmax'],e['ymin'],e['ymax']],origin='upper')
 ax.set_xlim(x(bbox[0]),x(bbox[1]));ax.set_ylim(y(bbox[2]),y(bbox[3]));ax.set_aspect('equal')
 ax.set_title('Guangxi · desembocadura del Nanliu / Beibu Bay' if panel=='B' else 'Hainan · Sanya',loc='left',fontsize=14,fontweight='bold',pad=10)
 for site in [s for s in r['sites'] if s['panel']==panel]:
  lat,lon=site['lat_unrounded'],site['lon_unrounded'];xx,yy=x(lon),y(lat);rad=1000/math.cos(math.radians(lat))
  ax.add_patch(Circle((xx,yy),rad,edgecolor='white',facecolor='#fbbf24',alpha=.30,linestyle='--',linewidth=1.5))
  ax.plot(xx,yy,'o',ms=7,color='#e11d48',mec='white',mew=1.2)
  label=site['label'].replace(' / ','\n');offset=(0,26)
  if label=='HN_17_bay':offset=(0,-35)
  ax.annotate(label,(xx,yy),xytext=offset,textcoords='offset points',ha='center',va='center',fontsize=10,color='white',bbox=dict(boxstyle='round,pad=.3',fc='#0f172a',alpha=.82,ec='none'))
 ax.xaxis.set_major_locator(MaxNLocator(5));ax.yaxis.set_major_locator(MaxNLocator(4))
 ax.xaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{math.degrees(v/6378137):.3f}°E'))
 ax.yaxis.set_major_formatter(FuncFormatter(lambda v,p:f'{math.degrees(2*math.atan(math.exp(v/6378137))-math.pi/2):.3f}°N'))
 ax.tick_params(labelsize=9);ax.grid(alpha=.25,color='white',linestyle=':')
fig.text(.10,.075,'Punto: centro del marcador del artículo. Círculo: radio operativo provisional de 1 km.',fontsize=11)
fig.text(.10,.052,'Las posiciones no son GPS de campo. El radio no es un intervalo de confianza ni un límite garantizado.',fontsize=10,color='#334155')
fig.text(.10,.028,'Fuente de muestreo: doi:10.1186/s40168-019-0675-9, Fig. 1. Imagen: Esri, Vantor, Earthstar Geographics y GIS User Community.',fontsize=8,color='#475569')
fig.savefig(args.output,dpi=180,facecolor='white')
