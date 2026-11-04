'use strict';const D=__DATA__;let geoMap=null;const liveMarkers=new Map();
const $=id=>document.getElementById(id),fmt=v=>v==null?'Unavailable':Number(v).toLocaleString('en-US',{maximumFractionDigits:2});
function color(v){if(v==null)return '#c6ccd2';if(v===0)return '#ffffff';return D.palette[Math.min(255,Math.max(0,Math.round(Math.log10(1+v)/D.vmax*255)))];}
function foreground(v){return v!=null&&v>0&&Math.log10(1+v)/D.vmax<.54?'white':'#202636';}
function option(parent,value,text){const o=document.createElement('option');o.value=value;o.textContent=text;parent.appendChild(o);}
option($('selection'),'global','All groups · median of __N_GROUPS__ study × coordinate summaries');
const gg=document.createElement('optgroup');gg.label='Study × coordinate groups';$('selection').appendChild(gg);
Object.entries(D.groups).sort((a,b)=>String(a[1].Country).localeCompare(String(b[1].Country))||a[1].lon-b[1].lon).forEach(([id,g])=>option(gg,'group:'+id,`${g.Country} · ${g.Locality} · ${g.study_id} · ${g.lat.toFixed(5)}, ${g.lon.toFixed(5)}`));
const ss=document.createElement('optgroup');ss.label='__N_PROFILES__ individual metagenomes';$('selection').appendChild(ss);Object.entries(D.samples).forEach(([id,s])=>option(ss,'sample:'+id,`${id} · ${s.Country} · ${s.study_id}`));
const genes=D.genes;Object.entries(genes).forEach(([k,g])=>option($('mapgene'),k,`${g} · ${k}`));$('mapgene').value=D.first_gene;
$('scale-description').textContent=`Plasma · log₁₀(1 + value), 0–${D.vmax} · fixed across all eight panels`;

function ruleText(rule,c){const th=$('threshold').value;const data=c.kind==='sample'?D.rule_samples[th][c.id][rule.id]:c.kind==='group'?D.rule_groups[th][c.id][rule.id]:D.rule_global[th][rule.id];if(c.kind==='sample')return ['Observed: rule not met','Required components co-detected','Indeterminate: a required row is unavailable'][data];return `Weighted sample evidence: ${(100*data[1]).toFixed(1)}% met · ${(100*data[0]).toFixed(1)}% not met · ${(100*data[2]).toFixed(1)}% indeterminate`;}
function updateDiagnostic(){const k=$('mapgene').value;const el=$('selected-diagnostic');el.replaceChildren();const h=document.createElement('h3');h.textContent=`Selected gene · ${genes[k]} · ${k}`;el.appendChild(h);const rows=D.diagnostics.filter(r=>r.KO===k);if(!rows.length){const p=document.createElement('p');p.textContent='No estimable positive-response diagnostic for this KO. See model-support data; this does not establish no environmental association.';el.appendChild(p);return;}const table=document.createElement('table');table.className='trace-table';const head=table.createTHead().insertRow();['Predictor','Slope / SD','Study-omission range','After recovery adjustment','Positive groups'].forEach(v=>{const th=document.createElement('th');th.textContent=v;head.appendChild(th);});const body=table.createTBody();rows.forEach(r=>{const tr=body.insertRow();[{'bio1':'Temperature','soc_mean':'Organic carbon','phh2o_mean':'pH'}[r.predictor],fmt(r.beta_per_sd),fmt(r.loo_min)+' to '+fmt(r.loo_max),fmt(r.beta_with_annotation_yield_sensitivity),r.n_groups].forEach(v=>{tr.insertCell().textContent=v;});});el.appendChild(table);appendConditionalSupport(el,rows);}

// Display audited support without changing the original estimates or their colour scale.
function appendConditionalSupport(parent,rows){
 const details=document.createElement('details');details.className='conditional-support';
 const summary=document.createElement('summary');summary.textContent='Where does the environmental contrast come from?';details.appendChild(summary);
 const p=document.createElement('p');p.textContent='Calculated on positive-response groups after removing study means and adjusting each predictor for the other covariates. A large study share indicates concentrated support; these are design diagnostics, not tests of significance.';details.appendChild(p);
 const table=document.createElement('table');table.className='trace-table';const head=table.createTHead().insertRow();
 ['Predictor','Studies contributing','Largest study share','Variation after adjustment','Information-equivalent study count'].forEach(v=>{const th=document.createElement('th');th.textContent=v;head.appendChild(th);});
 const body=table.createTBody(),pct=v=>v==null?'Unavailable':(100*v).toFixed(1)+'%';
 rows.forEach(r=>{const tr=body.insertRow();[{'bio1':'Temperature','soc_mean':'Organic carbon','phh2o_mean':'pH'}[r.predictor],fmt(r.n_studies_with_conditional_information),pct(r.max_conditional_study_information_share),pct(r.conditional_variation_fraction),fmt(r.effective_information_studies)].forEach(v=>tr.insertCell().textContent=v);});
 const note=document.createElement('p');note.className='subtle';note.textContent='Information-equivalent study count = 1 / Σ(study information share²). It describes concentration, not the number of independent replicates or model degrees of freedom. Variation after adjustment is the residual fraction of within-study predictor variation.';
 details.append(table,note);parent.appendChild(details);
}
function appendBlockEvidence(parent,rule,c){
 const data=D.block_evidence[rule.id],th=$('threshold').value;
 if(!data)return;
 const raw=c.kind==='sample'?data.samples[th][c.id]:c.kind==='group'?data.groups[th][c.id]:data.global[th];
 const states=c.kind==='sample'?raw.map(s=>[0,1,2].map(j=>s===j?1:0)):raw;
 const lower=states.reduce((sum,s)=>sum+s[1],0)/states.length;
 const upper=states.reduce((sum,s)=>sum+s[1]+s[2],0)/states.length;
 const coverage=document.createElement('p');coverage.className='block-coverage';
 coverage.textContent=`${c.kind==='sample'?'Logical block coverage':'Mean logical block coverage'}: ${(100*lower).toFixed(1)}–${(100*upper).toFixed(1)}% · ${states.length} block(s).`;
 const details=document.createElement('details');details.className='block-details';const summary=document.createElement('summary');summary.textContent='Inspect required components and unknowns';details.appendChild(summary);
 const table=document.createElement('table');table.className='trace-table';const head=table.createTHead().insertRow();['Block / required components','Evidence'].forEach(v=>{const th=document.createElement('th');th.textContent=v;head.appendChild(th);});
 const body=table.createTBody();states.forEach((state,i)=>{const tr=body.insertRow();const label=tr.insertCell();label.textContent=`${i+1}. `+data.expressions[i].replace(/K\d{5}/g,k=>`${D.genes[k]||k} [${k}]`);tr.insertCell().textContent=c.kind==='sample'?['Not met in observed annotations','Detected','Unresolvable: required row unavailable'][raw[i]]:`${(100*state[1]).toFixed(1)}% detected · ${(100*state[0]).toFixed(1)}% not met · ${(100*state[2]).toFixed(1)}% unresolvable`;});
 const note=document.createElement('p');note.className='subtle';note.textContent='“+” joins required subunits; “,” indicates alternatives; a minus sign marks an optional component. Bounds retain unknown blocks and are not confidence intervals, pathway abundance or same-organism completeness. Plus-linked complexes remain one logical block. '+(c.kind==='sample'?'':'Percentages weight BioSample labels within study–coordinate groups; the overview weights groups equally.');
 details.append(table,note);parent.append(coverage,details);
}

function current(){const v=$('selection').value;if(v==='global')return {values:D.global_values,meta:null,id:'All groups · median summary',kind:'global'};const i=v.indexOf(':'),kind=v.slice(0,i),id=v.slice(i+1);return {values:kind==='group'?D.group_values[id]:D.sample_values[id],meta:kind==='group'?D.groups[id]:D.samples[id],id,kind};}
const coords={};Object.entries(D.groups).forEach(([id,g])=>{if(!coords[g.coordinate_group])coords[g.coordinate_group]=[];coords[g.coordinate_group].push([id,g]);});
function updateMap(){updateDiagnostic();const ko=$('mapgene').value,c=current();$('map-gene-title').textContent=genes[ko]+' · '+ko;document.querySelectorAll('.point').forEach(el=>{const id=el.dataset.coord,v=D.coordinate_values[id][ko],g=coords[id][0][1];el.setAttribute('fill',color(v));el.classList.toggle('selected',Boolean(c.meta&&c.meta.coordinate_group===id));const t=`${g.Country} · ${g.Locality}\n${g.lat.toFixed(5)}, ${g.lon.toFixed(5)}\n${genes[ko]}: ${fmt(v)}${v==null?'':' assignments per million'}\n${coords[id].length} study group${coords[id].length===1?'':'s'}`;el.querySelector('title').textContent=t;el.setAttribute('aria-label',t);const marker=liveMarkers.get(id);if(marker){marker.setStyle({fillColor:color(v),color:c.meta&&c.meta.coordinate_group===id?'#0c1722':'#342440',weight:c.meta&&c.meta.coordinate_group===id?3:1});const label=document.createElement('span');label.style.whiteSpace='pre-line';label.textContent=t;marker.bindTooltip(label);}});}
function updateEnvironment(m){$('environment').replaceChildren();[['Temperature',m?m.bio1:null,'°C'],['Organic carbon',m&&m.soc_mean!=null?m.soc_mean/10:null,'g/kg'],['pH',m&&m.phh2o_mean!=null?m.phh2o_mean/10:null,'']].forEach(([label,value,unit])=>{const box=document.createElement('div');box.className='env-item';const l=document.createElement('span'),v=document.createElement('strong');l.textContent=label;v.textContent=m?(value==null?'Unavailable':fmt(value)):'—';if(value!=null&&unit){const u=document.createElement('small');u.textContent=' '+unit;v.appendChild(u);}box.append(l,v);$('environment').appendChild(box);});}
function update(){const c=current(),m=c.meta;$('selection-title').textContent=c.kind==='group'?`${m.Country} · ${m.Locality} · ${m.study_id}`:c.id;updateEnvironment(m);$('numbers').textContent=m?`${m.lat.toFixed(5)}, ${m.lon.toFixed(5)} · ${c.kind==='group'?m.n_runs+' runs · '+m.n_biosample_labels+' BioSample labels':'Individual metagenome'} · ${c.id}`:'Select a group or metagenome to view its environmental context.';
document.querySelectorAll('[data-ko-fill]').forEach(el=>el.setAttribute('fill',color(c.values[el.dataset.koFill])));document.querySelectorAll('[data-ko-text]').forEach(el=>el.setAttribute('fill',foreground(c.values[el.dataset.koText])));document.querySelectorAll('[data-ko-value]').forEach(el=>{const v=c.values[el.dataset.koValue];el.setAttribute('fill',foreground(v));el.textContent=v==null?'No row':Number(v).toFixed(1);});
const cards=$('cards');cards.replaceChildren();D.reactions.forEach((r,idx)=>{const card=document.createElement('article');card.className='card'+(['ambiguous','unresolved'].includes(r.status)?' ambiguous':'');const n=document.createElement('div');n.className='card-index';n.textContent=r.id+' / '+(r.status==='unresolved'?'UNRESOLVED':'GENETIC POTENTIAL');card.id='process-'+r.id;const h=document.createElement('h3');h.textContent=r.label;card.append(n,h);const rx=document.createElement('div');rx.className='reaction';rx.textContent=r.reaction;card.appendChild(rx);const tiles=document.createElement('div');tiles.className='tiles';r.kos.forEach((ko,i)=>{const v=c.values[ko],b=document.createElement('button');b.type='button';b.className='tile'+(v==null?' missing':'');if(v!=null)b.style.background=color(v);b.style.color=foreground(v);b.setAttribute('aria-label',`Show ${r.genes[i]}, ${ko}, on the map. ${v==null?'KO row unavailable':fmt(v)+' assignments per million'}`);const name=document.createElement('b');name.textContent=genes[ko];b.title=D.full_names[ko];const kid=document.createElement('small');kid.textContent=ko;const val=document.createElement('span');val.textContent=v==null?'No row':fmt(v);b.append(name,kid,val);b.onclick=()=>{$('mapgene').value=ko;updateMap();$('atlas').scrollIntoView?.({behavior:window.matchMedia?.('(prefers-reduced-motion: reduce)').matches?'auto':'smooth'});};tiles.appendChild(b);});card.appendChild(tiles);const note=document.createElement('p');note.className='note';note.textContent=r.note;card.appendChild(note);
if(r.rules.length){const details=document.createElement('details');const sum=document.createElement('summary');sum.textContent=r.rules.length+' operational rule(s) · inspect evidence';details.appendChild(sum);r.rules.forEach(rule=>{const p=document.createElement('div');p.className='rule-result';p.dataset.rule=rule.id;const title=document.createElement('strong');title.textContent=rule.id;const state=document.createElement('div');state.textContent=ruleText(rule,c);const expr=document.createElement('code');expr.textContent=rule.expression;p.append(title,state,expr);appendBlockEvidence(p,rule,c);details.appendChild(p);});card.appendChild(details);}else{const p=document.createElement('p');p.className='rule-result';p.textContent='No operational rule: no zero, prevalence or abundance is assigned to this process.';card.appendChild(p);}const a=document.createElement('a');a.href=r.source;a.textContent='Reference source ↗';a.target='_blank';a.rel='noopener';card.appendChild(a);cards.appendChild(card);});updateMap();}
// Leaflet supplies Web Mercator projection, tile selection, zoom, pan and touch handling.
// Point coordinates are never displaced to separate coincident samples.
let suppressMapClick=false;
const pointElements=[...document.querySelectorAll('#coordinate-records .point')];
function selectPoint(el){const found=coords[el.dataset.coord][0];$('selection').value='group:'+found[0];update();openPointChooser(el);}
function screenPosition(el){const g=coords[el.dataset.coord][0][1],p=geoMap.latLngToContainerPoint([g.lat,g.lon]),b=$('map').getBoundingClientRect();return {x:p.x+b.left,y:p.y+b.top};}
function markChoices(){document.querySelectorAll('.sample-choice').forEach(b=>b.classList.toggle('active',b.dataset.selection===$('selection').value));}
function openPointChooser(el){
 const anchor=screenPosition(el),nearby=pointElements.map(p=>{const q=screenPosition(p);return {el:p,d:Math.hypot(q.x-anchor.x,q.y-anchor.y)};}).filter(o=>o.d<=16).sort((a,b)=>a.d-b.d||a.el.dataset.coord.localeCompare(b.el.dataset.coord));
 const ids=new Set(nearby.map(o=>o.el.dataset.coord)),runs=Object.entries(D.samples).filter(([id,m])=>ids.has(m.coordinate_group));
 $('chooser-title').textContent=`${runs.length} metagenomes at ${ids.size} nearby coordinate${ids.size===1?'':'s'}`;
 const list=$('chooser-list');list.replaceChildren();
 function choice(parent,value,title,subtitle){const b=document.createElement('button');b.type='button';b.className='sample-choice';b.dataset.selection=value;b.textContent=title;const sub=document.createElement('span');sub.textContent=subtitle;b.appendChild(sub);b.addEventListener('click',()=>{$('selection').value=value;update();markChoices();});parent.appendChild(b);}
 nearby.forEach(({el:p})=>{const id=p.dataset.coord,g=coords[id][0][1],section=document.createElement('section');section.className='coordinate-choice';const h=document.createElement('h4');h.textContent=`${g.Country} · ${g.Locality} · ${g.lat.toFixed(5)}, ${g.lon.toFixed(5)}`;section.appendChild(h);const buttons=document.createElement('div');buttons.className='choice-buttons';
  coords[id].forEach(([gid,m])=>choice(buttons,'group:'+gid,'Group · '+m.study_id,m.n_runs+' runs'));
  runs.filter(([rid,m])=>m.coordinate_group===id).sort((a,b)=>a[0].localeCompare(b[0])).forEach(([rid,m])=>choice(buttons,'sample:'+rid,rid,m.study_id));
  section.appendChild(buttons);list.appendChild(section);
 });$('point-chooser').hidden=false;markChoices();
}
$('chooser-close').addEventListener('click',()=>{$('point-chooser').hidden=true;});

geoMap=L.map('map',{zoomControl:false,minZoom:2,maxZoom:17,scrollWheelZoom:true,preferCanvas:false,attributionControl:true}).setView([12,5],2);
// Embedded relief is tiled locally so high zoom never creates a giant image element.
// Its Web Mercator projection matches the online tiles; resampling adds no detail.
const reliefPane=geoMap.createPane('offlineRelief');reliefPane.style.zIndex='140';reliefPane.style.pointerEvents='none';
const reliefImage=new Image();
const reliefReady=new Promise((resolve,reject)=>{reliefImage.onload=resolve;reliefImage.onerror=reject;});
reliefImage.src=D.offline_relief.image;
const EmbeddedRelief=L.GridLayer.extend({
 createTile:function(coords,done){
  const tile=document.createElement('canvas');tile.width=tile.height=256;tile.setAttribute('aria-hidden','true');
  reliefReady.then(()=>{try{const span=reliefImage.naturalWidth/Math.pow(2,coords.z);tile.getContext('2d').drawImage(reliefImage,coords.x*span,coords.y*span,span,span,0,0,256,256);done(null,tile);}catch(error){done(error,tile);}},error=>done(error,tile));
  return tile;
 }
});
const offlineRelief=new EmbeddedRelief({pane:'offlineRelief',tileSize:256,maxZoom:17,opacity:1,updateWhenIdle:true,keepBuffer:1}).addTo(geoMap);
const overviewPane=geoMap.createPane('offlineOverview');overviewPane.style.zIndex='150';overviewPane.style.pointerEvents='none';
const overview=L.geoJSON(D.world,{pane:'offlineOverview',interactive:false,style:{color:'#a8b6b7',weight:.45,fill:false,fillOpacity:0}}).addTo(geoMap);
geoMap.attributionControl.setPrefix('<a href="https://leafletjs.com/" title="Interactive map library">Leaflet</a>');
geoMap.attributionControl.addAttribution('<a href="https://www.naturalearthdata.com/">Natural Earth</a> embedded relief');
const terrain=L.tileLayer('https://tile.opentopomap.org/{z}/{x}/{y}.png',{maxNativeZoom:17,maxZoom:17,opacity:1,updateWhenIdle:true,keepBuffer:1,attribution:'Map data © <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors, SRTM · Map style © <a href="https://opentopomap.org/">OpenTopoMap</a> (<a href="https://creativecommons.org/licenses/by-sa/3.0/">CC BY-SA 3.0</a>)'});
let loadedTiles=0,tileErrors=0;
function terrainStatus(){
 if(!$('relief-toggle').checked)return;
 $('map-status').textContent=loadedTiles?(tileErrors?'Detailed terrain partly loaded · embedded colour relief fills unavailable tiles.':'Detailed terrain · tiles update when you zoom · internet connected'):'Embedded colour relief · online detail unavailable. Reload terrain or open the downloaded HTML in a browser with internet.';
}
terrain.on('tileload',()=>{loadedTiles++;terrainStatus();});
terrain.on('tileerror',()=>{tileErrors++;terrainStatus();});
terrain.on('loading',()=>{loadedTiles=0;tileErrors=0;$('map-status').textContent='Embedded colour relief · loading online detail for this view…';});
terrain.on('load',terrainStatus);
if(navigator.onLine!==false)terrain.addTo(geoMap);else{$('map-status').textContent='Embedded colour relief · offline. High-resolution detail will load when connected.';}
function retryTerrain(){
 $('relief-toggle').checked=true;loadedTiles=0;tileErrors=0;updateRelief();terrain.redraw();
 $('map-status').textContent='Embedded colour relief · retrying online detail…';
}
$('terrain-retry').addEventListener('click',retryTerrain);
window.addEventListener('online',()=>{if($('relief-toggle').checked)retryTerrain();});
pointElements.forEach(el=>{const id=el.dataset.coord,g=coords[id][0][1];const marker=L.circleMarker([g.lat,g.lon],{radius:6,color:'#342440',weight:1,fillOpacity:1,bubblingMouseEvents:false}).addTo(geoMap);liveMarkers.set(id,marker);marker.on('click',()=>selectPoint(el));if(marker._path){marker._path.classList.add('sample-marker');marker._path.setAttribute('data-coordinate',id);marker._path.setAttribute('tabindex','0');marker._path.setAttribute('role','button');marker._path.setAttribute('aria-label','Select '+g.Country+', '+g.Locality);marker._path.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();selectPoint(el);}});}});
L.control.scale({imperial:false,maxWidth:140}).addTo(geoMap);
function renderMapView(){const z=geoMap.getZoom();$('zoom-level').textContent='Zoom '+z;$('zoom-in').disabled=z>=17;$('zoom-out').disabled=z<=2;}
function updateRelief(){if($('relief-toggle').checked){if(!geoMap.hasLayer(terrain))terrain.addTo(geoMap);terrain.setOpacity(Number($('relief-opacity').value)/100);}else{geoMap.removeLayer(terrain);$('map-status').textContent='Embedded colour relief · online detail disabled';}$('relief-opacity').disabled=!$('relief-toggle').checked;}
$('zoom-in').addEventListener('click',()=>geoMap.zoomIn());$('zoom-out').addEventListener('click',()=>geoMap.zoomOut());$('zoom-reset').addEventListener('click',()=>geoMap.setView([12,5],2));
$('zoom-selection').addEventListener('click',()=>{const m=current().meta;if(!m){$('map-status').textContent='Choose a sample or group first.';return;}geoMap.setView([m.lat,m.lon],Math.max(geoMap.getZoom(),12));});
geoMap.on('zoomend moveend',renderMapView);window.addEventListener('resize',()=>geoMap.invalidateSize());
$('relief-toggle').addEventListener('change',updateRelief);$('relief-opacity').addEventListener('input',updateRelief);
$('selection').addEventListener('change',()=>{update();markChoices();});$('mapgene').addEventListener('change',updateMap);$('threshold').addEventListener('change',update);
document.querySelectorAll('.diagram-gene').forEach(el=>{const select=()=>{$('mapgene').value=el.dataset.ko;updateMap();$('atlas').scrollIntoView?.({behavior:'smooth'});};el.addEventListener('click',select);el.addEventListener('keydown',e=>{if(e.key==='Enter')select();});});
$('heatmap-zoom').addEventListener('change',()=>{$('heatmap-image').style.width=$('heatmap-zoom').value+'%';});
document.querySelectorAll('a[href="#methods"]').forEach(a=>a.addEventListener('click',()=>{$('methods').open=true;}));
update();renderMapView();
// Expose a minimal testing interface; this does not change scientific data.
window.atlasMap={map:geoMap,terrain,overview,offlineRelief,reliefImage,markers:liveMarkers,tileSource:'https://tile.opentopomap.org/{z}/{x}/{y}.png'};
