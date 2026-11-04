const fs=require('fs'),path=require('path'),assert=require('assert');
const {JSDOM,VirtualConsole}=require(process.env.JSDOM_MODULE||'jsdom');
const root=process.env.ATLAS_TEST_ROOT||path.resolve(__dirname,'..');
const cycles=JSON.parse(fs.readFileSync(path.join(root,'curation/cycles.json')));
const results=[];
for(const [code,c] of Object.entries(cycles)){
 const errors=[],vc=new VirtualConsole();vc.on('jsdomError',e=>errors.push(e.message));
 const dom=new JSDOM(fs.readFileSync(path.join(root,'figures',c.slug,'viewer.html'),'utf8'),{
  runScripts:'dangerously',pretendToBeVisual:true,virtualConsole:vc,
  beforeParse(w){
   w.SVGSVGElement.prototype.createSVGRect=()=>({});
   w.SVGElement.prototype.getScreenCTM=()=>({a:1,inverse:()=>({a:1})});
   w.HTMLElement.prototype.scrollIntoView=()=>{};
   Object.defineProperty(w.HTMLElement.prototype,'clientWidth',{get(){return this.id==='map'?1000:1000}});
   Object.defineProperty(w.HTMLElement.prototype,'clientHeight',{get(){return this.id==='map'?520:500}});
   w.HTMLElement.prototype.getBoundingClientRect=function(){return {x:0,y:0,left:0,top:0,right:1000,bottom:520,width:1000,height:520}};
  }
 });
 const w=dom.window,d=w.document,api=w.atlasMap;
 if(!api){console.log(code,errors);process.exit(1);}
 const geographicLayers=Object.values(api.map._layers).filter(l=>l.feature&&l.getPane);
 const paneZ=el=>Number(w.getComputedStyle(el).zIndex);
 function correctLayerOrder(){const tz=paneZ(api.terrain.getPane());return geographicLayers.length>0&&geographicLayers.every(l=>paneZ(l.getPane())<tz)&&[...api.markers.values()].every(m=>paneZ(m.getPane())>tz);}
 const layerOrderAtStart=correctLayerOrder();
 const fallbackIgnoresClicks=geographicLayers.every(l=>l.options.interactive===false&&w.getComputedStyle(l.getPane()).pointerEvents==='none');
 const actualPaneOwnership=geographicLayers.every(l=>l._path?.closest('.leaflet-pane')===l.getPane())&&[...api.markers.values()].every(m=>m._path.closest('.leaflet-pane')===m.getPane());
 const terrainFullyOpaque=api.terrain.options.opacity===1&&d.getElementById('relief-opacity').value==='100';
 const layerLevels={overview:[...new Set(geographicLayers.map(l=>paneZ(l.getPane())))],terrain:paneZ(api.terrain.getPane()),samples:[...new Set([...api.markers.values()].map(m=>paneZ(m.getPane())))]};
 const sel=d.getElementById('selection');sel.value=sel.options[sel.options.length-1].value;sel.dispatchEvent(new w.Event('change'));
 const datum=w.eval('D');
 const sampleValues=datum.sample_values[sel.value.slice(7)];
 const cycleValuesUpdate=[...d.querySelectorAll('[data-ko-value]')].every(el=>el.textContent===(sampleValues[el.dataset.koValue]==null?'No row':Number(sampleValues[el.dataset.koValue]).toFixed(1)));
 const checks={cycle:code,layerOrderAtStart,fallbackIgnoresClicks,actualPaneOwnership,terrainFullyOpaque,layerLevels,cycleValuesUpdate,errors,options:sel.options.length,liveMarkers:api.markers.size,cards:d.querySelectorAll('#cards .card').length,sampleUpdates:d.getElementById('selection-title').textContent===sel.value.slice(7),cycleLinks:d.querySelectorAll('.cycle-switch a').length};
 d.getElementById('threshold').value='2';d.getElementById('threshold').dispatchEvent(new w.Event('change'));checks.rulesRendered=d.querySelectorAll('.rule-result').length;
 d.getElementById('zoom-in').click();checks.zoomButton=api.map.getZoom()===3;
 d.getElementById('zoom-selection').click();checks.zoomToSelection=api.map.getZoom()===12;
 [...api.markers.values()][0].fire('click');checks.overlapChooser=!d.getElementById('point-chooser').hidden&&d.querySelectorAll('.sample-choice').length>0;
 const sampleChoice=[...d.querySelectorAll('.sample-choice')].find(b=>b.dataset.selection.startsWith('sample:'));sampleChoice.click();checks.overlapSampleAccessible=sel.value===sampleChoice.dataset.selection;
 const fillsBefore=[...api.markers.values()].map(m=>m.options.fillColor);
 const tile=[...d.querySelectorAll('.tile')].find(t=>t.querySelector('small').textContent!==d.getElementById('mapgene').value);tile.click();checks.geneUpdates=d.getElementById('mapgene').value===tile.querySelector('small').textContent;
 checks.markerColoursUpdate=[...api.markers.values()].some((m,i)=>m.options.fillColor!==fillsBefore[i]);
 // Catch the CSS regression that a DOM-only event test previously missed.
 // Check sizing selectors in every media branch, including narrow-screen rules.
 api.map.setView([12,5],2,{animate:false});
 const overlay=d.querySelector('#map .leaflet-overlay-pane>svg');
 const computed=w.getComputedStyle(overlay),conflictingSizingRules=[];
 function inspectRules(rules){for(const rule of rules){if(rule.selectorText&&rule.style){let matches=false;try{matches=overlay.matches(rule.selectorText)}catch{};if(matches){for(const prop of ['width','height','min-width','min-height','max-width','max-height']){const value=rule.style.getPropertyValue(prop);if(value&&!(['max-width','max-height'].includes(prop)&&value==='none'))conflictingSizingRules.push({selector:rule.selectorText,property:prop,value});}}}if(rule.cssRules)inspectRules(rule.cssRules);}}
 for(const sheet of d.styleSheets)inspectRules(sheet.cssRules);
 checks.noLeafletSizingOverride=conflictingSizingRules.length===0&&computed.width!=='100%'&&computed.height!=='auto';
 checks.conflictingSizingRules=conflictingSizingRules;
 checks.nativeOverlayDimensions=Number(overlay.getAttribute('width'))>0&&Number(overlay.getAttribute('height'))>0&&overlay.getAttribute('viewBox').split(/\s+/).slice(2).every(v=>Number(v)>0);
 checks.compactAttribution=d.querySelectorAll('#map .leaflet-control-attribution svg').length===0&&d.querySelector('#map .leaflet-control-attribution').textContent.includes('Leaflet');
 checks.worldVisibleMarkers=[...api.markers.values()].filter(m=>api.map.getBounds().contains(m.getLatLng())).length;
 checks.markerPathsPresent=checks.worldVisibleMarkers===Object.keys(datum.coordinate_values).length&&d.querySelectorAll('#map path.sample-marker').length===Object.keys(datum.coordinate_values).length&&[...d.querySelectorAll('#map path.sample-marker')].every(p=>p.getAttribute('d')?.length>10&&w.getComputedStyle(p).display!=='none');
 api.map.setView([18.645,-91.83],17,{animate:false});checks.highZoom=api.map.getZoom()===17;checks.layerOrderAtHighZoom=correctLayerOrder();
 checks.nativeZoom=api.terrain.options.maxNativeZoom===17;
 const coords={x:32101,y:58620,z:17};const url=api.terrain.getTileUrl(coords);checks.highZoomTile=url.includes('/17/');
 checks.tilesAtCurrentZoom=Object.keys(api.terrain._tiles).some(k=>k.endsWith(':17'));
 api.terrain.fire('tileerror');checks.failureFallback=d.getElementById('map-status').textContent.includes('offline');
 api.terrain.fire('tileload');checks.terrainStatus=d.getElementById('map-status').textContent.includes('Detailed terrain');
 const toggle=d.getElementById('relief-toggle');toggle.checked=false;toggle.dispatchEvent(new w.Event('change'));checks.offlineOutline=api.map.hasLayer(api.terrain)===false&&Object.values(api.map._layers).some(l=>l.feature);
 toggle.checked=true;toggle.dispatchEvent(new w.Event('change'));
 checks.layerOrderAfterToggle=api.map.hasLayer(api.terrain)&&correctLayerOrder();
 const opacity=d.getElementById('relief-opacity');opacity.value='65';opacity.dispatchEvent(new w.Event('input'));checks.opacityControl=api.terrain.options.opacity===.65;opacity.value='100';opacity.dispatchEvent(new w.Event('input'));
 checks.embeddedCycleGenes=d.querySelectorAll('[data-ko-fill]').length;
 checks.passed=errors.length===0&&checks.options===1+Object.keys(datum.samples).length+Object.keys(datum.groups).length&&checks.liveMarkers===Object.keys(datum.coordinate_values).length&&checks.cycleLinks===8&&checks.cards>0&&checks.embeddedCycleGenes>0&&['sampleUpdates','zoomButton','zoomToSelection','overlapChooser','overlapSampleAccessible','geneUpdates','highZoom','nativeZoom','highZoomTile','tilesAtCurrentZoom','failureFallback','terrainStatus','offlineOutline','markerColoursUpdate','noLeafletSizingOverride','nativeOverlayDimensions','compactAttribution','markerPathsPresent','cycleValuesUpdate','layerOrderAtStart','fallbackIgnoresClicks','actualPaneOwnership','terrainFullyOpaque','layerOrderAtHighZoom','layerOrderAfterToggle','opacityControl'].every(k=>checks[k]);
 results.push(checks);w.close();
}
const out={runtime:'Real bundled Leaflet 1.9.4 in jsdom 26; element sizes and SVG support mocked. Interaction, all-media CSS isolation, SVG dimensions, actual layer panes/z-order, opacity, toggle and tile-URL tests; not a browser pixel-layout test. No tile downloads by this test.',results};
const name=process.env.ATLAS_QA_OUTPUT||path.join(__dirname,'DOM_validation.json');fs.writeFileSync(name+'.tmp',JSON.stringify(out,null,2));fs.renameSync(name+'.tmp',name);
console.log(JSON.stringify(results,null,2));if(results.some(r=>!r.passed))process.exit(1);
