"""Reference cycle diagrams. Arrows are biochemical references, never measured flux."""
import html,textwrap
G={
'C':(
 [('co2',100,290,'CO₂'),('bio',550,80,'Biomass precursors'),('poly',1030,120,'Complex organic matter'),('mono',1030,480,'Monomers'),('acet',550,500,'Central carbon intermediates'),('co',100,65,'CO')],
 [('co2','bio','Carbon fixation','K00855 K01601',False,270,150),('co2','acet','Wood–Ljungdahl','K00198 K14138',True,270,365),('poly','mono','CAZymes: unresolved','',True,950,285),('mono','acet','Selected aromatic reactions','K00448 K03381',False,680,430),('acet','bio','Anaplerosis / central metabolism','K01595',True,540,255),('co','co2','CO conversion','K03520',True,20,170)]),
'CH4_C1':(
 [('co2',110,90,'CO₂'),('acet',110,480,'Acetate'),('methyl',590,80,'Methyl-CoM'),('ch4',590,480,'CH₄'),('meoh',1050,480,'Methanol'),('c1',1050,90,'C1 intermediates / biomass')],
 [('co2','methyl','CO₂ branch components','K00200',True,270,45),('acet','methyl','Acetyl-CoA branch components','K00193',True,220,300),('methyl','ch4','Shared terminal Mcr / alkyl-CoM family','K00399 K00401 K00402',True,535,225),('ch4','meoh','Amo/Pmo unresolved; sMMO separate','K10944 K16157',True,735,435),('meoh','c1','Methanol conversion','K14028 K23995',False,950,255)]),
'S':(
 [('so4',100,80,'Sulfate'),('so3',580,80,'Sulfite'),('hs',1030,80,'Sulfide'),('s0',1030,300,'Sulfur / persulfides'),('tet',1030,490,'Tetrathionate'),('thio',580,490,'Thiosulfate'),('org',100,490,'Organic sulfur / DMS(P)')],
 [('so4','so3','AprAB: shared direction','K00394 K00395',True,250,45),('so3','hs','DsrAB: shared direction','K11180 K11181',True,720,45),('hs','s0','Sulfide oxidation','K17218',False,950,165),('thio','so4','Sox components','K17222',False,185,155),('org','so3','Organic-S utilization','K03119',False,135,390),('tet','thio','Ttr reduction / TsdA oxidation','K08357 K19713',True,720,445)]),
'P':(
 [('org',100,80,'Organic phosphate'),('phn',590,80,'Phosphonates'),('ptx',1050,80,'Phosphite'),('pi',590,310,'Phosphate pool (Pi)'),('cell',250,520,'Cellular phosphate'),('poly',950,520,'Polyphosphate')],
 [('org','pi','Phosphoester hydrolysis','K01077 K01113',False,165,170),('phn','pi','C–P cleavage / alternatives','K06163',False,545,165),('ptx','pi','Phosphite oxidation','K18916',False,850,170),('pi','cell','Pst phosphate uptake','K02040',False,285,385),('cell','poly','Polyphosphate turnover','K00937 K01524',True,535,475)]),
'Fe':(
 [('fe2',140,150,'External Fe(II)'),('fe3',990,150,'External Fe(III)'),('cell',420,480,'Cellular iron'),('mag',990,480,'Magnetite')],
 [('fe2','fe3','Iron oxidation / reduction: unresolved','',True,440,100),('fe2','cell','FeoAB acquisition','K04758 K04759',False,170,285),('cell','mag','Biomineralization: unresolved','',True,635,435)]),
'H2_electrons':(
 [('h2',140,120,'H₂'),('f420',990,120,'F₄₂₀-linked reducing equivalents'),('fd',140,480,'Ferredoxin'),('nad',990,480,'NAD(H) / ion coupling')],
 [('h2','f420','FrhABG: direction unresolved','K00440 K00441 K00443',True,420,70),('fd','nad','Rnf: electron-transfer components','K03617 K03616 K03615',True,420,430)]),
'O2_light':(
 [('o2',140,120,'O₂'),('water',990,120,'H₂O'),('light',140,480,'Light'),('charge',990,480,'Charge separation')],
 [('o2','water','Terminal oxidases: separate architectures','K00404 K00425',False,420,70),('light','charge','Reaction centres: separate architectures','K02689 K02703 K08928',False,420,430)])}

def make(cycle,inventory,values,color):
    if cycle=='N':
        import xml.etree.ElementTree as ET
        from pathlib import Path
        import re
        ET.register_namespace('', 'http://www.w3.org/2000/svg')
        svg=ET.fromstring((Path(__file__).resolve().parents[1]/'assets/nitrogen_reference_template.svg').read_text())
        for group in svg.iter():
            ident=group.attrib.get('id','')
            match=re.fullmatch(r'(box|label|value)_(K\d{5})',ident)
            if not match:continue
            kind,ko=match.groups();value=values.get(ko);bg=color(value)
            fg='#202636' if sum(int(bg[j:j+2],16)*w for j,w in [(1,.299),(3,.587),(5,.114)])>145 else 'white'
            if kind=='box':group.set('class','diagram-gene');group.set('data-ko',ko);group.set('tabindex','0');group.set('role','button')
            for element in group:
                style=re.sub(r';?\s*fill:\s*[^;]+','',element.attrib.get('style',''));element.set('style',style)
                element.set('data-ko-'+{'box':'fill','label':'text','value':'value'}[kind],ko)
                element.set('fill',bg if kind=='box' else fg)
                if kind=='value':element.text='No row' if value is None or value!=value else f'{value:.1f}'
        svg.set('width','100%');svg.set('height','auto')
        return ET.tostring(svg,encoding='unicode')
    nodes,edges=G[cycle];nd={k:(x,y,label) for k,x,y,label in nodes};parts=[]
    parts.append('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 630" role="img" aria-label="Reference transformations with individual gene representation"><rect width="1200" height="630" fill="white"/><defs><marker id="head" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto-start-reverse"><path d="M0 0 L8 4 L0 8 Z" fill="#6a6480"/></marker></defs>')
    arrows=[];annotations=[]
    for src,dst,label,ks,amb,tx,ty in edges:
        x1,y1,_=nd[src];x2,y2,_=nd[dst];dx=x2-x1;dy=y2-y1;length=(dx*dx+dy*dy)**.5;pad=min(97/(abs(dx)/length) if dx else float('inf'),29/(abs(dy)/length) if dy else float('inf'))+12
        xa,ya=x1+dx/length*pad,y1+dy/length*pad;xb,yb=x2-dx/length*pad,y2-dy/length*pad
        arrows.append(f'<path d="M{xa} {ya} L{xb} {yb}" stroke="{"#a7864a" if amb else "#6a6480"}" stroke-width="2" fill="none" {"stroke-dasharray=\"7,6\"" if amb else ""} {'' if amb else 'marker-end="url(#head)"'}/>')
        lines=textwrap.wrap(label,40);width=max(130,94*len(ks.split()),max(map(len,lines))*6.1+12);tx=min(tx,1180-width);boxheight=(25+len(lines)*15)+(54 if ks else 0)
        annotations.append(f'<rect x="{tx-8}" y="{ty-19}" width="{width+8}" height="{boxheight}" rx="7" fill="white" opacity=".98"/>')
        for j,line in enumerate(lines):annotations.append(f'<text x="{tx}" y="{ty+j*15}" font-family="Arial,sans-serif" font-size="12" fill="#3e3151">{html.escape(line)}</text>')
        ty+=(len(lines)-1)*15
        for i,k in enumerate(ks.split()):
            if k not in inventory.index:continue
            gene={'K10944':'Amo/Pmo A','K10945':'Amo/Pmo B','K10946':'Amo/Pmo C','K00198':'CooS/AcsA','K18916':'PtxD','K08352':'PhsA/PsrA'}.get(k,inventory.loc[k,'gene'].split(',')[0][:15]);v=values.get(k);xx=tx+i*94;yy=ty+12;bg=color(v);fg='#202636' if sum(int(bg[j:j+2],16)*w for j,w in [(1,.299),(3,.587),(5,.114)])>145 else 'white'
            annotations.append(f'<g class="diagram-gene" data-ko="{k}" tabindex="0" role="button" aria-label="Show {k} on map"><rect data-ko-fill="{k}" x="{xx}" y="{yy}" width="89" height="48" rx="5" fill="{bg}" stroke="#9b91a5"/><text data-ko-text="{k}" x="{xx+44.5}" y="{yy+18}" text-anchor="middle" fill="{fg}" font-family="Arial" font-size="11">{html.escape(gene)}</text><text data-ko-value="{k}" x="{xx+44.5}" y="{yy+36}" text-anchor="middle" fill="{fg}" font-family="Arial" font-size="10">{"No row" if v is None else f"{v:.1f}"}</text></g>')
    parts+=arrows+annotations
    for k,x,y,label in nodes:
        parts.append(f'<rect x="{x-97}" y="{y-29}" width="194" height="58" rx="22" fill="#f0eaf7" stroke="#c6b8d5"/>')
        lines=textwrap.wrap(label,25)
        for j,line in enumerate(lines):parts.append(f'<text x="{x}" y="{y+4-(len(lines)-1)*8+j*16}" text-anchor="middle" font-size="13" font-family="Arial" fill="#30203e">{html.escape(line)}</text>')
    parts.append('<text x="25" y="605" font-size="12" font-family="Arial" fill="#695b75">Reference schematic · individual gene tiles, no summed pathway abundance · dashed: direction, assignment or coverage unresolved</text></svg>')
    return ''.join(parts)
