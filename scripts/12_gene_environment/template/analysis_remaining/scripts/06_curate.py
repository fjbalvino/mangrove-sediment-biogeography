"""English interpretation layer; never changes the frozen KO rule expressions."""
from pathlib import Path
import json,re,pandas as pd
r=Path(__file__).resolve().parents[1]
# Process | readable name | reference transformation | inferential limit
text='''C01|Calvin–Benson–Bassham cycle|CO2 → organic carbon|Community components and RuBisCO/PRK do not resolve RuBisCO form, genomic linkage or net autotrophy.
C02|Reductive TCA cycle|CO2 → acetyl-CoA and biomass precursors|ATP citrate lyase is a marker, not a substitute for the pathway; shared central enzymes do not establish direction.
C03|Reverse oxidative TCA|CO2 → biomass precursors|Unresolved: the ordinary TCA enzyme inventory cannot identify reverse operation. Physiological and genomic evidence is needed.
C04|Wood–Ljungdahl pathway|CO2 ⇄ acetyl-CoA|The bacterial module shares enzymes with acetyl-CoA oxidation and CO metabolism. Direction remains unresolved.
C05|3-Hydroxypropionate bicycle|CO2 → biomass precursors|Community co-detection is not a complete pathway within one organism.
C06|3-HP/4-HB cycle|CO2 → biomass precursors|Shared carboxylases and reductases require the whole architecture and genomic context.
C07|Dicarboxylate/4-HB cycle|CO2 → biomass precursors|A component screen does not prove an operating carbon-fixation cycle.
C08|Reductive glycine pathway|CO2 + C1 units → glycine and biomass|Unresolved: glycine-cleavage and C1 enzymes also operate in other directions. Context is required.
C09|Anaplerotic carboxylation|Inorganic carbon → metabolic intermediates|Replenishing intermediates is not equivalent to autotrophic growth or net carbon sequestration.
C10|Cellulose processing|Cellulose → soluble sugars|Unresolved substrate-specific scope: sequence-based CAZyme subfamilies, domains and gene-cluster context are required.
C11|Hemicellulose processing|Xylan, mannan and xyloglucan → sugars|Unresolved: distinguish individual substrates with CAZyme subfamily/domain and cluster evidence.
C12|Other complex carbohydrates|Pectin, starch, chitin and algal glycans → monomers|Unresolved: these are different substrates, not one pathway. Broad KO assignments cannot replace substrate-level CAZyme annotation.
C13|Aromatic monomer processing|Protocatechuate or catechol → ring-cleavage products|PcaGH and CatA support selected aromatic-ring reactions, not lignin depolymerization.
C14|Aerobic hydrocarbon candidates|Alkanes → alcohols|AlkB is a candidate family requiring sequence-specific confirmation; it does not cover all hydrocarbons.
C15|Anaerobic aromatic processing|Benzoyl-CoA → reduced ring intermediates|BcrABCD indicates a dearomatization architecture, not hydrocarbon activation or complete degradation.
C16|Acetate and formate branches|Pyruvate → acetyl-CoA + formate; acetyl-CoA ⇄ acetate|Pfl requires activation context. Pta–AckA is reversible; neither identifies net acetate production.
C17|Lactate and butyrate branches|Pyruvate ⇄ lactate; butyryl-CoA ⇄ butyrate|Selected terminal reactions only. Ethanol, propionate and all alternative fermentation pathways are not resolved by these rules.
C18|Carbon monoxide conversion|CO ⇄ CO2|CoxL requires sequence-level specificity checks; CooS/AcsA overlaps Wood–Ljungdahl and cannot establish direction.
M01|Methyl/alkyl-CoM reductase family|Methyl-CoM ⇄ methane; related alkyl-CoM reactions|McrABG alone cannot distinguish methane production, anaerobic methane oxidation or related alkane metabolism.
M02|Hydrogenotrophic methanogenesis components|CO2 + reducing equivalents → methane|A reference pathway direction, not observed flux; community components may be distributed among organisms.
M03|Acetoclastic methanogenesis components|Acetate → methane + CO2|Shared acetyl-CoA machinery and community co-detection do not prove acetoclastic methane production.
M04|Methylotrophic methanogenesis components|Methanol or methylamines → methane|Alternative substrate architectures remain separate; MCR direction and organismal linkage remain unresolved.
M05|Aerobic methane oxidation candidates|Methane → methanol|Historical Amo/Pmo KOs cannot resolve methane versus ammonia substrate. Soluble MMO is a separate architecture requiring homolog/context checks; MmoD is accessory.
M06|Methanol oxidation families|Methanol → formaldehyde or downstream C1 products|MxaFI and XoxF support methylotrophic potential, not necessarily methane use. XoxF clade and substrate need sequence resolution.
M07|C1 assimilation pathways|C1 intermediates → biomass precursors|Serine, RuMP and XuMP are alternatives. Community pathway components do not imply that one organism carries the complete route.
S01|Dissimilatory sulfate pathway components|Sulfate ⇄ sulfite ⇄ sulfide|Sat/Apr/Dsr components are shared with oxidative pathways. The KO matrix does not establish sulfate-reduction direction or rate.
S02|Reverse-Dsr pathway candidates|Reduced sulfur → sulfite → sulfate|AprAB/DsrAB overlap S01. These are not independent markers and do not diagnose reverse-Dsr without sequence/context.
S03|Assimilatory sulfate reduction|Sulfate → sulfide for biosynthesis|Assimilation is distinct from respiratory sulfate reduction; community components do not identify flux.
S04|Sox thiosulfate oxidation|Thiosulfate → oxidized sulfur products|SoxAXBYZ and the SoxCD extension are separate rules. A core without SoxCD is not evidence of complete oxidation to sulfate.
S05|Sulfide oxidation families|Sulfide → sulfur intermediates|Sqr and FccAB are alternatives. The matrix cannot distinguish energy conservation from detoxification.
S06|Sulfite oxidation families|Sulfite → sulfate|SorAB and SoeABC are alternatives, not jointly required. SorAB is not sulfur oxygenase reductase or a sorbose enzyme.
S07|Anaerobic sulfur-intermediate reduction|Sulfite, thiosulfate or polysulfide → reduced products|Asr, Phs and Psr are separate architectures. PhsA/PsrA share a KO; community co-detection does not resolve substrate specificity.
S08|Tetrathionate and thiosulfate conversions|Tetrathionate ⇄ thiosulfate and other intermediates|Ttr reduction, TsdA oxidation and TetH hydrolysis are distinct. None alone demonstrates sulfur disproportionation.
S09|Persulfide processing|Persulfide → oxidized sulfur intermediates|ETHE1 represents one selected component, not a universal elemental-sulfur or disproportionation pathway.
S10|DMSP metabolism|DMSP → DMS or demethylation products; precursor → DMSP|DmdA, DddL and DsyB represent different steps. DddL does not cover other lyases; no atmospheric flux can be inferred.
S11|DMS, DMSO and methanethiol|DMSO → DMS; DMS → oxidized products; methanethiol → DMS|DmsABC, DdhABC and MddA are separate reactions; they do not establish a net DMS balance.
S12|Selected organic-sulfur utilization|Taurine or alkanesulfonates → released sulfur|TauD and SsuDE cover selected routes only, not all sulfonates, sulfate esters or sulfoquinovose pathways.
S13|Sulfur assimilation into amino acids|Sulfide + carbon skeleton → cysteine|CysE/CysK represent assimilation and biosynthesis, not sulfate respiration or complete methionine metabolism.
P01|Phosphate uptake|External phosphate → cellular phosphate|PstSCAB is transport potential. Presence or higher representation does not establish P limitation or active uptake.
P02|Organic phosphate hydrolysis|Phosphate esters → phosphate|PhoA/PhoD, PhoN/AphA, AppA and GlpQ have different substrates. The panel omits unvalidated PhoX/PafA mappings and does not identify extracellular localization.
P03|C–P lyase architecture|Phosphonate → released phosphorus and carbon products|PhnJ alone does not replace PhnGHIJKLM. Community co-detection does not prove a complete route or methane formation from methylphosphonate.
P04|Alternative phosphonate utilization|Aminoethylphosphonate or phosphonoacetate → phosphorus products|PhnWX, PhnYZ and PhnA are distinct substrate routes; homonymous PhnY enzymes cannot be substituted.
P05|Phosphite oxidation components|Phosphite → phosphate|PtxD and ApdA are different architectures. These rules do not establish energy-yielding phosphite oxidation or hypophosphite metabolism.
P06|Polyphosphate turnover|Cellular phosphate ⇄ polyphosphate|PPK direction depends on enzyme and conditions; PPX/GppA requires context. Do not substitute soluble pyrophosphatase Ppa or infer net storage.
P07|Central nucleotide metabolism|Phosphorylated intermediates ⇄ nucleotides|Context only: general nucleotide enzymes are not a specific index of environmental P mobilization.
Fe01|Ferrous iron oxidation|Fe(II) → Fe(III)|Unresolved: requires FeGenie-style sequence models, diagnostic motifs and genomic neighborhoods; generic cytochromes are insufficient.
Fe02|Ferric iron reduction and EET|Fe(III) → Fe(II)|Unresolved: extracellular electron-transfer systems require sequence and context. Methanogenesis Mtr proteins are not evidence for metal-reduction Mtr systems.
Fe03|Ferrous iron acquisition|External Fe(II) → cellular iron|Only FeoAB is operational here. Transport is not respiratory Fe oxidation/reduction, and does not cover siderophores, heme uptake or storage.
Fe04|Magnetite biomineralization|Cellular iron → magnetite|Unresolved: requires mam/mms sequence and genomic-context evidence.
H01|F420-linked hydrogenase components|H2 ⇄ protons + reducing equivalents|FrhABG is one selected family. Hydrogenase functional groups, physiological direction and coupling require sequence-based classification.
H02|Rnf electron-transfer components|Reduced ferredoxin ⇄ NAD(H)-linked electron transfer|Rnf component co-detection does not prove electron bifurcation, syntrophy or direct interspecies electron transfer.
O01|Terminal oxygen reductases|O2 → water|cbb3 and bd are separate architectures. Gene representation does not establish in situ oxygen concentration, exact affinity or respiratory rate.
O02|Phototrophic reaction centres|Light → charge separation|Photosystems I/II and PufLM are separate architectures. These markers do not establish photosynthetic activity or autotrophy.'''
ps={}
for line in text.splitlines():
    pid,label,reaction,note=line.split('|');ps[pid]={'id':pid,'label':label,'reaction':reaction,'note':note}
cycles={'C':{'name':'Carbon','slug':'carbon','symbol':'C','intro':'Carbon fixation, organic-carbon transformations and selected fermentation branches.'},'CH4_C1':{'name':'Methane & C1','slug':'methane_C1','symbol':'CH₄','intro':'Methane-associated enzyme families, methanol conversion and C1 assimilation.'},'S':{'name':'Sulfur','slug':'sulfur','symbol':'S','intro':'Inorganic and organic sulfur transformations, with shared reductive and oxidative markers exposed.'},'P':{'name':'Phosphorus','slug':'phosphorus','symbol':'P','intro':'Phosphate acquisition, organic-P utilization, phosphite oxidation and polyphosphate turnover.'},'Fe':{'name':'Iron','slug':'iron','symbol':'Fe','intro':'Iron acquisition evidence and the unresolved redox transformations that need sequence-level annotation.'},'H2_electrons':{'name':'Hydrogen & electron transfer','slug':'hydrogen','symbol':'H₂','intro':'Selected hydrogenase and Rnf components; contextual links between biogeochemical transformations.'},'O2_light':{'name':'Respiration & phototrophy','slug':'respiration_light','symbol':'O₂','intro':'Oxygen-reduction and phototrophic reaction-centre components; context for energy metabolism.'}}
src=pd.read_csv(r/'curation/process_scope_original.tsv',sep='\t');rules=pd.read_csv(r/'curation/operational_rules.tsv',sep='\t');inv=pd.read_csv(r/'results/KO_inventory.tsv',sep='\t').set_index('KO')
ref={'C':'10.1186/s40168-021-01213-8','CH4_C1':'10.1111/1755-0998.13589','S':'10.1111/1755-0998.13306','P':'10.1186/s40168-022-01292-1','Fe':'10.3389/fmicb.2020.00037','H2_electrons':'10.1038/srep34212','O2_light':'10.1186/s40168-021-01213-8'}
assert set(ps)==set(src.id)
items=[]
for row in src.itertuples():
    p=ps[row.id];p['cycle']=row.cycle;p['source']='https://doi.org/'+ref[row.cycle]
    rr=rules[rules.process_id==row.id];p['rules']=[]
    for q in rr.itertuples():
        kk=sorted(set(re.findall(r'K\d{5}',q.expression)))
        p['rules'].append({'id':q.rule_id,'expression':q.expression,'kos':kk,'source':q.source,'kind':q.evidence_kind,'primary':bool(q.primary)})
    p['kos']=list(dict.fromkeys(k for rule in p['rules'] for k in rule['kos']));p['genes']=[inv.loc[k,'gene'] for k in p['kos']]
    p['status']='unresolved' if len(rr)==0 else ('ambiguous' if row.id in ['C02','C04','C06','C09','C14','C16','C17','C18','M01','M02','M03','M04','M05','M06','S01','S02','S07','S09','P05','P06','H01','H02'] else 'reference_transformation')
    items.append(p)
(r/'curation/processes_english.json').write_text(json.dumps(items,ensure_ascii=False,indent=2)+'\n')
(r/'curation/cycles.json').write_text(json.dumps(cycles,ensure_ascii=False,indent=2)+'\n')
pd.DataFrame([{'cycle':p['cycle'],'process':p['id'],'label':p['label'],'reaction':p['reaction'],'interpretation':p['note'],'status':p['status'],'rule':rule['id'],'expression':rule['expression'],'KO':k,'name':inv.loc[k,'name'],'source':p['source']} for p in items for rule in (p['rules'] or [{'id':'unresolved','expression':'','kos':['']}]) for k in rule['kos'] if k]).to_csv(r/'results/reaction_gene_map.tsv',sep='\t',index=False)
print('Curated',len(items),'processes; frozen expressions unchanged')
