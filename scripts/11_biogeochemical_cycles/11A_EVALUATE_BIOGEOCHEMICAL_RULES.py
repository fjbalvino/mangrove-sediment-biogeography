"""Evaluate frozen KO rules at community level with explicit unknown values.

See docs/BIOGEOCHEMICAL_CYCLES.md for the input contract and commands.
No network calls. No module direction or same-genome completeness is inferred.
"""
from pathlib import Path
import argparse, csv, gzip, hashlib, json, re, sys, platform
import numpy as np
import pandas as pd
import scipy
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[2]
DEFAULT_CONFIG=ROOT/"config/biogeochemical_cycles"

def split_top(s,mode):
    depth=0; start=0; out=[]
    for i,c in enumerate(s):
        if c=='(': depth+=1
        elif c==')': depth-=1
        if depth<0: raise ValueError('Unbalanced parentheses')
        if depth==0 and ((mode=='space' and c.isspace()) or c==mode):
            if s[start:i].strip():out.append(s[start:i].strip())
            elif mode != 'space':raise ValueError('Empty alternative')
            start=i+1
    if depth:raise ValueError('Unbalanced parentheses')
    if s[start:].strip():out.append(s[start:].strip())
    elif mode != 'space':raise ValueError('Empty alternative')
    return out

def unwrap(s):
    while s.startswith('(') and s.endswith(')'):
        depth=0; closed=False
        for i,c in enumerate(s):
            depth+=(c=='(')-(c==')')
            if depth==0:
                closed=(i==len(s)-1);break
        if not closed:break
        s=s[1:-1].strip()
    return s

def parse(s,modules=None,stack=()):
    s=unwrap(s.strip())
    if not s:raise ValueError('Empty expression')
    # KEGG: top-level pathway blocks are each parenthesized, then ANDed.
    parts=split_top(s,'space')
    if len(parts)>1:return ('and',[parse(p,modules,stack) for p in parts])
    parts=split_top(s,',')
    if len(parts)>1:return ('or',[parse(p,modules,stack) for p in parts])
    depth=0; start=0; sign='+'; terms=[]
    for i,c in enumerate(s):
        depth+=(c=='(')-(c==')')
        if depth==0 and c in '+-':
            if i>start:terms.append((sign,s[start:i]))
            elif i!=0:raise ValueError('Consecutive operators')
            sign=c;start=i+1
    if start>0:
        terms.append((sign,s[start:]))
        children=[('optional',[parse(v,modules,stack)]) if op=='-' else parse(v,modules,stack) for op,v in terms]
        return ('and',children)
    if re.fullmatch(r'K\d{5}',s):return ('ko',s)
    if re.fullmatch(r'M\d{5}',s):
        if modules is None or s not in modules:raise ValueError('Unresolved module '+s)
        if s in stack:raise ValueError('Recursive module '+s)
        return parse(modules[s],modules,stack+(s,))
    raise ValueError('Unsupported token: '+s)

def kos(tree,include_optional=False):
    op,x=tree
    if op=='ko':return {x}
    if op=='optional' and not include_optional:return set()
    return set().union(*(kos(c,include_optional) for c in x))

def evaluate(tree,values,n):
    op,x=tree
    if op=='ko':return values.get(x,np.full(n,2,dtype=np.int8))
    if op=='optional':return np.ones(n,dtype=np.int8)
    v=np.stack([evaluate(c,values,n) for c in x])
    if op=='and':return np.where((v==0).any(axis=0),0,np.where((v==2).any(axis=0),2,1))
    if op=='or':return np.where((v==1).any(axis=0),1,np.where((v==2).any(axis=0),2,0))
    raise ValueError(op)

def self_tests():
    tests=[
        ('K00001+K00002-K00003',{'K00001':1,'K00002':1},1),
        ('K00001+K00002-K00003',{'K00001':1,'K00002':0},0),
        ('K00001,K00002',{'K00001':0},2),
        ('K00001,K00002',{'K00001':1},1),
        ('K00001 K00002',{'K00001':0},0),
        ('K00001 K00002',{'K00001':1},2),
        ('(K00001 K00002),K00003',{'K00001':0,'K00002':1,'K00003':1},1),
        ('K00001,K00002 K00003',{'K00001':1,'K00002':0,'K00003':0},0),
        ('(K00001+K00002-K00003)+(K00004,K00005)',{'K00001':1,'K00002':1,'K00005':1},1),
        ('K00001-(K00002+K00003)',{'K00001':1},1),
    ]
    for expr,d,expected in tests:
        actual=evaluate(parse(expr),{k:np.array([v]) for k,v in d.items()},1)[0]
        assert actual==expected,(expr,actual,expected)
    assert kos(parse('K00001+K00002-K00003'))=={'K00001','K00002'}
    assert kos(parse('K00001+K00002-K00003'),True)=={'K00001','K00002','K00003'}
    assert evaluate(parse('M00001',{'M00001':'K00001,K00002'}),{'K00002':np.array([1])},1)[0]==1
    for bad in ['K00001/K00002','K00001+','K00001++K00002','(K00001']:
        try:parse(bad)
        except ValueError:pass
        else:raise AssertionError('Invalid expression accepted: '+bad)
    # Exhaustively compare all ternary inputs with independent scalar logic.
    import itertools
    for a,b,c in itertools.product((0,1,2),repeat=3):
        ab=0 if 0 in (a,b) else 2 if 2 in (a,b) else 1
        want=1 if 1 in (ab,c) else 2 if 2 in (ab,c) else 0
        actual=evaluate(parse('(K00001+K00002),K00003'),dict(zip(['K00001','K00002','K00003'],map(lambda x:np.array([x]),[a,b,c]))),1)[0]
        assert want==actual
    return len(tests)+27+7

def read_counts(path, ids):
    """Reject duplicate headers before pandas can silently rename them."""
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as handle:
        header = next(csv.reader(handle, delimiter='\t'))
    if len(header) < 2 or len(set(header[1:])) != len(header[1:]):
        raise ValueError('Matrix needs unique sample columns')
    matrix = pd.read_csv(path, sep='\t', index_col=0)
    if matrix.empty or not matrix.index.is_unique:
        raise ValueError('Matrix needs nonempty, unique KO rows')
    if not set(ids) <= set(matrix.columns):
        raise ValueError('Sample manifest contains IDs missing from the matrix')
    matrix = matrix.loc[:, ids]
    if not all(re.fullmatch(r'K\d{5}', str(k)) for k in matrix.index):
        raise ValueError('First matrix column must contain KO identifiers K#####')
    try:
        values = matrix.to_numpy(dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError('Matrix must contain numeric protein-record counts') from error
    if not (np.isfinite(values).all() and (values >= 0).all()
            and (values == np.floor(values)).all()):
        raise ValueError('Counts must be finite, nonnegative integers; never round normalized inputs')
    return matrix.astype(float)


def run(a):
    if not __debug__:
        raise RuntimeError('Verification requires Python without -O')
    O=Path(a.out);O.mkdir(parents=True,exist_ok=True)
    tests=self_tests()
    ids=[x.strip() for x in Path(a.samples).read_text().splitlines() if x.strip()]
    if not ids or len(ids) != len(set(ids)):
        raise ValueError('Sample manifest must be nonempty with unique IDs')
    m = read_counts(a.matrix, ids)
    meta = pd.read_csv(a.metadata, keep_default_na=False, dtype=str)
    required = {'run_accession', 'study_id', 'sampling_site_id', 'REALM', 'lat', 'lon'}
    if not required <= set(meta.columns):
        raise ValueError('Missing metadata columns: ' + ', '.join(sorted(required-set(meta.columns))))
    if not meta.run_accession.is_unique:
        raise ValueError('Metadata run_accession must be unique')
    if not set(ids) <= set(meta.run_accession):
        raise ValueError('Sample manifest contains IDs missing from metadata')
    meta = meta.set_index('run_accession').loc[ids].copy()
    for c in meta:
        meta[c] = meta[c].replace(['NA','NaN','nan','N/A','null','None'], '')
    if (meta[['study_id', 'REALM']] == '').any().any():
        raise ValueError('Selected samples need documented study_id and REALM')
    config = Path(a.config_dir)
    rules=pd.read_csv(config/'reglas_operativas.tsv',sep='\t',keep_default_na=False)
    processes=pd.read_csv(config/'procesos_curados.tsv',sep='\t',keep_default_na=False).set_index('id')
    blocked=json.loads((config/'procesos_no_resueltos.json').read_text())
    md=pd.read_csv(config/'modulos_KEGG_definiciones.tsv',sep='\t',keep_default_na=False)
    modules=dict(zip(md.modulo,md.definicion_oficial))
    modules.update(json.loads((config/'module_dependencies.json').read_text()))
    if not rules.rule_id.is_unique or not processes.index.is_unique or not md.modulo.is_unique:
        raise ValueError('Configuration identifiers must be unique')
    if set(rules.rule_id) & set(md.modulo):
        raise ValueError('Rule and reference-module identifiers overlap')
    if not set(rules.process_id) <= set(processes.index) or not set(blocked) <= set(processes.index):
        raise ValueError('Rule/blocked category missing from the process configuration')
    if set(rules.process_id) & set(blocked):
        raise ValueError('Unresolved category must not also have an operational rule')
    module_rules=[]
    for _,row in md.iterrows():
        module_rules.append(dict(rule_id=row.modulo,process_id='referencia',label=row.nombre,expression=row.definicion_oficial,evidence_kind='modulo_referencia',limit='Repertorio comunitario; no capacidad individual, dirección ni flujo.',source=row.url,primary=False))
    combined=pd.concat([rules,pd.DataFrame(module_rules)],ignore_index=True)
    trees={r.rule_id:parse(r.expression,modules) for r in combined.itertuples()}
    # Every definition parses and becomes true when all required KOs are true.
    all_kos=set().union(*(kos(t,True) for t in trees.values()))
    for rid,t in trees.items():assert evaluate(t,{k:np.ones(1,dtype=int) for k in all_kos},1)[0]==1,rid
    # Names are optional display labels; no annotation database is bundled.
    names = {}
    if a.ko_labels:
        labels = pd.read_csv(a.ko_labels, sep='\t', keep_default_na=False)
        if not {'KO', 'name'} <= set(labels.columns) or not labels.KO.is_unique:
            raise ValueError('Optional label file needs unique KO and name columns')
        names = labels.set_index('KO')['name'].to_dict()
    unlabelled = sorted(all_kos-set(names))
    inventory=[]
    for k in sorted(all_kos):
        v=m.loc[k] if k in m.index else None
        inventory.append(dict(KO=k,nombre=names.get(k,k),url='https://www.kegg.jp/entry/'+k,en_matriz=v is not None,n_ge1=None if v is None else int((v>=1).sum()),n_ge2=None if v is None else int((v>=2).sum()),reglas=';'.join(rid for rid,t in trees.items() if k in kos(t,True)),solo_opcional=all(k not in kos(t) for t in trees.values())))
    pd.DataFrame(inventory).to_csv(O/'inventario_KO_panel.tsv',sep='\t',index=False,na_rep='NA')
    m.reindex(sorted(all_kos)).to_csv(O/'conteos_KO_panel.tsv.gz',sep='\t',na_rep='NA',compression='gzip')
    records=[];blocks=[]; rule_summary=[]; arrays={}
    total=m.sum(axis=0)
    for threshold in [1,2]:
        values={k:(m.loc[k].values>=threshold).astype(np.int8) for k in all_kos if k in m.index}
        for r in combined.itertuples():
            tree=trees[r.rule_id]; required=sorted(kos(tree));present=[k for k in required if k in m.index];missing=sorted(set(required)-set(present))
            logical=evaluate(tree,values,len(ids));arrays[(r.rule_id,threshold)]=logical
            detected=np.stack([values[k] for k in present]).sum(axis=0) if present else np.zeros(len(ids),int)
            for j,sample in enumerate(ids):
                observed=[k for k in present if values[k][j]==1]
                under=[k for k in present if values[k][j]==0]
                status=('regla_cumplida' if logical[j]==1 else 'sin_filas' if not present else 'componentes_parciales' if observed else 'sin_deteccion_con_vacios' if missing else 'sin_deteccion')
                records.append(dict(rule_id=r.rule_id,process_id=r.process_id,run_accession=sample,umbral_registros=threshold,resultado_logico={0:'FALSE',1:'TRUE',2:'UNKNOWN'}[int(logical[j])],estado=status,tipo_evidencia=r.evidence_kind,n_KO_requeridos=len(required),n_KO_con_fila=len(present),n_KO_observados=len(observed),KO_observados=';'.join(observed),KO_bajo_umbral=';'.join(under),KO_sin_fila=';'.join(missing)))
            # Blocks separated by top-level spaces; complexes keep their AND semantics.
            for bidx,expr in enumerate(split_top(r.expression,'space'),1):
                bt=parse(expr,modules);bv=evaluate(bt,values,len(ids))
                for j,sample in enumerate(ids):blocks.append(dict(rule_id=r.rule_id,bloque=bidx,expresion=expr,run_accession=sample,umbral_registros=threshold,resultado_logico={0:'FALSE',1:'TRUE',2:'UNKNOWN'}[int(bv[j])]))
            rule_summary.append(dict(rule_id=r.rule_id,process_id=r.process_id,etiqueta=r.label,umbral_registros=threshold,tipo_evidencia=r.evidence_kind,n_cumple=int((logical==1).sum()),n_no_cumple=int((logical==0).sum()),n_indeterminado=int((logical==2).sum()),n_muestras=len(ids),fraccion_cumple=float((logical==1).mean()),n_KO_requeridos=len(required),n_KO_con_fila=len(present),KO_sin_fila=';'.join(missing),limite=r.limit))
    long=pd.DataFrame(records); long.to_csv(O/'evidencia_por_regla_muestra.tsv.gz',sep='\t',index=False,compression='gzip')
    pd.DataFrame(blocks).to_csv(O/'evidencia_por_bloque_muestra.tsv.gz',sep='\t',index=False,compression='gzip')
    summary=pd.DataFrame(rule_summary);summary.to_csv(O/'resumen_reglas.tsv',sep='\t',index=False)
    for rid in combined.rule_id:
        assert not np.any((arrays[(rid,2)]==1)&(arrays[(rid,1)]!=1)),rid
    # Independent checks for the original catalytic sets.
    expected={'N01_nifHDK':['K02588','K02586','K02591'],'M01_MCR':['K00399','K00401','K00402'],'N09_nrfAH':['K03385','K15876'],'C02_ACL':['K15230','K15231']}
    independent_checks = 0
    for rid, ks in expected.items():
        if rid not in trees:
            continue
        for threshold in [1,2]:
            # Independent scalar oracle: missing rows remain UNKNOWN.
            want = []
            for sample in ids:
                parts = [None if k not in m.index else bool(m.at[k,sample] >= threshold) for k in ks]
                want.append(0 if False in parts else 2 if None in parts else 1)
            if not np.array_equal(arrays[(rid,threshold)], want):
                raise AssertionError('Independent catalytic-set check failed: '+rid)
            independent_checks += 1
    # Exact detection statements retained for each of the 65 parent categories.
    parentrows=[]; statuses={1:pd.DataFrame(index=processes.index,columns=ids),2:pd.DataFrame(index=processes.index,columns=ids)}
    for pid,p in processes.iterrows():
        pr=rules[(rules.process_id==pid)&(rules.primary==True)]
        if pid in blocked:
            parentrows.append(dict(process_id=pid,ciclo=p.ciclo,categoria=p.proceso,regla_representativa='',alcance_evaluado='',tipo='no_resuelto',n_reglas=0,n_cumple_ge1=None,n_cumple_ge2=None,limite=blocked[pid]))
            for t in [1,2]:statuses[t].loc[pid,:]='requiere_otra_evidencia'
        else:
            assert len(pr)==1
            r=pr.iloc[0]
            parentrows.append(dict(process_id=pid,ciclo=p.ciclo,categoria=p.proceso,regla_representativa=r.rule_id,alcance_evaluado=r.label,tipo=r.evidence_kind,n_reglas=int((rules.process_id==pid).sum()),n_cumple_ge1=int((arrays[(r.rule_id,1)]==1).sum()),n_cumple_ge2=int((arrays[(r.rule_id,2)]==1).sum()),limite=p.interpretacion_y_limite+' '+r.limit))
            for t in [1,2]:statuses[t].loc[pid,:]=long[(long.rule_id==r.rule_id)&(long.umbral_registros==t)].set_index('run_accession').loc[ids,'estado'].values
    for t,frame in statuses.items():frame.index.name='process_id';frame.to_csv(O/f'matriz_65_reglas_representativas_ge{t}.tsv',sep='\t')
    pd.DataFrame(parentrows).to_csv(O/'alcance_65_procesos.tsv',sep='\t',index=False,na_rep='NA')
    for t in [1,2]:
        table=pd.DataFrame({rid:arrays[(rid,t)] for rid in rules.rule_id},index=ids).T.replace(2,np.nan)
        table.index.name='rule_id';table.to_csv(O/f'matriz_97_reglas_ge{t}.tsv',sep='\t',na_rep='NA')
    # Descriptive design checks. Existing manual sample acceptance is immutable.
    meta['total_asignaciones_KO']=total
    meta['KO_detectados_ge1']=(m>=1).sum(axis=0)
    meta['n_reglas_cumplidas_ge1']=sum(arrays[(rid,1)]==1 for rid in rules.rule_id)
    meta['n_reglas_cumplidas_ge2']=sum(arrays[(rid,2)]==1 for rid in rules.rule_id)
    keep=['sample_accession','study_id','study_id_basis','sampling_site_id','sampling_site_basis','REALM','Country','Locality','lat','lon','total_asignaciones_KO','KO_detectados_ge1','n_reglas_cumplidas_ge1','n_reglas_cumplidas_ge2']
    meta[[k for k in keep if k in meta]].to_csv(O/'muestras_diseno_y_deteccion.tsv',sep='\t',na_rep='NA')
    cross=pd.crosstab(meta.study_id,meta.REALM);cross.to_csv(O/'estudio_por_realm.tsv',sep='\t')
    sites=meta[meta.sampling_site_id!='']
    design_by_realm=[]
    for realm,g in meta.groupby('REALM'):
        design_by_realm.append(dict(REALM=realm,metagenomas=len(g),estudios=g.study_id.nunique(),sitios_identificados=g.loc[g.sampling_site_id!='','sampling_site_id'].nunique(),sin_sitio=int((g.sampling_site_id=='').sum()),sin_coordenadas=int(((g.lat=='')|(g.lon=='')).sum()),mediana_asignaciones=float(g.total_asignaciones_KO.median())))
    pd.DataFrame(design_by_realm).to_csv(O/'diseno_por_realm.tsv',sep='\t',index=False)
    # Rank increment answers whether realm adds independent columns after study.
    study_dummies=pd.get_dummies(meta.study_id,dtype=float).values
    realm_dummies=pd.get_dummies(meta.REALM,dtype=float).values
    rank_study=int(np.linalg.matrix_rank(study_dummies));rank_both=int(np.linalg.matrix_rank(np.c_[study_dummies,realm_dummies]))
    sample_dup=meta[(meta.sample_accession!='')&meta.sample_accession.duplicated(False)] if 'sample_accession' in meta else meta.iloc[:0]
    sample_dup.to_csv(O/'accessions_biologicas_repetidas.tsv',sep='\t')
    # This is a detection diagnostic, not a diversity index or an inference test.
    rho=(float(spearmanr(meta.total_asignaciones_KO,meta.n_reglas_cumplidas_ge1).statistic)
         if meta.total_asignaciones_KO.nunique()>1 and meta.n_reglas_cumplidas_ge1.nunique()>1 else None)
    diagnostics=dict(n_muestras=len(ids),n_estudios=meta.study_id.nunique(),n_sitios_identificados=sites.sampling_site_id.nunique(),sin_sitio=int((meta.sampling_site_id=='').sum()),sin_coordenadas=int(((meta.lat=='')|(meta.lon=='')).sum()),estudios_con_mas_de_un_realm=int(((cross>0).sum(axis=1)>1).sum()),rango_estudio=rank_study,rango_estudio_mas_realm=rank_both,grados_libertad_realm_adicionales=rank_both-rank_study,sample_accessions_repetidas=sample_dup.sample_accession.nunique() if len(sample_dup) else 0,total_asignaciones_min=int(total.min()),total_asignaciones_max=int(total.max()),rho_rendimiento_anotacion_numero_reglas=rho,interpretacion_rho='Diagnóstico descriptivo entre totales de anotación y reglas cumplidas; reglas solapadas, sin P ni lectura causal.',n_reglas=len(rules),n_modulos_referencia=len(md),n_procesos=len(processes),n_procesos_con_regla=len(processes)-len(blocked),n_procesos_no_resueltos=len(blocked),KO_panel=len(all_kos),KO_panel_sin_fila=len(all_kos-set(m.index)),KO_sin_nombre_local=unlabelled,pruebas_logicas=tests,comparaciones_independientes=independent_checks,monotonicidad_umbral=True,unidad='Asignaciones KO de registros de proteínas predichas desde ensamblajes; no lecturas/cobertura/actividad.',versiones=dict(python=platform.python_version(),pandas=pd.__version__,numpy=np.__version__,scipy=scipy.__version__),sha256={key:hashlib.sha256(Path(path).read_bytes()).hexdigest() for key,path in dict(matrix=a.matrix, metadata=a.metadata, samples=a.samples, evaluator=__file__, **{f'config/{p.name}':p for p in sorted(config.iterdir()) if p.is_file()}, **({'ko_labels':a.ko_labels} if a.ko_labels else {})).items()})
    (O/'diagnostico_y_validacion.json').write_text(json.dumps(diagnostics,ensure_ascii=False,indent=2)+'\n')
    print(f"BIOGEOCHEMICAL_RULES_PASS samples={len(ids)} rules={len(rules)} reference_modules={len(md)}")

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--matrix', required=True, help='External KO-by-sample TSV, optionally gzipped')
    ap.add_argument('--metadata', required=True, help='External metadata CSV')
    ap.add_argument('--samples', required=True, help='Ordered sample IDs, one per line')
    ap.add_argument('--unit', required=True, choices=['predicted_protein_records'],
                    help='Declare input provenance; integer values alone do not establish this unit')
    ap.add_argument('--config-dir', default=str(DEFAULT_CONFIG))
    ap.add_argument('--ko-labels', help='Optional external display labels TSV: KO, name')
    ap.add_argument('--out', required=True, help='Local output directory, excluded from Git')
    try:
        run(ap.parse_args())
    except Exception as error:
        print(f'BIOGEOCHEMICAL_RULES_FAIL: {error}', file=sys.stderr)
        raise
