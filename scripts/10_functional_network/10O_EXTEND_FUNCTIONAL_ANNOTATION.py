#!/usr/bin/env python3
"""Frozen, many-to-many KEGG context for the entire observed KO repertoire.

This is identifier mapping, not sequence reannotation or network inference.
Downloads are explicit; annotation itself is offline. Generated references and
tables belong outside Git. KEGG API access is restricted to academic users.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import gzip
import json
import os
from pathlib import Path
import re
import sys
import threading
import tempfile
import time
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent / 'lib'))
from biogeochemical_context import read_matrix, logic, ROOT

BASE = 'https://rest.kegg.jp/'
ENDPOINTS = {
    'info_ko.txt': 'info/ko', 'info_module.txt': 'info/module',
    'info_pathway.txt': 'info/pathway', 'info_brite.txt': 'info/brite',
    'ko_names.tsv': 'list/ko', 'module_names.tsv': 'list/module',
    'pathway_names.tsv': 'list/pathway', 'reaction_names.tsv': 'list/reaction',
    'ko_module.tsv': 'link/module/ko', 'ko_pathway.tsv': 'link/pathway/ko',
    'ko_reaction.tsv': 'link/reaction/ko',
    'ko00001.json': 'get/br:ko00001/json', 'ko00002.json': 'get/br:ko00002/json',
}
BIOLOGICAL_ROOTS = {'09100', '09120', '09130', '09140'}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic_bytes(path, content):
    path = Path(path)
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix='.writing-', delete=False) as handle:
        tmp = Path(handle.name)
        handle.write(content)
        handle.flush()
        os.fsync(handle.fileno())
    tmp.replace(path)


def write_json(path, obj):
    atomic_bytes(path, (json.dumps(obj, indent=2) + '\n').encode())


def write_tsv(frame, path):
    content = frame.to_csv(sep='\t', index=False).encode()
    if str(path).endswith('.gz'):
        content = gzip.compress(content, mtime=0)
    atomic_bytes(path, content)


def names(path, prefix=''):
    result = {}
    for line in Path(path).read_text().splitlines():
        key, value = line.split('\t', 1)
        key = key.split(':')[-1]
        if prefix and key.startswith('map'):
            key = prefix + key[3:]
        if key in result:
            raise ValueError('Duplicate reference identifier: ' + key)
        result[key] = value
    return result


def protein_hierarchies(tree):
    ids = set()
    def walk(node, active=False):
        active = active or node['name'].startswith('09180 ')
        if active:
            ids.update(re.findall(r'\[BR:(ko\d{5})\]', node['name']))
        for child in node.get('children', []):
            walk(child, active)
    walk(tree)
    return sorted(ids)


def fetch_reference(out):
    """Two connections, globally <=2 starts/s; retain raw bytes and hashes."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'manifest.json').exists():
        verify_reference(out)
        print('Existing frozen reference verified; no live update performed.')
        return
    started = datetime.now(timezone.utc).isoformat()
    lock, last = threading.Lock(), [0.0]

    def fetch(item):
        filename, endpoint = item
        target = out / filename
        if target.exists():
            return filename, endpoint
        for attempt in range(3):
            with lock:
                time.sleep(max(0, .55 - (time.monotonic() - last[0])))
                last[0] = time.monotonic()
            try:
                req = Request(BASE + endpoint, headers={'User-Agent': 'AcademicFunctionalContext/1.0'})
                with urlopen(req, timeout=90) as response:
                    data = response.read()
                if not data or data.lstrip().startswith(b'<'):
                    raise ValueError('Empty or HTML response: ' + endpoint)
                target.parent.mkdir(parents=True, exist_ok=True)
                tmp = target.with_suffix(target.suffix + '.tmp')
                tmp.write_bytes(data)
                tmp.replace(target)
                print('Retrieved ' + filename, flush=True)
                return filename, endpoint
            except Exception:
                if attempt == 2:
                    raise
                time.sleep(2 * (attempt + 1))

    with ThreadPoolExecutor(max_workers=2) as pool:
        fetched = dict(pool.map(fetch, ENDPOINTS.items()))
        families = protein_hierarchies(json.loads((out / 'ko00001.json').read_text()))
        fetched.update(pool.map(fetch, [('protein_families/'+bid+'.json', 'get/br:'+bid+'/json') for bid in families]))
        modules = sorted(names(out / 'module_names.tsv'))
        batches = [('module_entries/%04d.txt' % i,
                    'get/' + '+'.join(modules[j:j+10]))
                   for i, j in enumerate(range(0, len(modules), 10))]
        fetched.update(pool.map(fetch, batches))
        # A changed release is not silently combined with the cached files.
        checks = [('end_' + key, value) for key, value in ENDPOINTS.items() if key.startswith('info_')]
        fetched.update(pool.map(fetch, checks))
    for key in ENDPOINTS:
        if key.startswith('info_') and (out / key).read_bytes() != (out / ('end_' + key)).read_bytes():
            raise ValueError('Reference changed during retrieval; use a fresh directory: ' + key)
    entries = module_entries(out)
    if set(entries) != set(modules):
        raise ValueError('Missing or unexpected module records in fetched reference')
    write_json(out / 'manifest.json', dict(schema_version=2, started_utc=started,
        finished_utc=datetime.now(timezone.utc).isoformat(), source='KEGG REST API',
        source_url='https://www.kegg.jp/kegg/rest/',
        files={key: dict(url=BASE+value, sha256=digest(out/key), bytes=(out/key).stat().st_size,
                        retrieved_utc=datetime.fromtimestamp((out/key).stat().st_mtime, timezone.utc).isoformat())
               for key, value in sorted(fetched.items())}))
    verify_reference(out)


def verify_reference(snapshot):
    manifest = json.loads((snapshot / 'manifest.json').read_text())
    for name, info in manifest['files'].items():
        if digest(snapshot / name) != info['sha256']:
            raise ValueError('Reference checksum mismatch: ' + name)
    missing = set(ENDPOINTS) - set(manifest['files'])
    if missing:
        raise ValueError('Incomplete reference manifest: ' + str(missing))
    expected_families = {'protein_families/'+bid+'.json' for bid in
        protein_hierarchies(json.loads((snapshot/'ko00001.json').read_text()))}
    if expected_families - set(manifest['files']):
        raise ValueError('Reference lacks detailed protein-family hierarchies')
    return manifest


def module_entries(snapshot):
    result = {}
    for path in sorted((snapshot / 'module_entries').glob('*.txt')):
        for record in path.read_text().split('///'):
            fields, field = {}, None
            for line in record.strip('\n').splitlines():
                label, value = line[:12].strip(), line[12:].strip()
                if label:
                    field = label
                if field:
                    fields.setdefault(field, []).append(value)
            if 'ENTRY' not in fields:
                continue
            mid = fields['ENTRY'][0].split()[0]
            if mid in result:
                raise ValueError('Duplicate module: ' + mid)
            result[mid] = {key: ' '.join(value) for key, value in fields.items()}
    return result


def hierarchy_rows(tree):
    """Retain every KO occurrence and its full reference path, without guessing."""
    result = []
    def walk(node, ancestors):
        label = node['name']
        if re.match(r'^K\d{5}(?:\s|$)', label):
            if node.get('children'):
                raise ValueError('Unexpected children below KO leaf')
            path = ancestors[1:]  # first item is the hierarchy name
            root = path[0].split()[0] if path else ''
            # Human disease/organismal context is retained but not an
            # environmental phenotype. Unclassified is a separate status.
            scope = ('biological_pathway' if root in BIOLOGICAL_ROOTS else
                     'protein_family' if root == '09180' else
                     'unclassified' if root == '09190' else
                     'organismal_or_disease_context')
            result.append(dict(KO=label[:6], root_id=root, scope=scope,
                root=path[0] if path else '', category=path[1] if len(path)>1 else '',
                term=path[2] if len(path)>2 else '', hierarchy_path=' > '.join(path)))
        else:
            for child in node.get('children', []):
                walk(child, ancestors + [label])
    walk(tree, [])
    return pd.DataFrame(result).drop_duplicates().reset_index(drop=True)


def protein_family_rows(snapshot):
    result = []
    for path in sorted((snapshot/'protein_families').glob('*.json')):
        bid = path.stem
        def walk(node, ancestors):
            label = node['name']
            if re.match(r'^K\d{5}(?:\s|$)', label):
                result.append(dict(KO=label[:6], brite_id=bid,
                    family_path=' > '.join(ancestors[1:]),
                    subfamily=ancestors[-1] if len(ancestors)>1 else '',
                    source='https://www.kegg.jp/brite/'+bid))
            for child in node.get('children', []):
                walk(child, ancestors+[label])
        walk(json.loads(path.read_text()), [])
    return pd.DataFrame(result).drop_duplicates().reset_index(drop=True)


def link_rows(path, id_name):
    frame = pd.read_csv(path, sep='\t', names=['KO', id_name], dtype=str)
    for column in frame:
        frame[column] = frame[column].str.split(':').str[-1]
    if not frame.KO.str.fullmatch(r'K\d{5}').all():
        raise ValueError('Invalid KO link source')
    if id_name == 'pathway_id':
        frame[id_name] = frame[id_name].str.replace(r'^map', 'ko', regex=True)
    return frame.drop_duplicates().reset_index(drop=True)


def description_status(description, current):
    if not description:
        return 'no_description_in_available_references'
    if re.search(r'uncharacterized|unknown function|hypothetical', description, re.I):
        return 'explicitly_uncharacterized_description'
    return 'named_in_current_reference' if current else 'historical_description_only'


def compile_modules(entries):
    definitions = {mid: row.get('DEFINITION', '') for mid, row in entries.items()}
    compiled, audit = {}, []
    for mid, row in sorted(entries.items()):
        expression = definitions[mid]
        try:
            tree = logic.parse(expression, definitions, (mid,))
            blocks = logic.split_top(logic.unwrap(expression.strip()), 'space')
            trees = [logic.parse(block, definitions, (mid,)) for block in blocks]
            if not logic.kos(tree):
                raise ValueError('Definition has no nonoptional KO requirement')
            compiled[mid] = (tree, blocks, trees)
            status, reason = 'evaluated_KO_logic', ''
        except ValueError as exc:
            status, reason = 'not_evaluated_definition', str(exc)
        audit.append(dict(module_id=mid, name=row.get('NAME', ''),
            source_class=row.get('CLASS', ''), definition=expression,
            evaluation_status=status, evaluation_limit=reason,
            source='https://www.kegg.jp/entry/' + mid))
    return compiled, pd.DataFrame(audit)


def summarize_membership(annotations, node_sets):
    """Unique-KO counts; overlapping reference terms never become disjoint bins."""
    rows = []
    for module, members in sorted(node_sets.items()):
        subset = annotations[annotations.KO.isin(members)]
        for (scope, category, term), group in subset.groupby(['scope', 'category', 'term'], dropna=False):
            kos = sorted(set(group.KO))
            rows.append(dict(statistical_module=module, scope=scope, category=category,
                term=term, n_KOs=len(kos), module_n_KOs=len(members),
                fraction_module_KOs=len(kos)/len(members), KOs=';'.join(kos)))
    return pd.DataFrame(rows)


def biological_pathways(hierarchy):
    eligible = hierarchy[hierarchy.scope == 'biological_pathway']
    ids = {'ko'+m for term in eligible.term.unique()
           for m in re.findall(r'\[PATH:ko(\d{5})\]', term)}
    # KEGG global/overview maps are too broad for edge-specific context.
    return {p for p in ids if not p.startswith(('ko011', 'ko012', 'ko013'))}


def evidence_states(compiled, values, n):
    tree, blocks, trees = compiled
    # A standalone optional item cannot inflate coverage merely by being
    # logically dispensable. Keep it in the exported definition, not the
    # required-block denominator.
    required_trees = [block for block in trees if logic.kos(block)]
    states = np.stack([logic.evaluate(block, values, n) for block in required_trees])
    return (logic.evaluate(tree, values, n), (states == 1).sum(axis=0),
            (states == 2).sum(axis=0), len(required_trees))


def export_annotation(args):
    snapshot, out = args.snapshot, args.outdir
    manifest = verify_reference(snapshot)
    if out.exists() and any(out.iterdir()):
        raise ValueError('Use an empty output directory; never overwrite a completed run')
    out.mkdir(parents=True, exist_ok=True)
    matrix = read_matrix(args.matrix, args.samples)
    universe = set(matrix.index)
    if not all(re.fullmatch(r'K\d{5}', k) for k in universe):
        raise ValueError('Invalid matrix KO identifier')
    nodes = pd.read_csv(args.nodes, sep='\t')
    candidates = pd.read_csv(args.candidates, sep='\t')
    edges = pd.read_csv(args.edges, sep='\t')
    old = pd.read_csv(args.dictionary, sep='\t').fillna('')
    for frame in [nodes, candidates, old]:
        if frame.KO.duplicated().any():
            raise ValueError('Duplicated KO identifier in annotation input')
    candidate_set, node_set = set(candidates.KO), set(nodes.KO)
    connected = set(edges.KO1) | set(edges.KO2)
    if not node_set <= connected <= candidate_set <= universe:
        raise ValueError('Inconsistent network, candidate and matrix KO universes')
    if nodes.module_label.isna().any():
        raise ValueError('Missing statistical module label')
    node_sets = {str(label): set(frame.KO) for label, frame in nodes.groupby('module_label')}
    memberships = dict(zip(nodes.KO, nodes.module_label))
    dictionary = names(snapshot / 'ko_names.tsv')
    historical = dict(zip(old.KO, old.KO_function))
    hierarchy = hierarchy_rows(json.loads((snapshot / 'ko00001.json').read_text()))
    hierarchy = hierarchy[hierarchy.KO.isin(universe)].copy()
    write_tsv(hierarchy, out / 'KO_hierarchy_membership.tsv.gz')
    families = protein_family_rows(snapshot)
    families = families[families.KO.isin(universe)].copy()
    write_tsv(families, out/'KO_protein_family_detail.tsv.gz')
    family_summary = families.rename(columns={'brite_id':'category', 'family_path':'term'}).assign(scope='protein_family_detail')
    write_tsv(summarize_membership(family_summary, node_sets), out/'statistical_module_protein_families.tsv')
    pathways = link_rows(snapshot / 'ko_pathway.tsv', 'pathway_id')
    pathways = pathways[pathways.KO.isin(universe)].copy()
    pathway_names = names(snapshot / 'pathway_names.tsv', 'ko')
    pathways['pathway_name'] = pathways.pathway_id.map(pathway_names).fillna('')
    pathways['source'] = 'https://www.kegg.jp/pathway/' + pathways.pathway_id
    write_tsv(pathways, out / 'KO_pathway_membership.tsv.gz')
    reactions = link_rows(snapshot / 'ko_reaction.tsv', 'reaction_id')
    reactions = reactions[reactions.KO.isin(universe)].copy()
    reactions['reference_reaction'] = reactions.reaction_id.map(names(snapshot / 'reaction_names.tsv')).fillna('')
    reactions['source'] = 'https://www.kegg.jp/entry/' + reactions.reaction_id
    write_tsv(reactions, out / 'KO_reaction_membership.tsv.gz')
    entries = module_entries(snapshot)
    if set(entries) != set(names(snapshot / 'module_names.tsv')):
        raise ValueError('Module record inventory differs from reference list')
    compiled, definition_audit = compile_modules(entries)
    write_tsv(definition_audit, out / 'KEGG_module_definitions.tsv')
    info_text = (snapshot/'info_ko.txt').read_text()
    info_match = re.search(r'([\d,]+) entries', info_text)
    info_count = int(info_match.group(1).replace(',', '')) if info_match else None
    write_json(out/'reference_consistency.json', dict(
        ko_names_list_entries=len(dictionary), ko_info_reported_entries=info_count,
        info_count_matches_list=info_count == len(dictionary),
        interpretation='Raw endpoint bytes are frozen independently; an info counter does not overwrite the retrieved KO list.',
        hierarchy_KOs_absent_from_current_name_list=sorted(set(hierarchy.KO)-set(dictionary)),
        protein_family_KOs_absent_from_current_name_list=sorted(set(families.KO)-set(dictionary))))
    # Raw cross-references and resolved logical memberships are distinct.
    module_links = link_rows(snapshot / 'ko_module.tsv', 'module_id')
    module_links = module_links[module_links.KO.isin(universe)].copy()
    module_links['module_name'] = module_links.module_id.map(names(snapshot / 'module_names.tsv')).fillna('')
    write_tsv(module_links, out / 'KO_module_crossreferences.tsv')
    block_rows, logical_rows = [], []
    for mid, (tree, blocks, trees) in compiled.items():
        required = logic.kos(tree)
        for ko in sorted(logic.kos(tree, True) & universe):
            logical_rows.append(dict(KO=ko, module_id=mid, module_name=entries[mid]['NAME'],
                optional_only=ko not in required, source='https://www.kegg.jp/entry/'+mid))
        for i, (block, block_tree) in enumerate(zip(blocks, trees), 1):
            req, all_kos = logic.kos(block_tree), logic.kos(block_tree, True)
            block_rows.append(dict(module_id=mid, block=i, expression=block,
                included_in_required_block_denominator=bool(req),
                nonoptional_KOs=';'.join(sorted(req)), optional_only_KOs=';'.join(sorted(all_kos-req)),
                unavailable_nonoptional_KOs=';'.join(sorted(req-universe)),
                interpretation='Logical requirement; not an inferred temporal or reaction order'))
    logical = pd.DataFrame(logical_rows)
    write_tsv(logical, out / 'KO_module_logical_membership.tsv')
    write_tsv(pd.DataFrame(block_rows), out / 'KEGG_module_blocks.tsv')
    targeted = set()
    rules_path = ROOT / 'config/gene_biogeography/operational_rules.tsv'
    for expr in pd.read_csv(rules_path, sep='\t').expression:
        targeted |= logic.kos(logic.parse(expr), True)
    mapped = set(hierarchy.loc[hierarchy.scope != 'unclassified', 'KO']) | set(families.KO)
    broad_biological = set(hierarchy.loc[hierarchy.scope.isin(['biological_pathway','protein_family']), 'KO']) | set(families.KO)
    hierarchy_counts = hierarchy.groupby('KO').size()
    function_counts = hierarchy[hierarchy.scope != 'unclassified'].groupby('KO').size()
    pcounts, mcounts = pathways.groupby('KO').size(), logical.groupby('KO').size()
    fcounts = families.groupby('KO').size()
    rows = []
    for ko in matrix.index:
        current, previous = dictionary.get(ko, ''), historical.get(ko, '')
        description = current or previous
        rows.append(dict(KO=ko, original_description=previous, current_description=current,
            description_status=description_status(description, bool(current)),
            description_changed=bool(current and previous and current != previous),
            current_KO_entry_available=bool(current),
            n_hierarchy_memberships=int(hierarchy_counts.get(ko, 0)),
            n_classified_memberships=int(function_counts.get(ko, 0)),
            n_protein_family_memberships=int(fcounts.get(ko, 0)),
            n_pathway_memberships=int(pcounts.get(ko, 0)),
            n_logical_module_memberships=int(mcounts.get(ko, 0)),
            original_targeted_catalogue=ko in targeted,
            detected_profiles=int(matrix.loc[ko].gt(0).sum()),
            annotation_records=int(matrix.loc[ko].sum()),
            network_candidate=ko in candidate_set, connected_node=ko in connected,
            statistical_module=memberships.get(ko, ''),
            source='https://www.kegg.jp/entry/'+ko))
    inventory = pd.DataFrame(rows)
    write_tsv(inventory, out / 'KO_annotation_inventory.tsv.gz')
    write_tsv(inventory[inventory.KO.isin(connected)], out / 'network_node_annotation.tsv')
    write_tsv(summarize_membership(hierarchy, node_sets), out / 'statistical_module_functional_composition.tsv')
    # Coverage denominators refer to KOs, not records or duplicated memberships.
    coverage = []
    for scope, kk in [('matrix', universe), ('detected_matrix', set(matrix.index[matrix.gt(0).any(axis=1)])),
                      ('network_candidates', candidate_set), ('connected_network', connected),
                      *[(label, members) for label, members in sorted(node_sets.items())]]:
        sub = inventory[inventory.KO.isin(kk)]
        coverage.append(dict(scope=scope, n_KOs=len(kk),
            current_named_KOs=int(sub.current_KO_entry_available.sum()),
            classified_KOs=len(kk & mapped), pathway_KOs=len(kk & set(pathways.KO)),
            biological_pathway_or_family_KOs=len(kk & broad_biological),
            logical_module_KOs=len(kk & set(logical.KO)),
            targeted_atlas_KOs=len(kk & targeted),
            classified_outside_targeted_atlas=len((kk & mapped)-targeted),
            explicitly_uncharacterized_KOs=int(sub.description_status.eq('explicitly_uncharacterized_description').sum())))
    coverage = pd.DataFrame(coverage)
    write_tsv(coverage, out / 'annotation_coverage.tsv')
    # Per-sample evaluation never combines KOs from different communities.
    community_rows, sample_frames, representation = [], [], []
    for threshold in (1, 2):
        values = {ko: matrix.loc[ko].ge(threshold).to_numpy(np.int8) for ko in universe}
        for mid, obj in compiled.items():
            state, found, unknown, nblocks = evidence_states(obj, values, matrix.shape[1])
            community_rows.append(dict(module_id=mid, name=entries[mid]['NAME'], threshold=threshold,
                n_profiles=matrix.shape[1], rule_met=int((state == 1).sum()),
                rule_not_met=int((state == 0).sum()), rule_unresolved=int((state == 2).sum()),
                logical_blocks=nblocks, median_block_coverage_lower=float(np.median(found/nblocks)),
                median_block_coverage_upper=float(np.median((found+unknown)/nblocks))))
            sample_frames.append(pd.DataFrame(dict(run_accession=matrix.columns, module_id=mid,
                threshold=threshold, rule_state=state, represented_blocks=found,
                unresolved_blocks=unknown, total_blocks=nblocks)))
    for label, members in sorted(node_sets.items()):
        values = {ko: np.array([int(ko in members)], dtype=np.int8) for ko in universe}
        for mid, obj in compiled.items():
            all_kos = logic.kos(obj[0], True)
            if not members & all_kos:
                continue
            state, found, unknown, nblocks = evidence_states(obj, values, 1)
            representation.append(dict(statistical_module=label, module_id=mid, name=entries[mid]['NAME'],
                rule_representation_state=int(state[0]), represented_blocks=int(found[0]),
                unresolved_blocks=int(unknown[0]), total_blocks=nblocks,
                represented_KOs=';'.join(sorted(all_kos & members)),
                measured_KOs_outside_statistical_module=';'.join(sorted((all_kos & universe)-members)),
                unavailable_KOs=';'.join(sorted(all_kos-universe))))
    write_tsv(pd.DataFrame(community_rows), out / 'community_module_rule_summary.tsv')
    write_tsv(pd.concat(sample_frames, ignore_index=True), out / 'community_module_sample_evidence.tsv.gz')
    write_tsv(pd.DataFrame(representation), out / 'statistical_module_biochemical_representation.tsv')
    # Shared reference membership is context, never another statistical edge.
    psets = pathways.groupby('KO').pathway_id.agg(set).to_dict()
    msets = logical.groupby('KO').module_id.agg(set).to_dict()
    eligible_pathways = biological_pathways(hierarchy)
    edge_rows = []
    for edge in edges.itertuples():
        shared_p = psets.get(edge.KO1, set()) & psets.get(edge.KO2, set()) & eligible_pathways
        shared_m = msets.get(edge.KO1, set()) & msets.get(edge.KO2, set())
        edge_rows.append(dict(KO1=edge.KO1, KO2=edge.KO2,
            shared_nonoverview_biological_pathways=';'.join(sorted(shared_p)),
            shared_KEGG_modules=';'.join(sorted(shared_m)),
            interpretation='Shared reference membership; no direct reaction, interaction or direction inferred'))
    write_tsv(pd.DataFrame(edge_rows), out / 'edge_functional_context.tsv')
    # All broad summaries are descriptive: overlapping categories are not an
    # enrichment family and their percentages need not sum to 100.
    categories = []
    for (scope, category), frame in hierarchy.groupby(['scope', 'category']):
        kk = set(frame.KO)
        categories.append(dict(scope=scope, category=category, matrix_KOs=len(kk),
            detected_KOs=len(kk & set(matrix.index[matrix.gt(0).any(axis=1)])),
            candidate_KOs=len(kk & candidate_set), connected_KOs=len(kk & connected),
            outside_targeted_atlas_KOs=len(kk-targeted)))
    write_tsv(pd.DataFrame(categories), out / 'functional_category_coverage.tsv')
    summary = dict(n_profiles=matrix.shape[1], n_matrix_KOs=len(universe),
        reference_KOs=len(dictionary), reference_KEGG_modules=len(entries),
        reference_protein_hierarchies=len(protein_hierarchies(json.loads((snapshot/'ko00001.json').read_text()))),
        evaluable_module_definitions=len(compiled),
        unresolved_module_definitions=len(entries)-len(compiled),
        biological_pathway_terms=len(eligible_pathways),
        KO_description_changes=int(inventory.description_changed.sum()),
        current_unlisted_matrix_KOs=sorted(universe-set(dictionary)),
        node_count=len(connected), edge_count=len(edges),
        edges_sharing_nonoverview_biological_pathway=sum(bool(r['shared_nonoverview_biological_pathways']) for r in edge_rows),
        edges_sharing_KEGG_module=sum(bool(r['shared_KEGG_modules']) for r in edge_rows),
        scope='Descriptive expansion; no new edges, partitions, enrichments, fluxes or sequence assignments',
        coverage=coverage.to_dict('records'))
    write_json(out / 'summary.json', summary)
    provenance = dict(created_utc=datetime.now(timezone.utc).isoformat(),
        inputs={key: dict(path=str(getattr(args,key).resolve()), sha256=digest(getattr(args,key)))
                for key in ['matrix', 'samples', 'nodes', 'candidates', 'edges', 'dictionary']},
        reference_manifest_sha256=digest(snapshot/'manifest.json'), reference=manifest,
        source_sha256={str(path.relative_to(ROOT)): digest(path) for path in [Path(__file__),
            ROOT/'scripts/12_gene_environment/template/evaluate.py',
            ROOT/'scripts/10_functional_network/lib/biogeochemical_context.py', rules_path]},
        runtime=dict(python=sys.version, numpy=np.__version__, pandas=pd.__version__))
    write_json(out / 'provenance.json', provenance)
    write_json(out/'decision_record.json', dict(
        D0=dict(status='viable', decision='Community DNA WGS KO annotation-record matrix; original assignments retained'),
        D1=dict(status='viable', decision='Describe the entire repertoire and biological context of existing network modules'),
        D2=dict(status='viable', decision='Same accession-selected cohort; no new sample-level inferential tests'),
        D3=dict(status='conditional', decision='Detection thresholds one/two records are operational; no sequence revalidation or calibrated abundance'),
        D4=dict(status='viable', decision='Offline identifier mapping to a hash-verified reference; no sequencing or annotation pipeline rerun'),
        D5=dict(status='conditional', decision='Many-to-many assignments, explicit unknowns, module rules and detailed families; no inferred carriers or flux'),
        D6=dict(status='viable', decision='Descriptive coverage and representation; no new enrichment P values'),
        D7=dict(status='viable', decision='Original/current descriptions kept; two detection thresholds; unsupported definitions listed'),
        D8=dict(status='viable', decision='Frozen references, source/input checksums and full tables; code-only repository'),
        not_estimable=['Activity or flux', 'Same-genome pathway completion', 'Taxonomic carriers from community KO counts alone',
                       'Resistance phenotype or pathogenicity', 'Significant functional specialization from descriptive counts alone']))
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    fetch = sub.add_parser('fetch-reference')
    fetch.add_argument('--outdir', type=Path, required=True)
    annotate = sub.add_parser('annotate')
    for flag in ['matrix', 'samples', 'nodes', 'candidates', 'edges', 'dictionary', 'snapshot', 'outdir']:
        annotate.add_argument('--' + flag, type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'fetch-reference':
        fetch_reference(args.outdir)
    else:
        export_annotation(args)


if __name__ == '__main__':
    main()
