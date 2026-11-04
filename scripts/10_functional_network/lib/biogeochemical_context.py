"""Overlay frozen biochemical evidence on an independently inferred KO network.

Membership and edges are never changed by biochemical annotation. Logical blocks
are not automatically reaction steps. Ordered steps come only from the curated
reference-chain file. Community detection evidence is evaluated per sample.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location(
    'atlas_rule_logic', ROOT / 'scripts/12_gene_environment/template/evaluate.py')
logic = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(logic)


def read_matrix(matrix, samples):
    ids = [x.strip() for x in Path(samples).read_text().splitlines() if x.strip()]
    x = pd.read_csv(matrix, sep='\t', index_col=0)
    if not ids or len(ids) != len(set(ids)) or x.index.has_duplicates or x.columns.has_duplicates:
        raise ValueError('Empty or duplicate sample/KO identifiers')
    if not set(ids).issubset(x.columns):
        raise ValueError('Selected accessions are missing from the matrix')
    x = x.loc[:, ids]
    a = x.to_numpy(float)
    if not np.isfinite(a).all() or (a < 0).any() or (a != np.floor(a)).any() or (a.sum(axis=0) <= 0).any():
        raise ValueError('Expected finite nonnegative predicted-protein record counts and positive profile totals')
    return x


def get_blocks(expression):
    return logic.split_top(logic.unwrap(expression.strip()), 'space')


def module_block_state(expression, members, matrix_kos):
    """Evidence of representation IN a statistical module, never sample absence."""
    required = logic.kos(logic.parse(expression))
    values = {k: np.array([int(k in members)], dtype=np.int8) for k in required & matrix_kos}
    state = int(logic.evaluate(logic.parse(expression), values, 1)[0])
    if state == 1:
        return 'logical_block_represented'
    if required & members:
        return 'partial_representation'
    if state == 2:
        return 'not_resolvable_from_matrix'
    return 'not_represented_in_module'


def chain_expression(step, rules):
    expr = str(rules.loc[step['rule'], 'expression'])
    if 'block' in step:
        blocks = get_blocks(expr)
        if not 1 <= step['block'] <= len(blocks):
            raise ValueError('Reference chain block outside frozen rule: ' + step['rule'])
        return blocks[step['block'] - 1]
    return expr


def export_context(matrix, nodes, candidates, edges, config, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rules = pd.read_csv(config / 'operational_rules.tsv', sep='\t').set_index('rule_id')
    processes = pd.read_csv(config / 'curated_processes.tsv', sep='\t').set_index('id')
    labels = pd.read_csv(config / 'KO_labels.tsv', sep='\t').set_index('KO')
    chains = json.loads((config / 'network_reference_chains.json').read_text())['chains']
    if rules.index.has_duplicates or processes.index.has_duplicates or nodes.KO.duplicated().any():
        raise ValueError('Nonunique rule, process or network node identifier')
    if not set(nodes.KO).issubset(matrix.index) or not set(nodes.KO).issubset(candidates):
        raise ValueError('Network nodes must be observed matrix rows and network candidates')
    if not set(rules.process_id).issubset(processes.index):
        raise ValueError('Unmapped rule process')
    connected = set(edges.KO1) | set(edges.KO2)
    membership = dict(zip(nodes.KO, nodes.module_label))
    annotation, block_rows, coverage = [], [], []
    for rule_id, rule in rules.iterrows():
        process = processes.loc[rule.process_id]
        tree = logic.parse(rule.expression)
        required = logic.kos(tree)
        all_kos = logic.kos(tree, include_optional=True)
        for k in sorted(all_kos):
            annotation.append(dict(KO=k, name=labels.loc[k, 'name'] if k in labels.index else '',
                cycle=process.cycle, process_id=rule.process_id, process=process.process,
                rule_id=rule_id, evidence_kind=rule.evidence_kind,
                optional_only_in_rule=k not in required,
                matrix_row_available=k in matrix.index, network_candidate=k in candidates,
                connected_in_reference=k in connected, statistical_module=membership.get(k, ''),
                annotation_limit=rule['limit'], process_limit=process.interpretation_and_limit,
                source=rule.source))
        for i, block in enumerate(get_blocks(rule.expression), 1):
            kk = logic.kos(logic.parse(block))
            block_rows.append(dict(rule_id=rule_id, process_id=rule.process_id, cycle=process.cycle,
                block=i, expression=block, required_or_alternative_KOs=';'.join(sorted(kk)),
                unavailable_KOs=';'.join(sorted(kk - set(matrix.index))),
                interpretation='Logical requirement; reaction order is not inferred from this field'))
            for module, n in nodes.groupby('module_label', sort=True):
                members = set(n.KO)
                coverage.append(dict(module=module, cycle=process.cycle, process_id=rule.process_id,
                    rule_id=rule_id, block=i, expression=block,
                    state=module_block_state(block, members, set(matrix.index)),
                    represented_KOs=';'.join(sorted(kk & members)),
                    matrix_KOs_outside_module=';'.join(sorted((kk & set(matrix.index)) - members)),
                    unavailable_KOs=';'.join(sorted(kk - set(matrix.index)))))
    annotation = pd.DataFrame(annotation)
    annotation.to_csv(out / 'KO_process_context.tsv', sep='\t', index=False)
    pd.DataFrame(block_rows).to_csv(out / 'logical_block_definitions.tsv', sep='\t', index=False)
    pd.DataFrame(coverage).to_csv(out / 'module_logical_block_representation.tsv', sep='\t', index=False)
    node_rows = []
    for row in nodes.itertuples():
        a = annotation[annotation.KO == row.KO]
        node_rows.append(dict(KO=row.KO, module=row.module_label,
            curated_cycles=';'.join(sorted(set(a.cycle))),
            curated_processes=';'.join(sorted(set(a.process_id))),
            context_status='mapped_to_curated_process' if len(a) else 'outside_targeted_catalogue',
            original_function=getattr(row, 'KO_function', ''),
            interpretation='Many-to-many membership; outside this targeted catalogue does not mean unknown function'))
    pd.DataFrame(node_rows).to_csv(out / 'network_nodes_process_context.tsv', sep='\t', index=False)
    attrition = []
    for cycle, a in annotation.groupby('cycle', sort=False):
        kk = set(a.KO)
        attrition.append(dict(cycle=cycle, queried_KOs=len(kk), matrix_KOs=len(kk & set(matrix.index)),
            candidate_KOs=len(kk & set(candidates)), connected_KOs=len(kk & connected),
            major_module_KOs=len(kk & set(nodes.KO))))
    pd.DataFrame(attrition).to_csv(out / 'cycle_selection_coverage.tsv', sep='\t', index=False)
    summaries, sample_rows, chain_membership = [], [], []
    for threshold in (1, 2):
        values = {k: (matrix.loc[k].to_numpy() >= threshold).astype(np.int8)
                  for k in set(annotation.KO) & set(matrix.index)}
        for chain in chains:
            states = []
            for role, steps in [('chain_step', chain['steps']), ('associated_machinery', chain.get('accessory_steps', []))]:
                for i, step in enumerate(steps, 1):
                    expression = chain_expression(step, rules)
                    state = logic.evaluate(logic.parse(expression), values, matrix.shape[1])
                    if role == 'chain_step':
                        states.append(state)
                    summaries.append(dict(chain=chain['id'], cycle=chain['cycle'], role=role,
                        step=i, label=step['label'], threshold=threshold, expression=expression,
                        detected=int((state == 1).sum()), not_met=int((state == 0).sum()),
                        unknown=int((state == 2).sum()), total=matrix.shape[1],
                        reference_direction=chain['direction'], source=chain['source'], limit=chain['limit']))
                    if threshold == 1:
                        kk = logic.kos(logic.parse(expression))
                        for k in sorted(kk):
                            chain_membership.append(dict(chain=chain['id'], role=role, step=i,
                                label=step['label'], KO=k, matrix_available=k in matrix.index,
                                network_candidate=k in candidates, connected=k in connected,
                                module=membership.get(k, ''), reference_direction=chain['direction']))
            s = np.stack(states)
            all_steps = np.where((s == 0).any(axis=0), 0, np.where((s == 2).any(axis=0), 2, 1))
            summaries.append(dict(chain=chain['id'], cycle=chain['cycle'], role='all_selected_chain_steps',
                step=0, label=chain['label'], threshold=threshold, expression='',
                detected=int((all_steps == 1).sum()), not_met=int((all_steps == 0).sum()),
                unknown=int((all_steps == 2).sum()), total=matrix.shape[1],
                reference_direction=chain['direction'], source=chain['source'], limit=chain['limit']))
            for j, accession in enumerate(matrix.columns):
                sample_rows.append(dict(run_accession=accession, chain=chain['id'], threshold=threshold,
                    selected_chain_state=int(all_steps[j]), detected_steps=int((s[:, j] == 1).sum()),
                    unknown_steps=int((s[:, j] == 2).sum()), total_selected_steps=len(states)))
    pd.DataFrame(summaries).to_csv(out / 'reference_chain_evidence.tsv', sep='\t', index=False)
    pd.DataFrame(sample_rows).to_csv(out / 'reference_chain_sample_evidence.tsv.gz', sep='\t', index=False)
    pd.DataFrame(chain_membership).to_csv(out / 'reference_chain_network_membership.tsv', sep='\t', index=False)
    process_scope = processes.reset_index().copy()
    process_scope['n_operational_rules'] = process_scope.id.map(rules.groupby('process_id').size()).fillna(0).astype(int)
    process_scope.to_csv(out / 'process_scope.tsv', sep='\t', index=False)
    report = dict(n_profiles=matrix.shape[1], n_processes=len(processes), n_rules=len(rules),
        n_queried_KOs=annotation.KO.nunique(), n_matrix_KOs=len(set(annotation.KO) & set(matrix.index)),
        n_statistical_module_nodes=len(nodes), n_module_nodes_mapped=len(set(annotation.KO) & set(nodes.KO)),
        n_reference_chains=len(chains), network_membership_modified=False,
        inference='Descriptive biochemical overlay; no new enrichment test or calibrated edge inference',
        caveats=['Logical blocks are not necessarily successive reactions',
                 'Outside-module is not missing from matrix or biological absence',
                 'Shared KOs cannot independently confirm multiple pathways',
                 'Community co-detection does not establish genomic linkage, expression, direction or flux'],
        curation_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in
                        [config / 'operational_rules.tsv', config / 'curated_processes.tsv',
                         config / 'KO_labels.tsv', config / 'network_reference_chains.json']})
    (out / 'context_audit.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['matrix', 'samples', 'nodes', 'candidates', 'edges', 'outdir']:
        p.add_argument('--' + key, type=Path, required=True)
    p.add_argument('--config', type=Path, default=ROOT / 'config/gene_biogeography')
    a = p.parse_args()
    print(json.dumps(export_context(read_matrix(a.matrix, a.samples), pd.read_csv(a.nodes, sep='\t'),
        set(pd.read_csv(a.candidates, sep='\t').KO), pd.read_csv(a.edges, sep='\t'), a.config, a.outdir), indent=2))


if __name__ == '__main__':
    main()
