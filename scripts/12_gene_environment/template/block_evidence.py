"""Report rule blocks for all cycles; never change strict rule decisions."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from evaluate import block_evidence


def main():
    root = Path(__file__).resolve().parent
    matrix = pd.read_csv(root/'inputs/functional_KO_counts_original.tsv.gz', sep='\t', index_col=0)
    ids = (root/'inputs/functional_accessions.txt').read_text().splitlines()
    matrix = matrix.loc[:, ids]
    rules = pd.read_csv(root/'inputs/operational_rules.tsv', sep='\t')
    rows, definitions = [], []
    for threshold in (1, 2):
        values = {k: (row.to_numpy() >= threshold).astype(np.int8) for k, row in matrix.iterrows()}
        for rule in rules.itertuples():
            result = block_evidence(rule.expression, values, len(ids))
            if threshold == 1:
                definitions.extend(dict(rule_id=rule.rule_id, block=i+1, expression=b)
                                   for i, b in enumerate(result['blocks']))
            for i, run in enumerate(ids):
                rows.append(dict(run_accession=run, rule_id=rule.rule_id, threshold=threshold,
                    n_blocks=result['n_blocks'],
                    **{k: result[k][i] for k in ('n_detected','n_not_met','n_unknown','coverage_lower','coverage_upper')}))
    out = root/'outputs'
    pd.DataFrame(rows).to_csv(out/'rule_block_evidence.tsv.gz', sep='\t', index=False)
    pd.DataFrame(definitions).to_csv(out/'rule_block_definitions.tsv', sep='\t', index=False)
    (out/'rule_block_interpretation.json').write_text(json.dumps({
        'unit': 'Top-level logical blocks in the frozen expression; not necessarily individual reactions',
        'bounds': 'Detected/total to (detected+unresolvable)/total; logical bounds, not confidence intervals',
        'complex': 'Mandatory plus-linked subunits are evaluated jointly within a block',
        'zero': 'Rule block not met in observed annotations; not biological absence',
        'one': 'Every block detected does not establish same-genome completeness, direction or activity',
        'thresholds': [1,2], 'strict_rule_outputs_changed': False
    }, indent=2)+'\n')


if __name__ == '__main__':
    main()
