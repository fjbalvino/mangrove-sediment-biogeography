#!/usr/bin/env python3
"""Independently check cohort, consensus edges and reported network source data."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--outdir', type=Path, required=True)
    args = p.parse_args()
    root = args.outdir.resolve()
    checks = {}

    def check(name, value):
        checks[name] = bool(value)
        if not value:
            raise AssertionError(name)

    def run(n):
        paths = {q.read_text().strip() for q in (root / 'resultados').glob('LATEST_' + n + '_*.txt')}
        if len(paths) != 1:
            raise ValueError('Ambiguous run ' + n)
        r = Path(paths.pop())
        check('stage_' + n + '_completed', (r / 'STATUS.txt').read_text().strip() == 'PASS')
        return r / 'tables'

    b, c, e, f = [run(n) for n in ['156', '171', '173', '174']]
    meta = pd.read_csv(b / '156_metadata_aligned.tsv', sep='\t')
    cohort = pd.read_csv(root / 'inputs/metadata_191.tsv', sep='\t')
    check('191_unique_samples', len(meta) == meta.run_accession.nunique() == 191)
    check('exact_cohort', set(meta.run_accession) == set(cohort.run_accession))
    check('no_air_controls', not set(meta.run_accession) & {'SRR31095011', 'SRR31095534', 'SRR31095535'})
    dimensions = pd.read_csv(b / '156_expected_dimensions_validation.tsv', sep='\t')
    check('input_dimensions_and_prevalence_classes', dimensions['pass'].all())
    transform = pd.read_csv(b / '156_transformation_validation.tsv', sep='\t')
    check('CLR_transform_checks', transform['pass'].all())
    candidates = pd.read_csv(c / '171_candidate_KO_selection.tsv', sep='\t')
    edges = pd.read_csv(c / '171_consensus_edges.tsv.gz', sep='\t')
    membership = pd.read_csv(c / '171_bootstrap_edge_membership.tsv.gz', sep='\t')
    design = pd.read_csv(c / '171_bootstrap_design.tsv.gz', sep='\t')
    check('1500_unique_candidates', len(candidates) == candidates.KO.nunique() == 1500)
    check('500_subsamples_of_81', len(design.groupby('bootstrap_id')) == 500 and
          design.groupby('bootstrap_id').run_accession.nunique().eq(81).all())
    check('at_most_four_per_stratum', design.groupby(['bootstrap_id', 'BALANCE_STRATUM']).size().le(4).all())
    check('2000_edges_per_subsample', membership.groupby('bootstrap_id').size().eq(2000).all())
    check('unique_edges_no_loops', not edges.duplicated(['KO1', 'KO2']).any() and edges.KO1.ne(edges.KO2).all())
    check('all_endpoints_are_candidates', set(edges.KO1) | set(edges.KO2) <= set(candidates.KO))
    check('consensus_thresholds', edges.edge_stability.ge(.8).all() and
          edges.sign_consistency.eq(1).all() and edges.mean_abs_weight.ge(.95).all())
    observed = membership.groupby(['KO1', 'KO2']).agg(n=('bootstrap_id', 'nunique'),
                                                     mean_abs=('abs_weight', 'mean'))
    paired = edges.set_index(['KO1', 'KO2']).join(observed)
    check('independent_edge_frequency', np.allclose(paired.edge_stability, paired.n / 500, atol=1e-12))
    check('independent_edge_weight', np.allclose(paired.mean_abs_weight, paired.mean_abs, atol=1e-12))
    nodes = pd.read_csv(e / '173_network_nodes_annotated.tsv', sep='\t')
    classes = pd.read_csv(root / 'inputs/KO_core_classes.tsv', sep='\t').set_index('KO')
    check('node_classes_match_atlas', nodes.functional_core_class.to_numpy().tolist() ==
          classes.loc[nodes.KO, 'functional_core_class'].to_numpy().tolist())
    check('major_modules_have_at_least_10_KOs', nodes.groupby('module_label').size().ge(10).all())
    enrich = pd.read_csv(f / '174_module_functional_category_enrichment_permutation.tsv', sep='\t')
    check('complete_seven_category_test_family', len(enrich)==nodes.module_label.nunique()*7 and enrich.groupby('module_label').category.nunique().eq(7).all())
    untestable = enrich[enrich.p_maxT_FWER.isna()]
    check('undefined_maxT_is_constant_zero_case', untestable.null_sd.eq(0).all() and untestable.observed_n.eq(0).all() and untestable.expected_n.eq(0).all())
    pv = enrich.p_empirical.to_numpy()
    ix = np.argsort(pv)
    raw = pv[ix] * len(pv) * np.sum(1 / np.arange(1, len(pv) + 1)) / np.arange(1, len(pv) + 1)
    adjusted = np.minimum(1, np.minimum.accumulate(raw[::-1])[::-1])
    check('BY_adjustment_recomputed', np.allclose(enrich.q_BY.to_numpy()[ix], adjusted))
    supported = pd.read_csv(f / '174_supported_module_functional_specialization.tsv', sep='\t')
    check('supported_categories_pass_both_corrections', supported.q_BY.le(.05).all() and supported.p_maxT_FWER.le(.05).all())
    annotation = pd.read_csv(f / '174_annotation_source_sensitivity_KOs.tsv', sep='\t')
    combined = pd.read_csv(f / '174_combined_annotation_enrichment.tsv', sep='\t')
    annotation_support = pd.read_csv(f / '174_annotation_source_support.tsv', sep='\t')
    check('annotation_sensitivity_exact_node_set',
          len(annotation) == annotation.KO.nunique() == len(nodes) and set(annotation.KO) == set(nodes.KO))
    check('annotation_change_flags', annotation.category_changed.equals(
        annotation.functional_category.ne(annotation.combined_category)))
    for label, table, column in [('primary', enrich, 'functional_category'),
                                  ('combined', combined, 'combined_category')]:
        observed = annotation.groupby(['module_label', column]).size()
        recomputed = [observed.get((row.module_label, row.category), 0) for row in table.itertuples()]
        check(label + '_category_counts_recomputed', np.array_equal(table.observed_n, recomputed))
        check(label + '_fold_enrichment_recomputed', np.allclose(
            table.fold_enrichment, (table.observed_n + .5) / (table.expected_n + .5)))
        pv = table.p_empirical.to_numpy()
        order = np.argsort(pv)
        raw = pv[order] * len(pv) * np.sum(1 / np.arange(1, len(pv) + 1)) / np.arange(1, len(pv) + 1)
        check(label + '_BY_family_recomputed', np.allclose(
            table.q_BY.to_numpy()[order], np.minimum(1, np.minimum.accumulate(raw[::-1])[::-1])))
    check('annotation_support_requires_both_sources',
          annotation_support.supported_both_sources.equals(
              annotation_support.q_BY.le(.05) & annotation_support.p_maxT_FWER.le(.05) &
              annotation_support.combined_q_BY.le(.05) & annotation_support.combined_maxT_P.le(.05)))
    perturb = pd.read_csv(f / '174_edge_perturbation_summary.tsv', sep='\t')
    check('500_completed_perturbations', perturb.perturbation.nunique() == 500 and np.isfinite(perturb.ARI_vs_primary_modules).all())
    mixing = pd.read_csv(f / '174_prevalence_degree_preserving_assortativity.tsv', sep='\t')
    check('999_valid_rewirings', mixing.n_valid_rewires.eq(999).all())
    algo = pd.read_csv(f / '174_community_algorithm_sensitivity.tsv', sep='\t')
    check('four_successful_community_algorithms', len(algo) == 4 and np.isfinite(algo.ARI_vs_primary_major_modules).all())
    threshold = pd.read_csv(f / '174_network_threshold_sensitivity.tsv', sep='\t')
    check('30_completed_threshold_combinations', len(threshold) == 30 and np.isfinite(threshold.ARI_vs_primary_major_modules).all())
    check('identical_threshold_graphs_use_identical_partitions', threshold.groupby('stability_min').ARI_vs_primary_major_modules.nunique().eq(1).all())
    check('primary_partition_replayed', np.isclose(algo.loc[algo.method.eq('leiden'), 'ARI_vs_primary_major_modules'].iloc[0], 1))
    result = {'status': 'PASS', 'checks': checks, 'n_checks': len(checks)}
    (root / 'network_independent_verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
