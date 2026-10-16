#!/usr/bin/env python3
"""Regenerate the canonical Figure 5 manuscript sources from verified tables."""
import argparse,json,re
from pathlib import Path
import numpy as np
import pandas as pd


def mdtable(df):
    return '| '+' | '.join(map(str,df.columns))+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in df.itertuples(index=False,name=None))

def fmt(x,d=3):return 'NA' if pd.isna(x) else f'{x:.{d}f}'

def portable_provenance(text):
    """Retain run identifiers without publishing host-specific execution paths."""
    text = re.sub(r'/workspace/[^\s"\t]+?/resultados/', 'source_runs/', text)
    return re.sub(r'/workspace/[^\s"\t]+?/inputs/', 'source_inputs/', text)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--analysis',type=Path,required=True);p.add_argument('--study',type=Path,required=True)
    p.add_argument('--parameter-audit',type=Path,required=True);p.add_argument('--repository',type=Path,default=Path(__file__).resolve().parents[1])
    a=p.parse_args();runs={}
    for num in ['144','156','171','172','173','174']:
        # Delivered tables use portable stage directories; fresh runs use LATEST pointers.
        if (a.analysis/num/'STATUS.txt').is_file():
            runs[num]=a.analysis/num
            assert (runs[num]/'STATUS.txt').read_text().strip()=='PASS'
            continue
        targets={v.read_text().strip() for v in (a.analysis/'resultados').glob('LATEST_'+num+'_*')}
        assert len(targets)==1
        runs[num]=Path(targets.pop());assert (runs[num]/'STATUS.txt').read_text().strip()=='PASS'
    read=lambda n,s:pd.read_csv(runs[n]/'tables'/s,sep='\t')
    nodes=read('173','173_network_nodes_annotated.tsv');edges=read('171','171_consensus_edges.tsv.gz')
    enrich=read('174','174_module_functional_category_enrichment_permutation.tsv')
    supported=read('174','174_supported_module_functional_specialization.tsv')
    annot=read('174','174_annotation_source_support.tsv')
    stability=read('174','174_module_membership_stability.tsv')
    perturb=read('174','174_edge_perturbation_summary.tsv')
    threshold=read('174','174_network_threshold_sensitivity.tsv')
    algorithms=read('174','174_community_algorithm_sensitivity.tsv')
    selection=read('174','174_backbone_selection_matched_permutation.tsv')
    study=pd.read_csv(a.study/'study_exclusion_results.tsv',sep='\t')
    study_scores=pd.read_csv(a.study/'study_exclusion_module_scores.tsv',sep='\t')
    known=study.query("unit_type=='BioProject'")
    partition=pd.read_csv(a.parameter_audit/'partition500/results/partition_grid_summary.tsv',sep='\t')
    audit=json.loads((a.parameter_audit/'results/audit_summary.json').read_text())
    metrics=read('171','171_network_global_metrics.tsv').set_index('metric').value
    near=int(nodes.functional_core_class.eq('near_global').sum());n=len(nodes)
    displayed=int((edges.KO1.isin(nodes.KO)&edges.KO2.isin(nodes.KO)).sum())
    # The manuscript wording is tied to this frozen cohort and evaluated settings.
    expected={('M01','Information processing'),('M03','Carbon / carbohydrate metabolism'),('M06','Transport / membrane systems')}
    assert set(zip(supported.module_label,supported.category))==expected,'Enrichment changed; review scientific wording before rebuilding.'
    assert annot.query('primary_supported and not combined_supported').module_label.to_list()==['M03']
    f=supported.set_index('module_label').fold_enrichment
    title='Widely distributed KOs form functional co-abundance modules'
    main_text=(f'Using gene annotations and CLR-transformed KO counts from 191 metagenomes, we identified groups of functions with similar abundance patterns across mangrove sediments. '
      f'The network contained {len(edges):,} associations among {int(metrics.n_connected_KOs)} connected KOs. Six major groups, termed modules, included {n} KOs, of which {near} ({100*near/n:.1f}%) were near-global (Fig. 5a; Supplementary Table S14). '
      f'M01 and M06 were enriched in information-processing and transport/membrane functions, respectively ({f.M01:.2f}- and {f.M06:.2f}-fold), with support under both annotation treatments. '
      f'Enrichment of carbon/carbohydrate functions in M03 ({f.M03:.2f}-fold) depended on the annotation source (Fig. 5b; Supplementary Tables S15–S22).\n\n'
      'KOs within the same module generally remained grouped when some network links were removed. To determine whether these functional groups also depended on particular sources of samples, we reconstructed the network after temporarily omitting all samples from each BioProject in turn. '
      'The largest change occurred when six metagenomes belonging to one BioProject were omitted: the number of major modules decreased from six to four, and the associations and module composition changed substantially (Fig. 5c; Supplementary Tables S23, S26 and S27). '
      'All 191 metagenomes were retained in the primary analysis. Thus, the modules describe groups of microbial functions whose abundances vary together in the sampled sediments, but their organization is sensitive to the samples included. '
      'These associations may reflect genes carried by the same microorganisms or similar environmental responses among their carriers; they do not establish direct biochemical interactions or metabolic activity.')
    caption=(f'Figure 5 | Functional co-abundance modules are dominated by widely distributed KOs. '
      f'(a) Six major modules inferred from 191 metagenomes and 1,500 candidate KOs. The displayed subnetwork contains {n} KOs and {displayed} associations; the full network contains {int(metrics.n_connected_KOs)} connected KOs and {len(edges):,} associations. '
      'Colours identify modules; labels give enriched functional categories and total KO counts. Diverse functions denotes no supported category enrichment. Squares, circles and triangles represent absolute-core, near-global and realm-core KOs, respectively. '
      '(b) Primary category enrichment supported by both Benjamini–Yekutieli q ≤ 0.05 and maxT family-wise P ≤ 0.05 across 4,999 blocked permutations and 42 module–category comparisons. Point sizes and labels give category KO counts; the dashed line denotes no enrichment. The open point and asterisk mark support sensitive to adding partial pathway annotations. '
      '(c) Membership sensitivity with 50-start Leiden re-estimation. Circles show mean within-module co-assignment across 500 edge perturbations; diamonds show the median across 43 BioProject omissions, with horizontal lines spanning the observed range. Edges were retained in perturbations according to their consensus selection frequencies. Higher co-assignment indicates that pairs remain together but does not penalize mergers. These ranges are not confidence intervals.')

    docs=a.repository/'docs';(docs/'FUNCTIONAL_NETWORK_RESULTS.md').write_text('**'+title+'**\n\n'+main_text+'\n\n'+caption+'\n')
    comp=pd.crosstab(nodes.module_label,nodes.functional_core_class).reindex(columns=['absolute_global_core','near_global','realm_core_only']).reset_index()
    comp.insert(1,'KOs',comp.iloc[:,1:].sum(axis=1));comp.columns=['Module','KOs','Absolute core','Near-global','Realm-core']
    sub=['# Supplementary results for the functional network',
      f'The network from 191 metagenomes contained {len(edges):,} associations among {int(metrics.n_connected_KOs)} connected KOs. Six major modules contained {n} KOs and {displayed} associations. Information-processing and transport/membrane enrichment retained support across annotation treatments; carbon/carbohydrate enrichment did not. Sensitivity analyses assess distinct sources of uncertainty and do not establish a unique biological partition.',
      '## Module composition',
      f'Major modules contained at least 10 KOs and were labelled by decreasing size, then internal edge count. Their persistence classes comprised 9 absolute-core, {near} near-global and 33 realm-core KOs. No peripheral KO occurred in these modules. The complete connected network comprised {int(metrics.n_components_among_connected_KOs)} components; the largest contained {int(metrics.largest_connected_component_size)} KOs. There were {int(metrics.n_community_modules_connected)} communities among connected nodes, including the six major modules. Module identifiers are size ranks, not conserved biological identities across parameter settings.',
      mdtable(comp),
      '## Functional category support and annotation sensitivity',
      'Primary categories were assigned from the frozen KO-description dictionary using the existing keyword rules and their deterministic tie order. Adding the incomplete historical KEGG MODULE text was evaluated with exactly the same permutation labels. Neither treatment is an independent biological validation. Fold enrichment is the observed count plus 0.5 divided by the mean permuted count plus 0.5. A supported label identifies an enriched category, not the function of every KO in the module.']
    s=annot.query('primary_supported').copy()
    labelmap={'Information processing':'Information processing','Carbon / carbohydrate metabolism':'Carbon/carbohydrate','Transport / membrane systems':'Transport/membrane'}
    t=pd.DataFrame({'Module':s.module_label,'Category':s.category.map(labelmap),'KOs':s.observed_n,'Fold':s.fold_enrichment.map(lambda x:fmt(x,2)),'BY q':s.q_BY.map(lambda x:fmt(x,4)),'maxT P':s.p_maxT_FWER.map(lambda x:fmt(x,4))})
    sub += [mdtable(t), 'The carbon/carbohydrate category in M03 contained five KOs: K01897, K21071, K00700, K08483 and K01912. The combined annotation assigned seven KOs to this category, but its support declined (BY q = '+fmt(s.loc[s.module_label.eq('M03'),'combined_q_BY'].iloc[0],4)+'; maxT P = '+fmt(s.loc[s.module_label.eq('M03'),'combined_maxT_P'].iloc[0],4)+'). Increased category counts alone therefore did not imply stronger conditional enrichment. No module passed both criteria for the anaerobic redox/energy category. M02, M04 and M05 are labelled Diverse functions because no primary category passed both criteria.',
      f'The partial pathway snapshot supplied annotations for {int(nodes.KEGG_module_ids.notna().sum())} of the {n} major-module KOs. Unmapped KOs cannot be interpreted as lacking pathway membership. KEGG MODULE enrichment was not tested because a complete, versioned mapping of the analysis universe was unavailable.',
      '## Sampling and network inference',
      'The accession manifest retained 191 profiles, including TX1_3 and all 23 Guadeloupe profiles, and excluded air controls SRR31095011, SRR31095534 and SRR31095535. This is the functional atlas cohort, not the balanced subset used for regional differentiation in Figure 6. Counts and metadata were joined by run_accession, and atlas persistence assignments were verified independently.',
      'Of 17,482 detected KOs, 13,308 met the association filter of occurrence in at least 10 profiles and at least 10 total assignments. A pseudocount of 1 was added before CLR transformation within this fixed reference universe. The network filter required occurrence in at least 20 profiles and at least 50 assignments; 10,614 KOs qualified. The 5,000 highest CLR-variance KOs were retained, followed by the top 1,500 network candidates, with KO identifiers resolving ties. Persistence labels and enrichment results were not used to select candidates.',
      'Five hundred subsamples without replacement were generated with seed 123, capped at four profiles per marine realm × country × locality stratum. Each subsample contained 81 profiles from 26 strata: 35 Central Indo-Pacific, 23 Tropical Atlantic, 20 Western Indo-Pacific and three Tropical Eastern Pacific profiles. Thirty-three profiles occurred in every subsample. This design limits local replication but neither equalizes realms nor removes study effects.',
      'For each subsample, proportionality was calculated as 2 cov(x,y)/(var(x)+var(y)) in CLR space. The 2,000 strongest absolute associations were retained. Consensus required selection frequency ≥ 0.80, complete sign consistency and mean absolute proportionality ≥ 0.95. Mean weights were calculated over subsamples in which each pair was retained. Selection frequency is not an edge false-discovery rate.',
      'Leiden modularity clustering was applied to the full candidate graph at resolution 1 until convergence, using seeds 2001–2050. The greatest modularity at this fixed graph and resolution was retained; ties retained the first maximum. The selected seed was 2024 and modularity was '+fmt(float(metrics.modularity_primary),6)+'. All connected components were optimized jointly; partition comparisons excluded isolated candidate KOs.',
      '## Conditional enrichment tests',
      'Four thousand nine hundred and ninety-nine permutations shuffled module labels within covariate blocks defined by occurrence, log total abundance, CLR variance, degree and absolute network strength. Blocks used robust scaling and k-means with 50 starts, a target size of 35 and a minimum size of 12. All six modules × seven categories were tested as one family, including zero counts. Support required both Benjamini–Yekutieli q ≤ 0.05 and a single-step maximum standardized permutation statistic maxT P ≤ 0.05. These are conditional post-selection tests on KOs; metagenomes are not treated as independent replicates in these tests. The smallest possible permutation P is 0.0002.',
      'Near-global functions were predominant descriptively but were not enriched after conditioning on selection covariates. For major modules versus connected nodes, the observed near-global proportion was '+fmt(near/n,3)+' versus '+fmt(selection.query("analysis=='Major modules versus connected nodes' and functional_core_class=='near_global'").null_expected_proportion.iloc[0],3)+' under the null, with BY q = 1. Candidate filtering and graph selection therefore matter when interpreting a core functional profile.',
      '## Stability of the inferred partition']
    primary=partition.query("stability==0.8 and resolution==1 and algorithm=='leiden'").iloc[0]
    sub += [f'All 50 Leiden starts on the primary graph produced six major modules, but memberships varied: median pairwise adjusted Rand index (ARI) within the largest connected component was {primary.ARI_largest_median:.3f} (5th–95th percentiles {primary.ARI_largest_q05:.3f}–{primary.ARI_largest_q95:.3f}). The best partition is a reproducible reference at this resolution, not evidence that six modules is a unique biological answer.',
      f'Across 500 edge perturbations, each consensus edge was retained independently with its selection frequency and every graph was reclustered with the same 50-start Leiden policy. ARI across the {n} reference major-module KOs had median {perturb.ARI_vs_primary_modules.median():.3f} and 5th–95th percentiles {perturb.ARI_vs_primary_modules.quantile(.05):.3f}–{perturb.ARI_vs_primary_modules.quantile(.95):.3f}. Co-assignment measures whether pairs remain together and does not penalize mergers; the whole-partition ARI is reported alongside it.',
      mdtable(pd.DataFrame({'Module':stability.module_label,'KOs':stability.n_KOs,'Mean co-assignment':stability.mean_stability.map(fmt)}))]
    unique=threshold.query('mean_abs_min==0.95 and sign_consistency_min==1').sort_values('stability_min')
    sub += ['Thirty threshold combinations reduced to five distinct edge sets. Identical edge sets reused the same partition, so stochastic variation was not attributed to a threshold change. The mean-strength and sign cutoffs did not distinguish graphs within the tested grid. The following values refer to the best of 50 Leiden starts at resolution 1.',
      mdtable(pd.DataFrame({'Frequency':unique.stability_min.map(lambda x:fmt(x,2)),'Edges':unique.n_edges,'Connected KOs':unique.n_connected_nodes,'Major modules':unique.n_large_modules,'ARI':unique.ARI_vs_primary_major_modules.map(fmt)})),
      'The alternative algorithms were applied to the same full candidate graph. Leiden and Louvain used 50 starts; Infomap used a fixed seed and Walktrap was deterministic. ARI was evaluated across reference major-module nodes. Scores are agreement measures, not evidence of biological truth.',
      mdtable(pd.DataFrame({'Algorithm':algorithms.method,'Major modules':algorithms.n_large_modules,'ARI':algorithms.ARI_vs_primary_major_modules.map(fmt)})),
      '## Convergence and parameter scale',
      'The primary 500-subsample run reproduced the audited consensus network. Edge Jaccard relative to the 500-subsample graph was 0.910 at 100 subsamples, 0.942 at 200, 0.972 at 300 and 0.992 at 400. These values assess Monte Carlo convergence conditional on this cohort. Independent batches of 100 draws had median edge Jaccard 0.881; the two 250-draw halves had Jaccard 0.916. They do not establish reproducibility in an external cohort.',
      'In the prior edge-budget audit, 1,000, 2,000 and 4,000 retained associations per subsample produced 448, 1,060 and 2,435 consensus edges using 100 subsamples. Increasing density changed the connected set and partition stability. Resolution 1 and the 2,000-edge budget are reference choices supported by this audit, not a demonstrated global optimum. Enrichment was not used to optimize either parameter. Pseudocount and candidate-universe alternatives were not evaluated in this revision.',
      '## Sensitivity to metadata study units',
      f'Each of 43 BioProject units was removed in turn; the unresolved project label was evaluated separately. Every omission used 500 newly drawn locality-capped subsamples and the same 50-start Leiden policy. The CLR reference and 1,500 candidate KOs were fixed, so these tests evaluate sensitivity to sample removal conditional on candidate selection. Different BioProjects may belong to one publication; this metadata grouping does not establish 43 independent studies.',
      f'Excluding the six PRJNA437562 profiles produced the greatest change: 565 edges among 200 connected KOs, four major modules, edge Jaccard 0.111 and partition ARI 0.099. Those profiles came from two recorded Chinese localities, Hainan and Guanxi, and all four realms remained represented. This dependence does not identify its biological or technical cause and prevents a general claim of study-independent modular organization. The canonical cohort was retained; influence alone is insufficient grounds for post hoc exclusion, and any eligibility revision would require a separate protocol audit. Across BioProject omissions, edge Jaccard ranged from {known.edge_Jaccard.min():.3f} to {known.edge_Jaccard.max():.3f}, with median {known.edge_Jaccard.median():.3f}. ARI across the {n} reference major-module nodes ranged from {known.ARI_reference_major.min():.3f} to {known.ARI_reference_major.max():.3f}, with median {known.ARI_reference_major.median():.3f}. Networks contained {known.n_edges.min()}–{known.n_edges.max()} edges and {known.n_major.min()}–{known.n_major.max()} major modules. The table shows the six omissions with lowest partition agreement; the complete 44-case results and accession lists accompany the archive.']
    low=known.nsmallest(6,'ARI_reference_major')
    sub += [mdtable(pd.DataFrame({'BioProject':low.study_id,'Removed':low.n_removed,'Subsample n':low.n_subsample,'Edge Jaccard':low.edge_Jaccard.map(fmt),'ARI':low.ARI_reference_major.map(fmt)}))]
    sc=study_scores[study_scores.case.isin(known.case)].groupby('module_label').agg(jac_min=('best_overlap_Jaccard','min'),jac_med=('best_overlap_Jaccard','median'),co_min=('within_module_coassignment','min'),ret_min=('connected_KO_fraction','min')).reset_index()
    sub += ['Some BioProject omissions also change geographic coverage; removing PRJNA1215819 eliminates the Tropical Eastern Pacific from the retained cohort. Such effects cannot be separated from loss of that realm in this sensitivity analysis. Best-overlap Jaccard compares each reference module with its closest reconstructed module and penalizes both losses and mergers. The minimum summarizes the most influential BioProject omission for each metric, which need not be the same omission across columns.',
      mdtable(pd.DataFrame({'Module':sc.module_label,'Minimum Jaccard':sc.jac_min.map(fmt),'Median Jaccard':sc.jac_med.map(fmt),'Minimum co-assignment':sc.co_min.map(fmt),'Minimum KO retention':sc.ret_min.map(fmt)}))]
    unresolved=study.query("unit_type!='BioProject'").iloc[0]
    sub += [f'Removing the unresolved project-label unit ({unresolved.n_removed} profile) gave edge Jaccard {unresolved.edge_Jaccard:.3f} and partition ARI {unresolved.ARI_reference_major:.3f}. This case is not described as exclusion of a known study.',
      '## Complete primary category tests',
      'O is the observed count; E is the permutation mean. P is the upper-tail empirical probability. BY q and maxT P refer to the complete 42-test family. Categories with no enrichment remain visible to avoid reporting only significant comparisons. NA denotes an undefined standardized maxT statistic: for M05 in Signal and stress and metals, the observed and all permuted counts were zero. BY includes all 42 comparisons; maxT uses the 41 comparisons with positive permutation variance.']
    for cat,t in enrich.groupby('category',sort=True):
      t=t.sort_values('module_label')
      sub += ['### '+cat.replace(' / ',' and '),mdtable(pd.DataFrame({'Module':t.module_label,'O':t.observed_n,'E':t.expected_n.map(lambda x:fmt(x,2)),'Fold':t.fold_enrichment.map(lambda x:fmt(x,2)),'P':t.p_empirical.map(lambda x:fmt(x,4)),'BY q':t.q_BY.map(lambda x:fmt(x,4)),'maxT P':t.p_maxT_FWER.map(lambda x:fmt(x,4))}))]
    sub += ['## Data and reproducibility',
      'The archive contains original inputs, the frozen accession manifest and annotations, scripts, source tables, network objects, the assembled Figure 5 and individual panels, and both Word documents. Input checksums and error logs accompany the analysis. The numerical sensitivity engine reproduced all 1,000,000 primary edge selections across 500 subsamples before exclusions were evaluated. R 4.3.3, igraph 2.2.1 and data.table 1.17.8 define the inferential environment. The code repository contains the current scripts, text results and numerical summaries. Images, binary reports, raw matrices and fitted objects remain in the companion archive.']
    # Follow the accepted manuscript order and continuous supplementary numbering.
    first = sub.index('## Complete primary category tests')
    last = sub.index('## Data and reproducibility')
    category_tests = sub[first:last]
    del sub[first:last]
    insertion = sub.index('## Stability of the inferred partition')
    sub[insertion:insertion] = category_tests
    table_titles = [
        'Composition and persistence classes of the six major modules',
        'Supported functional categories and enrichment estimates',
        'Complete tests for amino acid and cofactor biosynthesis',
        'Complete tests for anaerobic redox and energy',
        'Complete tests for carbon and carbohydrate metabolism',
        'Complete tests for information processing',
        'Complete tests for signal and stress and metals',
        'Complete tests for transport and membrane systems',
        'Complete tests for uncharacterized and other functions',
        'Module membership consistency under edge perturbation',
        'Sensitivity to the consensus selection frequency',
        'Sensitivity to the community detection algorithm',
        'BioProject omissions with the lowest partition agreement',
        'Module sensitivity across BioProject omissions',
    ]
    numbered = []
    counter = 0
    for block in sub:
        if block.startswith('| '):
            numbered.append('**Supplementary Table S' + str(14 + counter) + '. ' + table_titles[counter] + '**')
            counter += 1
        numbered.append(block)
    assert counter == len(table_titles)
    (docs/'FUNCTIONAL_NETWORK_SUPPLEMENTARY_RESULTS.md').write_text('\n\n'.join(numbered)+'\n')
    summary={'n_samples':191,'n_candidates':1500,'n_subsamples':500,'n_connected':int(metrics.n_connected_KOs),'n_edges':len(edges),'n_major_modules':6,'n_major_KOs':n,'n_displayed_edges':displayed,'near_global_major':near,'supported_categories':supported.to_dict('records'),'module_stability':stability.to_dict('records'),'n_study_units':len(study),'study_ARIs_range':[float(known.ARI_reference_major.min()),float(known.ARI_reference_major.max())],'manuscript_words':len(main_text.split()),'caption_words':len(caption.split())}
    out=a.repository/'results/functional_network';out.mkdir(parents=True,exist_ok=True)
    # Refresh compact tables from the same verified runs used for the manuscript.
    # This prevents a stale module summary from surviving a revised partition.
    for table in out.glob('*.tsv'):
        stage = table.name.split('_', 1)[0]
        if stage in runs:
            source = runs[stage] / 'tables' / table.name
            if source.is_file() and source.resolve() != table.resolve():
                table.write_text(portable_provenance(source.read_text()))
    (out/'network_results_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
