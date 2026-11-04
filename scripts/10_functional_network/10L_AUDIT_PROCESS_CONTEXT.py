#!/usr/bin/env python3
"""Run the descriptive network redesign pilot and biochemical context overlay.

This is a diagnostic stage, not a replacement inference of the final network.
No observed P values, edge FDR, new modules or causal interpretations are produced.
"""
import os
for _key in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[_key] = '1'
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent / 'lib'))
from biogeochemical_context import read_matrix, export_context, ROOT


def clr(x, pseudocount):
    if pseudocount <= 0:
        raise ValueError('Pseudocount must be positive')
    z = np.log(np.asarray(x, float) + pseudocount)
    return z - z.mean(axis=0, keepdims=True)


def candidate_selection(x, cap=1500, pseudocount=1):
    reference = x.index[(x.gt(0).sum(axis=1) >= 10) & (x.sum(axis=1) >= 10)]
    y = pd.DataFrame(clr(x.loc[reference], pseudocount), index=reference, columns=x.columns)
    eligible = x.index[(x.gt(0).sum(axis=1) >= 20) & (x.sum(axis=1) >= 50)]
    eligible = sorted(set(eligible) & set(reference))
    variance = y.loc[eligible].var(axis=1, ddof=1)
    selected = variance.sort_values(ascending=False, kind='stable').head(cap).index.tolist()
    return reference, selected, y


def grouped_profiles(y, meta):
    a = y.T.join(meta[['analysis_group', 'biosample_group']], validate='one_to_one')
    first = a.groupby(['analysis_group', 'biosample_group'], sort=True).mean()
    grouped = first.groupby(level='analysis_group', sort=True).mean()
    groups = meta.groupby('analysis_group', sort=True)[['study_id', 'coordinate_group', 'lat', 'lon']].first()
    return grouped, groups.loc[grouped.index]


def prepare_covariance(x, studies=None, within=False):
    """Equal study weights, equal coordinate groups within study; descriptive."""
    x = np.asarray(x, float)
    if studies is None:
        w = np.ones(len(x)) / len(x)
        z = x - x.mean(axis=0)
        n_studies = None
    else:
        studies = np.asarray(studies)
        levels, inv, counts = np.unique(studies, return_inverse=True, return_counts=True)
        if within:
            keep = counts[inv] >= 2
            x, studies = x[keep], studies[keep]
            levels, inv, counts = np.unique(studies, return_inverse=True, return_counts=True)
            if len(levels) == 0:
                raise ValueError('No within-study contrast across coordinate groups')
        w = 1.0 / (len(levels) * counts[inv])
        if within:
            means = np.stack([x[inv == j].mean(axis=0) for j in range(len(levels))])
            z = x - means[inv]
        else:
            z = x - np.sum(x * w[:, None], axis=0)
        n_studies = len(levels)
    return z, w, len(x), n_studies


def pair_rho(x, left, right, studies=None, within=False):
    z, w, n, n_studies = prepare_covariance(x, studies, within)
    variance = np.sum(z * z * w[:, None], axis=0)
    cov = np.sum(z[:, left] * z[:, right] * w[:, None], axis=0)
    den = variance[left] + variance[right]
    out = np.full(len(left), np.nan)
    np.divide(2 * cov, den, out=out, where=den > 1e-15)
    return out, n, n_studies


def plot_pilot(out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                         'axes.spines.top':False,'axes.spines.right':False})
    figdir = out / 'figures'
    figdir.mkdir(exist_ok=True)
    def save_figure(fig, stem):
        # Publish complete files so synchronisation never observes partial PNGs.
        for suffix in ('png', 'pdf'):
            destination = figdir / f'{stem}.{suffix}'
            temporary = destination.with_name(destination.name + '.tmp')
            fig.savefig(temporary, format=suffix, dpi=240, bbox_inches='tight')
            temporary.replace(destination)
    a = pd.read_csv(out / 'process_context/cycle_selection_coverage.tsv', sep='\t')
    names = {'C':'Carbon','CH4_C1':'Methane / C1','N':'Nitrogen','S':'Sulfur',
             'P':'Phosphorus','Fe':'Iron acquisition','H2_electrons':'Hydrogen / electron transfer',
             'O2_light':'Respiration / phototrophy'}
    fig, ax = plt.subplots(figsize=(12,6.4))
    series = [('matrix_KOs','Available in matrix','#c5b7db'),
              ('candidate_KOs','Selected as network candidates','#9c179e'),
              ('connected_KOs','Connected in historical network','#ed7953'),
              ('major_module_KOs','In six major modules','#f4cc33')]
    y = np.arange(len(a)); width = .18
    for j,(col,label,color) in enumerate(series):
        bars=ax.barh(y+(j-1.5)*width,a[col],height=width,color=color,label=label)
        ax.bar_label(bars,padding=3,fontsize=8)
    ax.set_yticks(y,[names[x] for x in a.cycle]); ax.invert_yaxis()
    ax.set_xlabel('Distinct KOs within each cycle / context panel')
    ax.set_xlim(0,float(a.matrix_KOs.max())*1.20)
    ax.legend(loc='lower right',frameon=False,fontsize=9)
    ax.set_title('Biochemical coverage is lost during network selection',loc='left',fontweight='bold',pad=16)
    fig.text(.015,.02,'191 profiles | 19 of 236 major-module KOs map to the targeted catalogue. KOs shared between cycles are counted in each panel.',fontsize=9,color='#505564')
    fig.tight_layout(rect=(0,.055,1,1))
    save_figure(fig, '01_cycle_selection_coverage')
    plt.close(fig)
    evidence = pd.read_csv(out / 'process_context/reference_chain_evidence.tsv', sep='\t')
    d=evidence[(evidence.chain=='N_denitrification_Nap')&(evidence.threshold==1)&(evidence.role=='chain_step')]
    total=int(d.total.iloc[0]); joint=int(evidence[(evidence.chain=='N_denitrification_Nap')&(evidence.threshold==1)&(evidence.role=='all_selected_chain_steps')].detected.iloc[0])
    fig,ax=plt.subplots(figsize=(13,4.5)); ax.set_xlim(0,4); ax.set_ylim(0,2.3); ax.axis('off')
    genes=['napA + napB','nirK OR nirS','norB + norC','nosZ']
    chemistry=['Nitrate → nitrite','Nitrite → nitric oxide','Nitric oxide → nitrous oxide','Nitrous oxide → nitrogen']
    for j,row in enumerate(d.itertuples()):
        color=plt.cm.plasma(row.detected/total)
        ax.add_patch(FancyBboxPatch((j+.055,.9),.85,.73,boxstyle='round,pad=0.025',facecolor=color,edgecolor='#30243b',lw=.7))
        ax.text(j+.48,1.43,genes[j],ha='center',va='center',weight='bold',fontsize=12)
        ax.text(j+.48,1.13,f'{row.detected} / {total} profiles',ha='center',va='center',fontsize=11)
        ax.text(j+.48,.72,chemistry[j],ha='center',va='center',fontsize=9)
        if j<3:ax.annotate('',xy=(j+1.025,1.25),xytext=(j+.93,1.25),arrowprops={'arrowstyle':'->','lw':1.3,'color':'#77707c'})
    ax.text(.03,2.1,'Reference chain reveals evidence absent from the historical network',fontweight='bold',fontsize=14)
    ax.text(.03,1.86,f'All four selected steps co-detected in {joint}/{total} profiles; no nitrogen-catalogue KO is connected in the historical network.',fontsize=10)
    ax.text(.03,.34,'Arrows: biochemical reference only. Colour: fraction of profiles meeting each component rule (≥1 record).',fontsize=10)
    ax.text(.03,.12,'Community co-detection does not establish the same organism, enzyme activity, net nitrogen loss or a complete pathway.',fontsize=10,color='#505564')
    save_figure(fig, '02_nitrogen_chain_context')
    plt.close(fig)
    e=pd.read_csv(out / 'historical_edge_diagnostics.tsv.gz',sep='\t')
    fig,ax=plt.subplots(figsize=(7,6))
    ax.scatter(e.run_global__all__pc1,e.within_study_coordinates__all__pc1,s=10,alpha=.35,color='#8b0aa5',rasterized=True)
    ax.plot([0,1],[0,1],color='#999999',lw=1,ls='--')
    ax.set(xlim=(.94,1.005),ylim=(-.1,1.01),xlabel='Global proportionality across 191 runs',ylabel='Pooled within-study proportionality across coordinates')
    ax.set_title('Historical edges: global versus within-study association',fontsize=11)
    fig.text(.02,.015,'1,062 preselected edges; 51 coordinate-study groups in 10 study labels. Descriptive comparison, not edge significance.',fontsize=8)
    fig.tight_layout(rect=(0,.035,1,1))
    save_figure(fig, '03_association_design_sensitivity')
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['matrix', 'metadata', 'samples', 'reference-nodes', 'reference-edges',
                'reference-candidates', 'reference-design', 'outdir']:
        p.add_argument('--' + key, type=Path, required=True)
    p.add_argument('--influential-study', default='PRJNA437562')
    p.add_argument('--config', type=Path, default=ROOT / 'config/gene_biogeography')
    a = p.parse_args()
    out = a.outdir.resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError('Use an empty pilot run directory')
    out.mkdir(parents=True, exist_ok=True)
    x = read_matrix(a.matrix, a.samples)
    meta = pd.read_csv(a.metadata)
    if meta.run_accession.isna().any() or meta.run_accession.duplicated().any():
        raise ValueError('Metadata requires unique, nonmissing run_accession')
    meta = meta.set_index('run_accession').loc[x.columns].copy()
    if meta.study_id.isna().any() or not np.isfinite(meta[['lat', 'lon']].to_numpy(float)).all():
        raise ValueError('Study labels and verified coordinates required for every profile')
    if (meta.lat.abs() > 90).any() or (meta.lon.abs() > 180).any():
        raise ValueError('Invalid coordinates')
    meta['coordinate_group'] = [f'coord_{lat:.5f}_{lon:.5f}' for lat, lon in zip(meta.lat, meta.lon)]
    meta['analysis_group'] = meta.study_id + '::' + meta.coordinate_group
    meta['biosample_group'] = [str(s) if pd.notna(s) and str(s).strip() else 'run:' + r
                               for r, s in zip(meta.index, meta.sample_accession)]
    nodes = pd.read_csv(a.reference_nodes, sep='\t')
    edges = pd.read_csv(a.reference_edges, sep='\t')
    candidates = pd.read_csv(a.reference_candidates, sep='\t')
    reference, selected, y = candidate_selection(x, len(candidates))
    if selected != candidates.KO.tolist():
        raise ValueError('Reference candidates not reproduced; stop before comparing historical modules')
    np.testing.assert_allclose(y.loc[selected].var(axis=1), candidates.recalculated_CLR_variance,
                               rtol=1e-10, atol=1e-12)
    context = export_context(x, nodes, set(selected), edges, a.config, out / 'process_context')
    print('Reference candidate identities and variances reproduced; biological overlay exported', flush=True)
    sample_audit = meta[['study_id','sample_accession','coordinate_group','analysis_group','Country','Locality','REALM']].copy()
    sample_audit['KO_assignments'] = x.sum(axis=0)
    sample_audit['detected_KOs'] = x.gt(0).sum(axis=0)
    sample_audit['zero_fraction'] = x.eq(0).mean(axis=0)
    sample_audit.to_csv(out / 'sample_recovery_audit.tsv', sep='\t')
    design = meta.groupby('study_id').agg(n_runs=('study_id','size'),
        n_biosample_labels=('biosample_group','nunique'), n_coordinates=('coordinate_group','nunique'),
        n_localities=('Locality','nunique'), n_realms=('REALM','nunique'))
    design['within_coordinate_contrast_available'] = design.n_coordinates >= 2
    design.to_csv(out / 'study_design_audit.tsv', sep='\t')
    sample_groups = []
    omit = meta.study_id.eq(a.influential_study)
    if not omit.any() or omit.all():
        raise ValueError('Influential study must identify a proper subset')
    for label, ix in [(a.influential_study, omit), ('Other profiles', ~omit)]:
        d = sample_audit.loc[ix]
        sample_groups.append(dict(group=label, n=len(d), median_assignments=float(d.KO_assignments.median()),
            min_assignments=float(d.KO_assignments.min()), max_assignments=float(d.KO_assignments.max()),
            median_detected_KOs=float(d.detected_KOs.median()), median_zero_fraction=float(d.zero_fraction.median())))
    pd.DataFrame(sample_groups).to_csv(out / 'influential_study_recovery.tsv', sep='\t', index=False)
    rules = pd.read_csv(a.config / 'operational_rules.tsv', sep='\t')
    queried = set(rules.expression.str.findall(r'K\d{5}').explode())
    features = sorted(set(selected) | (queried & set(reference)))
    lookup = {k:i for i,k in enumerate(features)}
    left = np.array([lookup[k] for k in edges.KO1]); right = np.array([lookup[k] for k in edges.KO2])
    diagnostics = edges[['KO1','KO2','mean_weight','mean_abs_weight','edge_stability']].copy()
    diagnostics = diagnostics.rename(columns={'mean_weight':'historical_selected_draw_mean'})
    scenarios = []
    for pseudocount in (0.5, 1.0, 2.0):
        z = pd.DataFrame(clr(x.loc[reference], pseudocount), index=reference, columns=x.columns).loc[features]
        for retained in ('all', 'omit_influential'):
            keep = np.ones(len(meta), bool) if retained == 'all' else (~omit).to_numpy()
            zz, mm = z.loc[:, keep], meta.loc[keep]
            grouped, groups = grouped_profiles(zz, mm)
            for mode in ('run_global', 'study_coordinate_global', 'within_study_coordinates'):
                key = f'{mode}__{retained}__pc{pseudocount:g}'
                if mode == 'run_global':
                    values, n, ns = pair_rho(zz.T, left, right)
                else:
                    values, n, ns = pair_rho(grouped, left, right, groups.study_id,
                                             within=mode == 'within_study_coordinates')
                diagnostics[key] = values
                finite = np.isfinite(values)
                scenarios.append(dict(scenario=key, n_profiles=int(keep.sum()), n_units=n,
                    n_study_labels=ns, finite_edges=int(finite.sum()), median_rho=float(np.nanmedian(values)),
                    edges_rho_ge_095=int((values >= .95).sum()), negative_edges=int((values < 0).sum()),
                    interpretation='Sensitivity of historical edges; 0.95 is a diagnostic threshold, not calibrated significance'))
    # Reference-universe sensitivity is separate from a change in candidate IDs.
    broad = x.index[x.sum(axis=1) > 0]
    z = pd.DataFrame(clr(x.loc[broad], 1), index=broad, columns=x.columns).loc[features]
    diagnostics['run_global__all__pc1__all_detected_reference'], _, _ = pair_rho(z.T, left, right)
    # Unconditional summaries for the same historical edges over the exact R draws.
    design_draws = pd.read_csv(a.reference_design, sep='\t')
    if set(design_draws.run_accession) != set(x.columns):
        raise ValueError('Historical sampling design does not match the selected cohort')
    all_rho = []
    yfeat = y.loc[features]
    for _, d in design_draws.groupby('bootstrap_id', sort=True):
        all_rho.append(pair_rho(yfeat.loc[:, d.run_accession].T, left, right)[0])
    all_rho = np.array(all_rho)
    diagnostics['all_historical_draws_mean'] = np.nanmean(all_rho, axis=0)
    diagnostics['all_historical_draws_p05'] = np.nanquantile(all_rho, .05, axis=0)
    diagnostics['all_historical_draws_p95'] = np.nanquantile(all_rho, .95, axis=0)
    diagnostics['conditional_minus_unconditional_mean'] = diagnostics.historical_selected_draw_mean - diagnostics.all_historical_draws_mean
    diagnostics.to_csv(out / 'historical_edge_diagnostics.tsv.gz', sep='\t', index=False)
    pd.DataFrame(scenarios).to_csv(out / 'association_sensitivity_summary.tsv', sep='\t', index=False)
    selection = []
    for pc in (.5, 1., 2.):
        for label, data in [('all', x), ('omit_influential', x.loc[:, ~omit.to_numpy()])]:
            rr, ss, _ = candidate_selection(data, len(candidates), pc)
            selection.append(dict(cohort=label, pseudocount=pc, n_reference=len(rr), n_selected=len(ss),
                candidate_Jaccard=len(set(ss)&set(selected))/len(set(ss)|set(selected)),
                curated_candidates=len(set(ss)&queried)))
    pd.DataFrame(selection).to_csv(out / 'candidate_selection_sensitivity.tsv', sep='\t', index=False)
    # Existing gene pairs in each curated rule; no new significant network is selected.
    pairs = set()
    for ks in rules.expression.str.findall(r'K\d{5}'):
        present = sorted(set(ks) & set(features))
        pairs.update((present[i], present[j]) for i in range(len(present)) for j in range(i+1,len(present)))
    pairs = sorted(pairs)
    pi = np.array([lookup[k] for k,_ in pairs]); pj = np.array([lookup[k] for _,k in pairs])
    pr = pd.DataFrame(pairs, columns=['KO1','KO2'])
    grouped, groups = grouped_profiles(yfeat, meta)
    pr['run_global_rho'],_,_ = pair_rho(yfeat.T, pi, pj)
    pr['study_coordinate_global_rho'],_,_ = pair_rho(grouped, pi, pj, groups.study_id)
    pr['within_study_coordinate_rho'],_,_ = pair_rho(grouped, pi, pj, groups.study_id, True)
    pr.to_csv(out / 'curated_rule_pair_diagnostics.tsv.gz', sep='\t', index=False)
    report = dict(status='DESCRIPTIVE PILOT COMPLETE; final network calibration not executed',
        profiles=len(meta), matrix_KOs=len(x), detected_KOs=int((x.sum(axis=1)>0).sum()),
        reference_KOs=len(reference), historical_candidates_reproduced=True,
        study_labels=meta.study_id.nunique(), coordinate_groups=meta.coordinate_group.nunique(),
        study_coordinate_groups=meta.analysis_group.nunique(), biosample_labels=meta.biosample_group.nunique(),
        study_labels_with_multiple_coordinates=int((design.n_coordinates>=2).sum()),
        reference_edges=len(edges), reference_negative_edges=int((edges.consensus_sign<0).sum()),
        reference_draws=len(all_rho), reference_modules=int(nodes.module_label.nunique()),
        historical_draws_interpretation='Overlapping subsamples, not independent replication or confidence intervals',
        within_study_estimand='Equal-study weighted pooled covariance of study-coordinate CLR deviations from their study mean',
        weighting='Run means within BioSample label, then label means within study-coordinate; equal study weights and equal coordinate-group weights within each study',
        independent_units='Study labels, BioSample labels and coordinate groups do not by themselves certify independence',
        total_count_detected_KOs_Spearman=float(spearmanr(x.sum(axis=0),x.gt(0).sum(axis=0)).statistic),
        context=context, sample_groups=sample_groups,
        versions=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__),
        inputs={k:dict(path=str(v.resolve()),sha256=hashlib.sha256(v.read_bytes()).hexdigest())
                for k,v in vars(a).items() if isinstance(v,Path) and v.is_file()})
    (out / 'pilot_summary.json').write_text(json.dumps(report, indent=2) + '\n')
    plot_pilot(out)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
