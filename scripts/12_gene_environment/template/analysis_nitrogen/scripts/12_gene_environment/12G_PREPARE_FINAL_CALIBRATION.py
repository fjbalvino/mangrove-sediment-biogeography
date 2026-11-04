"""Freeze the final conditional-design experiment before generating responses."""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd

root=Path(__file__).resolve().parents[2]
out=root/'calibration_final';out.mkdir(exist_ok=True)
g=pd.read_csv(root/'results/study_coordinate_groups.tsv',sep='\t').set_index('analysis_group')
g=g.dropna(subset=['bio1','soc_mean','phh2o_mean']).copy()
x=g[['bio1','soc_mean','phh2o_mean']].to_numpy()/[1,10,10]
scale=x.std(axis=0,ddof=1);x=(x-x.mean(axis=0))/scale
for j in range(3):g[f'x{j+1}']=x[:,j]
g['recovery']=np.log10(g.total_KO_assignments)
g['recovery']=(g.recovery-g.recovery.mean())/g.recovery.std(ddof=1)
for col in ['x1','x2','x3','recovery']:
    g[col]=g[col]-g.groupby('study_id')[col].transform('mean')
u=pd.read_csv(root/'design_review/unit_design_manifest.tsv',sep='\t')
u=u[u.analysis_group.isin(g.index)].copy()
assert len(g)>0 and len(u)>0 and not u.run_accession.duplicated().any()
d=pd.read_csv(root/'results/group_detection_ge1.tsv',sep='\t').set_index('analysis_group').loc[g.index]
a=pd.read_csv(root/'results/group_relative_representation.tsv',sep='\t').set_index('analysis_group').loc[g.index]
kos=[k for k in a if a[k].notna().any()];assert len(kos)>0
alpha=d[kos].copy();mu=a[kos].copy()
for k in kos:
    positive=np.log2(a[k].where(a[k]>0))
    for lab,ix in g.groupby('study_id').groups.items():
        p=(d.loc[ix,k].sum()+.5)/(len(ix)+1)
        alpha.loc[ix,k]=np.log(p/(1-p))
        m=positive.loc[ix].median() if positive.loc[ix].notna().any() else positive.median()
        mu.loc[ix,k]=m
assert np.isfinite(mu.to_numpy()).all(), 'Undefined positive nuisance margin: all-zero/unsupported genes require design review'
assert np.isfinite(g[['x1','x2','x3','recovery']].to_numpy()).all(), 'Undefined predictor/recovery scales'
markers=[k for k in ['K02588','K02567','K00368','K15864','K00376','K03385'] if k in kos]
assert markers, 'No prespecified signal markers available'
dominant=g.groupby('study_id').x1.apply(lambda z: float((z*z).sum())).idxmax()
(out/'signal_genes.txt').write_text('\n'.join(markers)+'\n')
(out/'dominant_study.txt').write_text(dominant+'\n')
g.to_csv(out/'design.tsv',sep='\t')
u[['run_accession','analysis_group','study_id','descriptive_run_weight_within_group']].to_csv(out/'run_weights.tsv',sep='\t',index=False)
alpha.to_csv(out/'detection_intercepts.tsv',sep='\t')
mu.to_csv(out/'positive_intercepts.tsv',sep='\t')
plan={
 'status':'Frozen before final Monte Carlo responses; no real gene association tests',
 'seed':2026110307,'null_replicates_per_scenario':1000,'alternative_replicates_per_scenario':500,
 'scenarios':['independent','study_dependence','dominant_heterogeneity','spatial_50km','heavy_tails'],
 'scope':'Conditional calibration of two derived response means, not a joint hurdle likelihood or sequence/assembly simulation',
 'cohort_profiles':sum(1 for x in (root/'inputs/functional_accessions.txt').read_text().splitlines() if x.strip()),'eligible_profiles':len(u),'eligible_groups':len(g),'study_labels':g.study_id.nunique(), 'family_genes':kos,'family_tests':len(kos)*6,'dominant_study':dominant,
 'predictor_SD':dict(zip(['Temperature','Organic carbon','pH'],scale.tolist())),
 'covariates':'Three environmental variables, complete-design SD scaled, study centred; fixed study intercepts. Detection also includes study-centred standardized log10 total KO assignments as a recovery proxy, not read depth.',
 'detection':'Group fraction from run-level Bernoulli outcomes using unchanged BioSample/run averaging weights. Logit marginal probability = smoothed study-KO intercept + 0.75 recovery + environmental slopes. No invented binomial denominator.',
 'positive':'Synthetic log2 positive group-mean relative representation generated only for groups with detection fraction >0; conditional mean uses study-KO median and environmental slopes. Separate detection and positive innovations preserve known conditional truth.',
 'gene_dependence':'Gaussian copula common component correlation 0.4 across the available panel genes, independently in detection and positive components.',
 'independent':'Independent run Gaussian latent gates and independent group Gaussian positive errors.',
 'study_dependence':'Gate latent variance: 0.4 shared study + 0.3 shared group + 0.3 independent run. Positive errors: Gaussian unit noise plus independent unit-variance study-specific slopes on each predictor.',
 'dominant_heterogeneity':'Independent gates; positive Gaussian SD 4 in the study with the largest within-study temperature sum of squares and 1 elsewhere; study fixed before simulation.',
 'spatial_50km':'Gate and positive Gaussian copulas have correlation 0.5 I + 0.5 exp(-great-circle distance/50 km), allowing dependence across studies.',
 'heavy_tails':'Independent gates; positive Gaussian copula transformed through t3 quantiles / sqrt(3), variance one.',
 'truth':'All environmental slopes zero under null. Mixed alternative sets all three detection slopes to log(2), positive slopes to 0.5 for available prespecified markers, zero for other family genes.',
 'signal_genes':markers,
 'fitting':'R lm for positives; quasibinomial(logit) GLM for fractional detection; conventional residual-t reference and official clubSandwich CR2/Satterthwaite.',
 'boundary_policy':'Profile saturated one-group study intercepts out; detection also profiles all-zero/all-one study strata out. Nonconvergence, rank deficiency, residual df <=1, fitted detection mean outside (1e-8,1-1e-8), nonfinite or nonpositive uncertainty gives no inference. No substitute estimator.',
 'family':'Fixed family = available panel genes x3 predictors x2 response components, separately per uncertainty estimator. BH q=0.05; failed P values assigned 1 in the fixed family and failures reported separately.',
 'metrics':'Conditional-on-valid bias, empirical SD, SE, coverage, rejection and MC uncertainty; unconditional BH FDR and power; valid-fit fraction. Store sufficient statistics and every replicate FDR, not redundant individual fitted objects.',
 'decision_rule':'For every null gene/predictor/component/scenario cell, >=95% numerically valid fits, upper 95% Wilson bound of type-I error <=0.075, lower 95% Wilson bound of coverage >=0.925. Across scenarios and signals, upper FDR Monte Carlo bound <=0.075. These are study-specific tolerance bounds, not universal standards or simultaneous confidence guarantees. Conservatism and power reported; if no candidate passes, retain descriptive/exploratory inference, with no further candidate search in this cycle.',
 'limitations':['Conditional on observed predictors and nuisance margins; no validation of coordinate or environmental measurement error','Annotation uncertainty, compositional denominator changes and missing environmental data mechanism not regenerated','Study labels are assumed blocks, not certified independent biological replicates','No MaAsLin3 execution or claim of full likelihood calibration'],
 'input_hashes':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'results/study_coordinate_groups.tsv',root/'design_review/unit_design_manifest.tsv',root/'results/group_detection_ge1.tsv',root/'results/group_relative_representation.tsv',root/'inputs/metadata_functional.csv']},
 'references':['10.1002/sim.8086','10.1002/(SICI)1099-1255(199611)11:6<619::AID-JAE418>3.0.CO;2-1','10.1080/07350015.2016.1247004','10.1080/07350015.2023.2174123','10.1111/j.2517-6161.1995.tb02031.x']
}
(out/'final_simulation_plan.json').write_text(json.dumps(plan,indent=2)+'\n')
print(json.dumps({'groups':len(g),'profiles':len(u),'genes':len(kos),'plan_written':True}))
