"""Freeze the extension family, DGPs and tolerances before simulation."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
r=Path(__file__).resolve().parents[1];o=r/'calibration_final';o.mkdir(exist_ok=True)
g=pd.read_csv(r/'results/study_coordinate_groups.tsv',sep='\t').set_index('analysis_group').dropna(subset=['bio1','soc_mean','phh2o_mean']).copy()
x=g[['bio1','soc_mean','phh2o_mean']].to_numpy()/[1,10,10];sd=x.std(0,ddof=1);x=(x-x.mean(0))/sd
for j in range(3):g[f'x{j+1}']=x[:,j]
g['recovery']=np.log10(g.total_KO_assignments);g.recovery=(g.recovery-g.recovery.mean())/g.recovery.std(ddof=1)
for c in ['x1','x2','x3','recovery']:g[c]-=g.groupby('study_id')[c].transform('mean')
u=pd.read_csv(r/'results/unit_design_manifest.tsv',sep='\t');u=u[u.analysis_group.isin(g.index)]
a=pd.read_csv(r/'results/group_relative_representation.tsv.gz',sep='\t',index_col=0).loc[g.index]
d=pd.read_csv(r/'results/group_detection_ge1.tsv.gz',sep='\t',index_col=0).loc[g.index]
ks=[k for k in a if a[k].notna().any()];a=a[ks];d=d[ks];alpha=d.copy();mu=a.copy()
for k in ks:
    pos=np.log2(a[k].where(a[k]>0))
    for s,ix in g.groupby('study_id').groups.items():
        p=(d.loc[ix,k].sum()+.5)/(len(ix)+1);alpha.loc[ix,k]=np.log(p/(1-p))
        mu.loc[ix,k]=pos.loc[ix].median() if pos.loc[ix].notna().any() else pos.median()
assert len(g)>0 and len(u)>0 and len(ks)>0 and np.isfinite(mu.to_numpy()).all(), 'Undefined positive nuisance margin: audit all-zero/unsupported genes before simulation'
g.to_csv(o/'design.tsv',sep='\t');u[['run_accession','analysis_group','study_id','descriptive_run_weight_within_group']].to_csv(o/'run_weights.tsv',sep='\t',index=False)
alpha.to_csv(o/'detection_intercepts.tsv',sep='\t');mu.to_csv(o/'positive_intercepts.tsv',sep='\t')
markers=['K00855','K15230','K00399','K14028','K11180','K17218','K06163','K00937','K04759','K00440','K00425','K08928']
markers=[k for k in markers if k in ks]
assert markers, 'No prespecified signal markers available'
assert np.isfinite(g[['x1','x2','x3','recovery']].to_numpy()).all(), 'Undefined predictor/recovery scales'
dominant=g.groupby('study_id').x1.apply(lambda z:float((z*z).sum())).idxmax()
plan={'created_before_simulation':True,'seed':2026110401,'null_replicates_per_scenario':1000,'mixed_replicates_per_scenario':500,'family_genes':ks,'family_tests':len(ks)*6,'family_scope':'All unique available KOs in the frozen remaining-cycle rules; overlapping KOs tested once. Nitrogen-only genes are outside this new family. No process scores or post-hoc per-cycle BH families.',
'dominant_study':dominant,'signal_genes':markers,'signal_rationale':'Prespecified representatives of carbon fixation, C1, sulfur, phosphorus, iron acquisition, hydrogen and respiratory/phototrophic context; all three predictor slopes receive the same planted signal. Selection is biochemical, not based on observed slopes.',
'scenarios':['independent','study_dependence','dominant_heterogeneity','spatial_50km','heavy_tails'],
'DGP':'Same conditional design as the frozen nitrogen procedure: smoothed study-KO detection intercepts, logit slope 0.75 recovery, run Bernoulli gates aggregated with original weights; study-KO positive medians and independent positive innovations. Separate common gene Gaussian components correlate genes 0.4. Study dependence adds 0.4 study/0.3 group/0.3 run gate variance and unit study-specific positive slopes. Dominant heterogeneity: positive SD4 in the frozen largest within-study temperature contributor. Spatial: 0.5I+0.5 exp(-great-circle/50km). Heavy tails: positive t3/sqrt(3). Mixed alternative: log(2) detection slopes, 0.5 log2-positive slopes for available prespecified markers.',
'fitting':'Fractional quasibinomial logit with study fixed effects, three predictors and recovery; OLS positive log2 group means with study fixed effects and three predictors. Conventional residual-t and CR2/Satterthwaite. NumPy implementation cross-checked against official R/clubSandwich before simulation.',
'boundary_policy':'Omit saturated singleton strata; detection also omits all-zero/all-one strata. Require >=2 studies, residual df>1, full rank, converged IRLS60 iterations (relative deviance epsilon1e-11), fitted means strictly between1e-8 and1-1e-8, finite positive SE>1e-10, finite confidence width and df>0. Invalid P=1 for fixed-family BH; failures reported.',
'decision_rule':'Every null gene/predictor/component/scenario cell: >=95% valid; upper95% Wilson type-I bound<=0.075; lower95% Wilson coverage bound>=0.925. Across scenarios/signals: upper95% Monte Carlo FDR bound<=0.075. BH q=0.05. Study-specific tolerances, not universal standards. No further candidate search if neither passes.',
'predictor_SD_ddof1':dict(zip(['Temperature','Organic carbon','pH'],sd.tolist())),
'limits':['Conditional on current environmental design and annotation margins','No annotation, coordinate, environmental measurement-error or missingness validation','Study labels are provenance blocks, not certified independent biological replicates','No full hurdle likelihood or MaAsLin3 execution','No observed-gene P or Q exported unless calibration passes'],
'input_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [o/'design.tsv',o/'run_weights.tsv',o/'detection_intercepts.tsv',o/'positive_intercepts.tsv']}}
(o/'simulation_plan.json').write_text(json.dumps(plan,indent=2)+'\n')
print('Frozen:',len(ks),'genes;',len(ks)*6,'tests;',len(markers)*6,'planted nonnulls')
