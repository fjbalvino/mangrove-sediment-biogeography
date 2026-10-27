# Taxonomic analyses

## Input and sampling units

All repeated taxonomic analyses start from `data/phyloseq_taxonomic.rds` and
`data/metadata_taxonomic.csv`, containing the 341 accessions listed in
`config/samples/taxonomic_accessions.txt`. Counts are already rarefied to
18,145 reads per profile. They are not rarefied again when library sizes are
equal. Input feature identifiers and counts remain unchanged before each
analysis applies its explicitly recorded taxonomic filters.

There are 323 georeferenced profiles and 122 sampling sites. Site identifiers
combine coordinate pairs rounded to five decimal places, except that the eight
Mud Island profiles are grouped as one published sampling site across campaigns
and depths. The Central Indo-Pacific, Tropical Atlantic and Western Indo-Pacific
have 69, 39 and 12 sites. Temperate Northern Atlantic and Tropical Eastern Pacific
each have one site and do not enter primary realm tests requiring at least five
sites. Sites must have a consistent regional assignment.

Core analysis uses named localities within marine realms, not each coordinate
pair. Its 339 profiles form 47 realm-locality units. Two Fujian profiles have
supported realm assignments but unresolved exact localities and are omitted
from locality-based calculations. Missing environmental observations are not
imputed. The metadata retain source-material annotations and unresolved source
questions; the reference is not restricted to verified unmanipulated sediment.

## Hill diversity and dominance

Scripts 03A and 03E calculate richness (q0), exp(Shannon entropy) (q1) and inverse
Simpson concentration (q2). Profile-level plots and Dunn–Holm letters are
**descriptive**. Primary Kruskal–Wallis tests compare site medians across the
three adequately replicated realms, with Holm adjustment across the three Hill
orders. Pairwise Dunn contrasts use Holm adjustment within each order.

At 120 sites, q0 differs among realms (H(2) = 17.6741, adjusted P = 0.0004358), as does
q1 (H(2) = 10.7926, adjusted P = 0.0090668). The q2 test is not significant
(H(2) = 5.27644, adjusted P = 0.0714883). CIP–TA contrasts support q0 and q1 differences.
TNA is not the least-diverse realm in this reference and has no independent
realm-level inference from its one site.

Script 03F calculates the fraction of retained reads in the three leading
features or identical taxonomic annotation groups within each profile, then
summarizes by site. The median is 1.7443% (IQR 1.2028–4.1430%); the maximum site
median is 90.394%. Annotation-group aggregation preserves counts. Dominance is
not inferred from the diversity tests alone.

## Environmental GAM

Scripts 03B–D fit a Gaussian identity-link GAM to log1p(Shannon entropy), with
14 environmental smooths and a locality random-effect smooth. Estimation uses
REML, shrinkage selection, basis dimension k = 6 and gamma = 1.4. Complete cases
number 277; the prespecified four-SD response rule retains 274 profiles from 33
localities. The model explains 73.2487% of deviance (adjusted R² = 0.718219).

The largest centered partial-effect amplitudes over each predictor's central 95%
range are precipitation seasonality (BIO15), soil organic carbon, soil pH and
precipitation of the wettest month (BIO13). Sand content is the other environmental
smooth retained as significant after BH adjustment. The full term table,
concurvity, k-index diagnostics and residual plots are exported. Near-complete
concurvity and residual diagnostics limit separate attribution: curves describe
conditional associations, not independent environmental effects or causation.

## Composition and geographic distance

Script 04A computes Aitchison distances from CLR-transformed counts with a
pseudocount of 1 and fits NMDS to 323 profiles (stress = 0.054323). Longitude contours
are a GAM-smoothed descriptive surface restricted to the observed convex hull.

Script 04B compares Bray–Curtis similarity from relative abundances and Sørensen
similarity from taxon occurrence against log10 Haversine distance. The displayed
50,063 pairs exceed 0.1 km. This descriptive pair set includes 15 comparisons within
the documented Mud Island site; they do not provide independent site replication.
Mantel matrix-label tests use 9,999 two-sided permutations. Correlations of
dissimilarity with log10 distance are 0.319358 and 0.286979 (both P = 0.0001).

The profile-level slope contrast is 0.00050645 (Bray–Curtis minus Sørensen;
P = 0.7446). In 1,000 draws selecting one profile per site, both slopes remain
negative in every draw. Their median values are −0.112519 and −0.127053.
The slope contrast changes magnitude with balancing; its sign does not change.
Nonsignificance of the profile-level contrast does not establish equal decay rates.

## Regional PERMANOVA and variation partitioning

Script 04C uses site-centroid CLR coordinates for primary Aitchison PERMANOVA,
with complementary Bray–Curtis checks and dispersion tests. Regional levels
require at least five sites. Tests use 9,999 permutations; 199 equal-site
subsamples evaluate imbalance. Primary Aitchison R² values are 0.159272 for
realms, 0.430132 for provinces, 0.463381 for ecoregions and 0.480333 for countries
(all P = 0.0001). Dispersion also differs, so separation cannot be attributed
exclusively to shifts in multivariate location.

Ecoregion comparisons are estimable with permutations blocked by study
(108 sites, 21 studies; three studies span ecoregions; P = 0.0001). The corresponding
realm, province and country blocked designs are not estimable under the stated
replication checks. Leave-one-country-out results show coverage-sensitive realm
effect sizes. The eight Texas BioProjects are represented by one shared study ID.

Script 04D uses 99 sites with complete environmental information. Seven environmental
PCs, two positive dbMEM axes, three spatial trend terms and marine realm jointly
explain 45.0919% adjusted variation (F(16, 82) = 6.0300; P = 0.0001). Unlike the primary
realm contrast, this descriptive partition includes all five represented realms.
Unique fractions are 1.1270% environmental, 8.0014% spatial and 3.6237% realm.
Negative shared adjusted fractions remain reported, rather than being truncated.

After conditioning on study, the joint adjusted fraction is 1.9206%
(F(12, 61) = 1.22331; P = 0.001). Unique environmental and spatial tests are not significant
(P = 0.6003 and 0.1253); the unique realm effect has zero constrained rank and is
reported as not estimable. Study–geography overlap limits causal attribution.

## Core-associated taxa and realm restriction

Script 05A pools counts within each realm-locality unit. Globally eligible
features require at least 50 reads and occurrence in at least 5 units. A feature
is core-associated if it occurs in every locality of at least one realm with
at least 3 localities and has mean relative abundance at least 0.0001 (0.01%)
in that realm. This threshold uses the locality proportions of all retained
features. It defines realm-supported persistence, not an absolute global core.

The 2,773 core-associated features are 27.4093% of 10,117 eligible features and
contribute 3,815,419 of 6,094,344 eligible reads (62.6059%). Figure 3a uses this
eligible-read denominator. Script 05C gives equal weight to each locality after
normalizing its core counts. Pseudomonadota and Actinomycetota contribute
53.7523% and 24.3269% of these locality-balanced proportions.

Script 05B independently describes realm restriction among all 13,022 detected
features in the same 339 profiles. The 476 features detected in one sampled realm
account for 1,522 of 6,151,155 retained reads (0.0247433%). One passes global
eligibility. Realm restriction does not establish endemism, and the fraction
must not be described as uniformly ineligible or of low occupancy.

Script 05D normalizes core and eligible non-core counts separately. Geographic
comparisons use 43 units and 903 pairs; climatic comparisons use 41 units and 820
pairs, with Euclidean distance across standardized BIO1, BIO2, BIO3, BIO13 and
BIO15. Mantel tests and paired slope contrasts use 9,999 two-sided locality-label
permutations. Geographic Mantel tests use raw kilometres, whereas geographic
slopes use log10 distance. Both climatic analyses use Euclidean climatic distance.

Core similarity declines more steeply with geographic distance (slope difference
−0.0303245; P = 0.0007). The climatic slope contrast is unsupported
(−0.0100552; P = 0.0780). The non-core climatic Mantel test is not significant
(P = 0.0559). These spatial comparisons do not measure temporal stability.

## Retained iCAMP estimates

The separate completed iCAMP analysis is preserved, not repeated. It uses 400
metagenomes, 1,892 retained bacterial genera matched to an external GTDB R226
bac120 phylogeny, 24,692 within-realm pairs and 999 null randomizations.
The abundance-weighted contributions are 9.812573% heterogeneous selection,
3.372473% homogeneous selection, 26.726199% dispersal limitation,
12.946459% homogenizing dispersal and 47.142296% drift, weak selection or other
undominated processes. Categories not assigned to heterogeneous or homogeneous
selection sum to 86.814954%. This is not a direct measurement of a purely
stochastic fraction. Estimates apply to the original retained, tree-matched
cohort and do not inherit the current reference's metadata or sample selection.

## Outputs

Each stage writes one current directory under `results/`, with full-precision
tables, method settings and model/session information. Figures 2 and 3 are composed
from those current outputs. The companion archive holds numerical summaries, full results and input objects;
the code repository contains no result tables. Script 06A accepts
its independent input through `ICAMP_PHYLOSEQ`; the taxonomic workflow does not
invoke it. Functional analyses are outside this taxonomic rerun.

The current manuscript places iCAMP in Figure 3. The retained plot assemblers still use their historical layout; see the [manuscript guide](MANUSCRIPT_GUIDE.md) before assembling final artwork.
