# Functional atlas and taxonomic–functional coupling

## Analytical units and counts

The atlas contains 191 metagenomes with both functional and taxonomic profiles,
17,482 nonzero KOs and 17,906 nonzero taxonomic features. All taxonomic profiles
contain 18,145 rarefied reads. The functional profiles quantify representation
among assembled predicted proteins. A protein assigned to multiple KOs
contributes its full count to each KO; percentages therefore use the sum of KO
annotation counts as their denominator.

The accession manifest is `config/samples/functional_accessions.txt`. Exact input
filenames, selection and checksums are in [DATA_AVAILABILITY.md](../DATA_AVAILABILITY.md).
The same paired profiles are used for taxonomic–functional comparisons.
Geographic analyses use 184 profiles with valid coordinates, representing 57
coordinate-defined sites. Named localities for core pooling are a distinct unit:
26 realm–locality combinations (CIP, 10; TA, 7; WIP, 8; TEP, 1).

## Functional persistence classes

Classes are mutually exclusive and assigned in this order:

1. Absolute core: detection in all 191 metagenomes.
2. Near-global: detection in at least 95% of metagenomes, excluding the absolute core.
3. Realm-core: remaining KOs detected in every locality of at least one realm
   represented by at least three localities, with mean relative abundance at
   least 0.0001 (0.01%) in that realm.
4. Periphery: all remaining KOs.

Counts are pooled within realm–locality units before conversion to proportions.
Localities receive equal weight when calculating the realm-specific mean.
No taxonomic read-count or prevalence filter is transferred to the KO matrix.

## PCA and regional differentiation

PCA uses CLR(count + 1), centers each feature, and does not scale feature
variances. The inset summarizes summed absolute PC1 loadings for the 100 KOs
with the largest absolute loadings, grouped by persistence class. These weights
do not measure the proportion of a pairwise realm contrast explained by a class.

Primary regional inference uses one mean CLR vector per coordinate-defined site.
Realms require at least five sites: Tropical Atlantic (TA), 23; Central
Indo-Pacific (CIP), 22; Western Indo-Pacific (WIP), 11. The single Tropical
Eastern Pacific site is retained descriptively. The omnibus test uses 56 sites;
the TA–CIP contrast uses 45 sites. All three pairwise contrasts are tested and
their P values Holm-adjusted within each model family.

Aitchison PERMANOVA is implemented as the equivalent multivariate linear-model
partition of Euclidean sums of squares. Tests use 9,999 site-label permutations.
PERMDISP uses distances to group spatial medians, a sqrt[n/(n−1)] small-sample
correction and residual permutations. A nonsignificant dispersion test indicates
no detected difference; it does not establish identical dispersion.

The adjusted model includes the mean within-site log10 total KO annotation count
and uses Freedman–Lane permutations of reduced-model residuals. Unique R² is the
realm sum of squares divided by original total variation; partial R² uses realm
plus residual variation and is reported separately. Hellinger sensitivity uses
square-root relative counts, averaged within sites. Equal-site resampling uses
1,000 draws and reports effect-size distributions rather than repeated independent
significance tests.

For TA–CIP, unadjusted pseudo-F(1,43) = 16.2965, R² = 0.27483 and Holm P = 0.0003;
PERMDISP P = 0.7932. After accounting for annotation totals, pseudo-F(1,42) =
2.3421, unique R² = 0.02503, partial R² = 0.05282 and Holm P = 0.0044. PC1 closely
tracks annotation totals (Spearman ρ = 0.9960), limiting interpretation as a
purely geographic gradient. Each study occurs within one realm. A realm effect
conditional on study is not estimable, and within-study permutations cannot
exchange realm labels. One site contains multiple study identifiers, all retained
in the site metadata.

## Taxonomic–functional coupling and geographic distance decay

Pearson Mantel tests compare complete taxonomic and functional Aitchison distance
matrices for 191 matched profiles and for 57 mean-CLR site profiles. Upper-tail
and two-sided matrix-label permutation probabilities are both saved. Each test
uses 9,999 permutations.

Geographic decay uses Bray–Curtis similarity calculated separately from the
taxonomic and functional relative-abundance profiles, regressed on log10
Haversine distance in kilometres. Pairs separated by at most 0.1 km are omitted.
The 184 georeferenced profiles contribute 15,824 eligible pairs. Site analyses
average relative compositions within sites with equal sample weight, giving
1,583 eligible site pairs. The slope contrast is taxonomic minus functional.
Taxonomic and functional labels are permuted together against the fixed
geographic matrix; individual pairwise comparisons are not treated as independent
permutation units. Slope tests are two-sided and use 9,999 permutations.

Sampling sensitivity selects one paired metagenome per site in each of 1,000
draws. The resulting percentile ranges describe sensitivity to sample selection
within the observed sites; they are not population confidence intervals.
Functional decay is slower than taxonomic decay at both profile and site levels
and in all 1,000 draws. This supports partial retention of encoded functional
potential; it does not directly measure gene expression or functional redundancy.

## Reproduction and outputs

Install `config/functional_requirements.txt` and use
`workflow/04_RUN_FUNCTIONAL_ANALYSES.sh`. The default invocation previews the
commands; `--execute` runs them. Paths can be set in `config/functional.env.example`.
The scripts use relative repository defaults and write current filenames in
place. The seed is 20261013, with deterministic offsets recorded for each test.

Script 08A computes all atlas statistics and writes input checksums, software
versions, classifications, PCA scores/loadings, regional tests, pair matrices,
sample/site metadata and numerical summaries. Script 08B draws Figure 4 and
Supplementary Figure S4 directly from those outputs. Compact summaries and the
figure-generating scripts are committed; PNG images are omitted from this update.
The workflow generates the figures locally under `figures/functional/`. Full
matrices and large outputs remain in the companion archive. The site-level tests correspond to Supplementary Table S12;
coupling and geographic decay correspond to Supplementary Table S13.

The implementation reproduces the supplied source matrices exactly. Before
updating sample selection, the source classification counts and Mantel statistic
were reproduced. Repository tests independently check Euclidean sums of squares,
conditional model degrees of freedom, prevalence-class priority, and the
0.0001-fraction realm-core threshold.

Implementation references:

- [PERMANOVA](https://vegandevs.github.io/vegan/reference/adonis.html)
- [Multivariate dispersion](https://vegandevs.github.io/vegan/reference/betadisper.html)
- [Residual permutation implementation](https://rdrr.io/cran/vegan/src/R/permutest.betadisper.R)
