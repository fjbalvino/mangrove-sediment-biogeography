# Taxonomic diversity, composition and community assembly

This guide describes the analyses supporting the manuscript section
“Taxonomic diversity is structured by environmental and spatial gradients”,
Figure 2, the taxonomic core results in Figure 3, and their supplements.
Sample membership is specified by fixed sample-identifier manifests. Each analysis retains
the observations required by its design, so its sample count must be reported
with its results.

## Inputs and configuration

Run the workflow from the repository root. R dependencies are listed in
`config/R_packages.tsv`. The derived data and their distribution status are
described in [DATA_AVAILABILITY.md](../DATA_AVAILABILITY.md); the accession
manifests are included in this repository, but they do not replace the count
matrices, environmental metadata or reference tree.

| Input | Default path | Configuration variable |
|---|---|---|
| Phyloseq taxonomic count object, 405 samples | `data/taxonomic_counts.rds` | `MAGENTA_PHYLOSEQ_405` |
| Sample and environmental metadata | `data/environmental_metadata.csv` | `MAGENTA_METADATA_CSV` |
| Diversity-analysis identifiers, 345 samples | `config/samples/taxonomic_diversity_accessions.txt` | `ALPHA_SAMPLE_MANIFEST` |
| Core-analysis identifiers, 401 samples | `config/samples/taxonomic_core_accessions.txt` | `CORE_SAMPLE_MANIFEST` |
| Derived core count object | `data/taxonomic_core_counts.rds` | `MAGENTA_PHYLOSEQ_401` |
| Bacterial genus reference phylogeny | `data/genus_phylogeny_GTDB_R226_bac120.newick` | `TREE_FILE` |

`BIOGEOGRAPHY_BASE` sets the project directory and `BIOGEOGRAPHY_RESULTS` sets
the results directory. The default results directory is `results/` beneath
the project root. `config/taxonomic.env.example` shows these settings.
Metadata and count objects are joined using `run_accession`. This field contains
public run accessions and project sample identifiers; every manifest entry must
match it exactly, including identifiers that are not ENA/SRA accessions.

The supplied 405-sample count matrix already contains 18,145 retained reads
per profile. The diversity scripts standardize profiles to this common depth;
sampling without replacement at the existing depth leaves these counts
unchanged. The scripts export supplied library sizes. These inputs should not
be described as unprocessed or unrarefied read counts.

Script `01H_SELECT_CORE_ANALYSIS_SAMPLES.R` selects the 401 sample identifiers in the
core manifest without changing retained count values or the feature space. It
creates the derived object used by the core, community-assembly and beta
scripts. The beta scripts subsequently restrict that object to the diversity
sample table. An existing derived object is not overwritten by 01H.

## Running the analyses

The workflow driver prints its commands by default:

```bash
bash workflow/03_RUN_TAXONOMIC_ANALYSES.sh
```

To run diversity, GAM, dominance, ordination, distance-decay, regional
PERMANOVA and variance partitioning:

```bash
bash workflow/03_RUN_TAXONOMIC_ANALYSES.sh --execute
```

Add `--core` for the taxonomic core analyses and `--assembly` for iCAMP and
its sampling-robustness analysis. Both flags can be supplied together. A
configuration file can be supplied with `--config FILE`. The driver creates
the 401-sample object if it is absent and stops when a stage fails.

| Order | Script | Principal outputs |
|---|---|---|
| 1 | `03A_CALCULATE_HILL_DIVERSITY_BY_REALM.R` | Samplewise Hill diversity, Fig. 2a and descriptive realm tests |
| 2 | `03E_HILL_DIVERSITY_SITE_INFERENCE.R` | Coordinate-site tests, q0/q2 figure and cohort accounting |
| 3 | `03F_TAXONOMIC_DOMINANCE.R` | Three-group concentration, realm-leading annotations and taxonomy audit |
| 4 | `03B_FIT_SHANNON_ENVIRONMENTAL_GAM.R` | GAM fit, model input, smooth statistics, diagnostics and effect curves |
| 5 | `03C_PLOT_TOP4_SHANNON_GAM_EFFECTS.R` | Fig. 2b and GAM summary tables |
| 6 | `03D_PLOT_REMAINING_SIGNIFICANT_GAM_EFFECTS.R` | Additional GAM effects meeting the adjusted significance criterion |
| 7 | `04A_NMDS_AITCHISON_LONGITUDE.R` | Fig. 2d, sample manifest and distance objects |
| 8 | `04B_DISTANCE_DECAY_AND_MANTEL.R` | Fig. 2c, matrix-permutation tests and site-balanced sensitivity |
| 9 | `04C_PERMANOVA_REGIONAL_STRUCTURE.R` | Site-level regional PERMANOVA, dispersion and sampling sensitivities |
| 10 | `04D_DBRDA_VARIANCE_PARTITIONING.R` | dbRDA, adjusted variation fractions and study-conditioned sensitivity |

Scripts 03A–03F are in `scripts/03_alpha_diversity/`; scripts 04A–04D are in
`scripts/04_beta_diversity/`. `_gam_outputs.R` is a shared plotting helper,
called by 03C and 03D, rather than an independent analysis stage.

## Analysis subsets and units

| Analysis | Manuscript subset | Unit used for inference or summary |
|---|---|---|
| Hill diversity distributions | 345 metagenomes | Individual profiles for descriptive plots |
| Primary realm comparisons of Hill diversity | 309 metagenomes at 128 sites in four realms | One median per coordinate site |
| GAM | 279 metagenomes from 34 named localities | Metagenomes with a locality random effect |
| Ordination, distance-decay and dominance | 327 metagenomes at 130 sites | Profile distances; additional site-balanced summaries |
| Regional PERMANOVA | Coordinate sites represented by sufficiently replicated regional groups | Site centroids; counts depend on the grouping |
| Environmental, spatial and realm partitioning | 104 sites with complete covariates | Site centroids |
| Core classification | 401 metagenomes | Counts pooled within each realm–locality unit |
| Primary iCAMP analysis | 400 metagenomes; 1,892 genera | Abundance-weighted turnover over 24,692 within-realm pairs |

A coordinate site is defined by latitude and longitude formatted to five
decimal places. A named locality is a separate metadata grouping and must not
be substituted for a coordinate site. Eighteen diversity profiles lack valid
coordinates. The Indo-Pacific and Tropical Eastern Pacific each have one
georeferenced site; they remain in descriptive plots but are excluded from
primary realm tests requiring at least five sites per group.

## Hill diversity and dominance

Hill q0 is observed richness, q1 is the exponential of Shannon entropy, and
q2 is inverse Simpson diversity. Figure 2a shows q1; the supplementary figure
shows q0 and q2. Half-violins show distributions, boxes show medians and
interquartile ranges, and points represent metagenomes. Displayed letters
summarize metagenome-level Dunn tests with Holm correction and are descriptive
because sites can contribute several profiles.

Primary comparisons use coordinate-site medians. Omnibus Kruskal–Wallis
P values are Holm-adjusted across the three diversity orders. Pairwise Dunn
P values are Holm-adjusted within each order. Script 03E exports the retained
site and sample manifests as well as full-precision statistics.

Script 03F tests concentration directly using C3: the sum of the three largest
relative abundances in each profile. Features are combined only when their
complete taxonomic annotations are identical; unresolved labels remain
explicit. The resulting entities are taxonomic annotation groups, not
necessarily identified species. Read totals are checked before and after
aggregation. C3 is summarized first by the median within each coordinate
site, then across sites, giving each site one value.

The identities of the three leading groups in a realm are calculated
separately: relative abundances are averaged within sites, and those site
means receive equal weight within the realm. Their combined realm mean is
not the same statistic as the median of profile-specific C3 values. These
dominance summaries are descriptive and use the supplied standardized matrix
without another rarefaction.

## Environmental GAM

The response is `log1p(Shannon entropy)`, where Shannon entropy is `log(q1)`.
The Gaussian identity model fits 14 environmental smooths jointly and
includes a random-effect smooth for sampling locality. Environmental terms
use thin-plate bases with `k = 6`; estimation uses REML, outer/BFGS
optimization, `select = TRUE` and `gamma = 1.4`.

The environmental inputs are human footprint, BIO1, BIO2, BIO3, BIO13, BIO15,
GBIF Shannon diversity and species richness within 2 km, soil organic carbon,
pH, sand, silt, clay and bulk density. Of the 345 diversity profiles, 281 have
complete model inputs. The response rule removes observations at least four
standard deviations from the complete-case mean on the transformed response
scale, leaving 279 profiles. Both removed sample identifiers and complete-case
decisions are exported.

The four main-figure predictors are ranked by the amplitude of their centered
partial effects across the central 95% of observed predictor values. They are
not ranked by F statistics or P values. Supplementary curves exclude these
four terms and retain environmental smooths with Benjamini–Hochberg-adjusted
P < 0.05 across the 14 environmental tests. Plotting reads the saved model and
curves without refitting; curve values and uncertainty are checked against
predictions from that model.

Full-precision output includes effective degrees of freedom, reference
degrees of freedom, F statistics, approximate P values, adjusted P values,
concurvity and basis-dimension diagnostics. Numerical convergence is reported
separately from these diagnostics. Strong overlap among the predictors and
residual patterns limit separate attribution: partial effects describe
conditional associations with regional environmental setting.

## Composition and geographic distance

Script 04A calculates Aitchison distances from centered log-ratio-transformed
counts with pseudocount 1. The NMDS shows metagenomes colored by marine realm;
longitude contours are a descriptive smooth restricted to observed ordination
space. Regional inference uses site centroids reconstructed from the same
Aitchison geometry, rather than treating plotted metagenomes as independent
sites.

Script 04B compares Bray–Curtis similarity from relative abundance with
Sørensen similarity from presence–absence. Similarity is `1 − dissimilarity`
and is plotted against log10 Haversine distance in kilometres. Pairs separated
by 0.1 km or less are excluded, leaving 51,366 pairs. Statistical tests use
9,999 sample-label matrix permutations; ordinary regression P values treating
all pairs as independent are not used. A paired metric contrast and 1,000
draws retaining one metagenome per coordinate site assess sensitivity to
unequal replication. The rates of decline should not be described as
equivalent merely because the full-data slope contrast is nonsignificant.

Script 04C reports site-level PERMANOVA alongside multivariate dispersion,
equal-site subsampling and coverage sensitivities. At least five sites per
regional level are required. Significant dispersion means regional separation
can reflect both location and within-region heterogeneity. Comparisons that
cannot be estimated under study restrictions are marked accordingly.

Script 04D partitions adjusted variation among environmental PCA, spatial
predictors and marine realm. Environmental components retain at least 90% of
predictor variation, subject to the configured component limit. Spatial terms
include geographic trends and selected distance-based Moran eigenvectors.
The 104 complete sites are also evaluated after conditioning on study. The
combined adjusted fraction is 44.3% without study conditioning and 2.8% after
it; overlap between study coverage and geography limits attribution to
environmental or biogeographic effects alone. Negative adjusted fractions are
sampling-adjusted estimates, not negative ecological contributions.

## Taxonomic core and locality balancing

Run 05A before 05B, 05C and 05D in
`scripts/05_taxonomic_core_biogeography/`. These stages classify the core,
describe realm restriction, summarize core-associated phyla, and compare
core/non-core distance-decay, respectively. They use the 401-sample core
object and named locality metadata.

Counts are pooled within each `REALM || Locality` unit. Global eligibility
requires at least 50 reads and detection in at least five such units. A taxon
is core-associated if it occurs in every sampled locality within at least one
realm represented by three or more localities and its mean relative abundance
within that realm is at least 0.01%. The abundance threshold uses each
locality's full retained read total before the global eligibility filter;
locality proportions receive equal weight within the realm.

The core read share and the phylum summary have different denominators:

- The 71.6% core read share pools reads assigned to the 3,758 core-associated
  features and divides by reads assigned to all 10,346 eligible features.
  It is not a percentage of all sequenced reads or an equal-locality average.
- The phylum composition first divides each phylum's core-associated reads by
  all core-associated reads in that locality, then averages the proportions
  across localities with equal weight. This is the meaning of
  **locality-balanced** in Figure 3b. Pseudomonadota and Actinomycetota account
  for 52.6% and 24.9% under this summary.

Core association describes consistent detection across localities within at
least one qualifying realm. It does not establish temporal persistence,
occurrence in every realm, or invariant relative abundance. A collectively
abundant core can still vary in composition across geographic and climatic
gradients.

## Phylogeny-based community assembly

Script `06A_ICAMP_GENUS_COMMUNITY_ASSEMBLY.R` combines genus-level abundance
profiles with the external GTDB R226 bac120 genus phylogeny. This bacterial
reference tree is an independent input, not a tree inferred from Kraken
assignments. The exact tree used for the manuscript analysis has SHA-256:

```text
15761713cbfbe4edcf0600bb65afbea0f39cfe6a01a7b111a92aeeba59ffd3fd
```

The primary settings are 999 null randomizations, seed 42, genera detected in
at least 10 metagenomes, no top-abundance genus cap, minimum bin size 24,
`phylo_metric = bMNTD` and `sig_index = RC`. Samples must have at least 50%
reference-tree coverage, calculated from genus-assigned counts before the
prevalence filter; this is not the fraction of all sequenced reads. The manuscript analysis contains
400 metagenomes and 1,892 genera. The main summary averages abundance-weighted
assignments over the 24,692 within-realm pairs; between-realm pairs are
excluded from that summary.

The script provides `--help` and `--preflight_only 1` for inspection of its
inputs and settings. A small `--smoke_test 1` run is a software check and does
not reproduce the primary estimates. An explicitly selected tree may be
supplied through `--tree PATH` or `TREE_FILE`; its provenance and coverage must
be retained with the results.

`06B_ICAMP_ROBUSTNESS_REALM_STUDY_LOCALITY.R` evaluates equal-realm,
equal-study, equal-locality and balanced-sampling summaries. It accepts
`--run PATH` or `ICAMP_RUN` to identify a completed primary run. Otherwise it
uses the most recent unlabelled primary-run pointer. Use an explicit run path
when several analyses coexist.

The five reported categories are null-model assignments. In particular, the
undominated category includes drift, weak selection and other undominated
processes; it is not a direct measurement of demographic drift. The 86.8%
not assigned to heterogeneous or homogeneous selection must not be relabelled
as exclusively stochastic turnover. Results apply to the retained bacterial
genera represented by the reference tree. Matching labels and documenting
tree coverage do not by themselves demonstrate a validated taxonomic synonym
crosswalk or phylogenetic niche signal.

## Outputs and reproducibility records

Each main run writes a timestamped directory. The alpha and GAM stages write
`LATEST_TAXONOMIC_ALPHA.txt` and `LATEST_SHANNON_GAM.txt` beneath the results
directory; downstream stages read these pointers. The beta stages use
`LATEST_TAXONOMIC_NMDS.txt`, `LATEST_TAXONOMIC_DISTANCE_DECAY.txt`,
`LATEST_TAXONOMIC_PERMANOVA.txt` and
`LATEST_TAXONOMIC_VARIATION_PARTITIONING.txt`.

For individual alpha/GAM stages, `ALPHA_INPUT_DIR` and `GAM_INPUT_DIR` select
existing runs explicitly. Optional output overrides are `ALPHA_OUTPUT_DIR`,
`GAM_OUTPUT_DIR`, `ALPHA_SUPPLEMENT_OUTPUT_DIR`, `DOMINANCE_OUTPUT_DIR` and
`GAM_FIGURE_OUTPUT_DIR`. Supplementary diversity and dominance outputs normally
reside beneath their alpha run; GAM figures and summary tables reside in the
GAM run's `manuscript/` directory.

Retain sample manifests, input checksums, full-precision tables, saved model
and distance objects, parameter logs and R session information with each
reported run. Rounded manuscript values should be traced to those exports.
The identifier manifests specify exact sample membership; they should not be
reconstructed from geographic labels or the results of significance tests.
