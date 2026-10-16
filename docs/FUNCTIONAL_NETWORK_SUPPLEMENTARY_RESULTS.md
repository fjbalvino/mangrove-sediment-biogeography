# Supplementary results for the functional network

The network from 191 metagenomes contained 1,062 associations among 327 connected KOs. Six major modules contained 236 KOs and 999 associations. Information-processing and transport/membrane enrichment retained support across annotation treatments; carbon/carbohydrate enrichment did not. Sensitivity analyses assess distinct sources of uncertainty and do not establish a unique biological partition.

## Module composition

Major modules contained at least 10 KOs and were labelled by decreasing size, then internal edge count. Their persistence classes comprised 9 absolute-core, 194 near-global and 33 realm-core KOs. No peripheral KO occurred in these modules. The complete connected network comprised 35 components; the largest contained 250 KOs. There were 43 communities among connected nodes, including the six major modules. Module identifiers are size ranks, not conserved biological identities across parameter settings.

**Supplementary Table S14. Composition and persistence classes of the six major modules**

| Module | KOs | Absolute core | Near-global | Realm-core |
| --- | --- | --- | --- | --- |
| M01 | 57 | 7 | 47 | 3 |
| M02 | 52 | 0 | 45 | 7 |
| M03 | 49 | 0 | 39 | 10 |
| M04 | 27 | 0 | 22 | 5 |
| M05 | 27 | 1 | 26 | 0 |
| M06 | 24 | 1 | 15 | 8 |

## Functional category support and annotation sensitivity

Primary categories were assigned from the frozen KO-description dictionary using the existing keyword rules and their deterministic tie order. Adding the incomplete historical KEGG MODULE text was evaluated with exactly the same permutation labels. Neither treatment is an independent biological validation. Fold enrichment is the observed count plus 0.5 divided by the mean permuted count plus 0.5. A supported label identifies an enriched category, not the function of every KO in the module.

**Supplementary Table S15. Supported functional categories and enrichment estimates**

| Module | Category | KOs | Fold | BY q | maxT P |
| --- | --- | --- | --- | --- | --- |
| M01 | Information processing | 28 | 2.15 | 0.0182 | 0.0002 |
| M03 | Carbon/carbohydrate | 5 | 3.45 | 0.0242 | 0.0020 |
| M06 | Transport/membrane | 11 | 2.62 | 0.0182 | 0.0020 |

The carbon/carbohydrate category in M03 contained five KOs: K01897, K21071, K00700, K08483 and K01912. The combined annotation assigned seven KOs to this category, but its support declined (BY q = 0.1575; maxT P = 0.0980). Increased category counts alone therefore did not imply stronger conditional enrichment. No module passed both criteria for the anaerobic redox/energy category. M02, M04 and M05 are labelled Diverse functions because no primary category passed both criteria.

The partial pathway snapshot supplied annotations for 102 of the 236 major-module KOs. Unmapped KOs cannot be interpreted as lacking pathway membership. KEGG MODULE enrichment was not tested because a complete, versioned mapping of the analysis universe was unavailable.

## Sampling and network inference

The accession manifest retained 191 profiles, including TX1_3 and all 23 Guadeloupe profiles, and excluded air controls SRR31095011, SRR31095534 and SRR31095535. This is the functional atlas cohort, not the balanced subset used for regional differentiation in Figure 6. Counts and metadata were joined by run_accession, and atlas persistence assignments were verified independently.

Of 17,482 detected KOs, 13,308 met the association filter of occurrence in at least 10 profiles and at least 10 total assignments. A pseudocount of 1 was added before CLR transformation within this fixed reference universe. The network filter required occurrence in at least 20 profiles and at least 50 assignments; 10,614 KOs qualified. The 5,000 highest CLR-variance KOs were retained, followed by the top 1,500 network candidates, with KO identifiers resolving ties. Persistence labels and enrichment results were not used to select candidates.

Five hundred subsamples without replacement were generated with seed 123, capped at four profiles per marine realm × country × locality stratum. Each subsample contained 81 profiles from 26 strata: 35 Central Indo-Pacific, 23 Tropical Atlantic, 20 Western Indo-Pacific and three Tropical Eastern Pacific profiles. Thirty-three profiles occurred in every subsample. This design limits local replication but neither equalizes realms nor removes study effects.

For each subsample, proportionality was calculated as 2 cov(x,y)/(var(x)+var(y)) in CLR space. The 2,000 strongest absolute associations were retained. Consensus required selection frequency ≥ 0.80, complete sign consistency and mean absolute proportionality ≥ 0.95. Mean weights were calculated over subsamples in which each pair was retained. Selection frequency is not an edge false-discovery rate.

Leiden modularity clustering was applied to the full candidate graph at resolution 1 until convergence, using seeds 2001–2050. The greatest modularity at this fixed graph and resolution was retained; ties retained the first maximum. The selected seed was 2024 and modularity was 0.569437. All connected components were optimized jointly; partition comparisons excluded isolated candidate KOs.

## Conditional enrichment tests

Four thousand nine hundred and ninety-nine permutations shuffled module labels within covariate blocks defined by occurrence, log total abundance, CLR variance, degree and absolute network strength. Blocks used robust scaling and k-means with 50 starts, a target size of 35 and a minimum size of 12. All six modules × seven categories were tested as one family, including zero counts. Support required both Benjamini–Yekutieli q ≤ 0.05 and a single-step maximum standardized permutation statistic maxT P ≤ 0.05. These are conditional post-selection tests on KOs; metagenomes are not treated as independent replicates in these tests. The smallest possible permutation P is 0.0002.

Near-global functions were predominant descriptively but were not enriched after conditioning on selection covariates. For major modules versus connected nodes, the observed near-global proportion was 0.822 versus 0.810 under the null, with BY q = 1. Candidate filtering and graph selection therefore matter when interpreting a core functional profile.

## Complete primary category tests

O is the observed count; E is the permutation mean. P is the upper-tail empirical probability. BY q and maxT P refer to the complete 42-test family. Categories with no enrichment remain visible to avoid reporting only significant comparisons. NA denotes an undefined standardized maxT statistic: for M05 in Signal and stress and metals, the observed and all permuted counts were zero. BY includes all 42 comparisons; maxT uses the 41 comparisons with positive permutation variance.

### Amino acid and cofactor biosynthesis

**Supplementary Table S16. Complete tests for amino acid and cofactor biosynthesis**

| Module | O | E | Fold | P | BY q | maxT P |
| --- | --- | --- | --- | --- | --- | --- |
| M01 | 0 | 6.05 | 0.08 | 1.0000 | 1.0000 | 1.0000 |
| M02 | 8 | 5.84 | 1.34 | 0.1664 | 1.0000 | 1.0000 |
| M03 | 8 | 4.59 | 1.67 | 0.0510 | 1.0000 | 0.7856 |
| M04 | 1 | 2.15 | 0.57 | 0.9086 | 1.0000 | 1.0000 |
| M05 | 4 | 2.96 | 1.30 | 0.3362 | 1.0000 | 1.0000 |
| M06 | 2 | 1.41 | 1.31 | 0.4360 | 1.0000 | 1.0000 |

### Anaerobic redox and energy

**Supplementary Table S17. Complete tests for anaerobic redox and energy**

| Module | O | E | Fold | P | BY q | maxT P |
| --- | --- | --- | --- | --- | --- | --- |
| M01 | 0 | 6.46 | 0.07 | 1.0000 | 1.0000 | 1.0000 |
| M02 | 5 | 6.89 | 0.74 | 0.8936 | 1.0000 | 1.0000 |
| M03 | 15 | 8.42 | 1.74 | 0.0048 | 0.2181 | 0.1832 |
| M04 | 4 | 4.69 | 0.87 | 0.7434 | 1.0000 | 1.0000 |
| M05 | 4 | 2.17 | 1.69 | 0.1502 | 1.0000 | 1.0000 |
| M06 | 6 | 5.37 | 1.11 | 0.4592 | 1.0000 | 1.0000 |

### Carbon and carbohydrate metabolism

**Supplementary Table S18. Complete tests for carbon and carbohydrate metabolism**

| Module | O | E | Fold | P | BY q | maxT P |
| --- | --- | --- | --- | --- | --- | --- |
| M01 | 0 | 1.75 | 0.22 | 1.0000 | 1.0000 | 1.0000 |
| M02 | 1 | 0.87 | 1.09 | 0.6302 | 1.0000 | 1.0000 |
| M03 | 5 | 1.09 | 3.45 | 0.0004 | 0.0242 | 0.0020 |
| M04 | 0 | 0.66 | 0.43 | 1.0000 | 1.0000 | 1.0000 |
| M05 | 0 | 1.00 | 0.33 | 1.0000 | 1.0000 | 1.0000 |
| M06 | 0 | 0.62 | 0.45 | 1.0000 | 1.0000 | 1.0000 |

### Information processing

**Supplementary Table S19. Complete tests for information processing**

| Module | O | E | Fold | P | BY q | maxT P |
| --- | --- | --- | --- | --- | --- | --- |
| M01 | 28 | 12.73 | 2.15 | 0.0002 | 0.0182 | 0.0002 |
| M02 | 5 | 8.35 | 0.62 | 0.9644 | 1.0000 | 1.0000 |
| M03 | 3 | 8.94 | 0.37 | 0.9986 | 1.0000 | 1.0000 |
| M04 | 1 | 4.35 | 0.31 | 0.9944 | 1.0000 | 1.0000 |
| M05 | 6 | 5.92 | 1.01 | 0.5738 | 1.0000 | 1.0000 |
| M06 | 0 | 2.71 | 0.16 | 1.0000 | 1.0000 | 1.0000 |

### Signal and stress and metals

**Supplementary Table S20. Complete tests for signal and stress and metals**

| Module | O | E | Fold | P | BY q | maxT P |
| --- | --- | --- | --- | --- | --- | --- |
| M01 | 0 | 0.07 | 0.88 | 1.0000 | 1.0000 | 1.0000 |
| M02 | 0 | 0.19 | 0.73 | 1.0000 | 1.0000 | 1.0000 |
| M03 | 0 | 0.26 | 0.66 | 1.0000 | 1.0000 | 1.0000 |
| M04 | 1 | 0.21 | 2.10 | 0.2128 | 1.0000 | 0.8734 |
| M05 | 0 | 0.00 | 1.00 | 1.0000 | 1.0000 | NA |
| M06 | 0 | 0.27 | 0.65 | 1.0000 | 1.0000 | 1.0000 |

### Transport and membrane systems

**Supplementary Table S21. Complete tests for transport and membrane systems**

| Module | O | E | Fold | P | BY q | maxT P |
| --- | --- | --- | --- | --- | --- | --- |
| M01 | 5 | 6.92 | 0.74 | 0.8826 | 1.0000 | 1.0000 |
| M02 | 1 | 4.19 | 0.32 | 0.9930 | 1.0000 | 1.0000 |
| M03 | 3 | 6.53 | 0.50 | 0.9836 | 1.0000 | 1.0000 |
| M04 | 6 | 3.43 | 1.66 | 0.1024 | 1.0000 | 1.0000 |
| M05 | 2 | 3.05 | 0.70 | 0.8604 | 1.0000 | 1.0000 |
| M06 | 11 | 3.89 | 2.62 | 0.0002 | 0.0182 | 0.0020 |

### Uncharacterized and other

**Supplementary Table S22. Complete tests for uncharacterized and other functions**

| Module | O | E | Fold | P | BY q | maxT P |
| --- | --- | --- | --- | --- | --- | --- |
| M01 | 24 | 23.02 | 1.04 | 0.4352 | 1.0000 | 1.0000 |
| M02 | 32 | 25.67 | 1.24 | 0.0204 | 0.7414 | 0.5370 |
| M03 | 15 | 19.17 | 0.79 | 0.9524 | 1.0000 | 1.0000 |
| M04 | 14 | 11.51 | 1.21 | 0.2060 | 1.0000 | 1.0000 |
| M05 | 11 | 11.90 | 0.93 | 0.7286 | 1.0000 | 1.0000 |
| M06 | 5 | 9.73 | 0.54 | 0.9964 | 1.0000 | 1.0000 |

## Stability of the inferred partition

All 50 Leiden starts on the primary graph produced six major modules, but memberships varied: median pairwise adjusted Rand index (ARI) within the largest connected component was 0.754 (5th–95th percentiles 0.604–0.940). The best partition is a reproducible reference at this resolution, not evidence that six modules is a unique biological answer.

Across 500 edge perturbations, each consensus edge was retained independently with its selection frequency and every graph was reclustered with the same 50-start Leiden policy. ARI across the 236 reference major-module KOs had median 0.852 and 5th–95th percentiles 0.649–0.944. Co-assignment measures whether pairs remain together and does not penalize mergers; the whole-partition ARI is reported alongside it.

**Supplementary Table S23. Module membership consistency under edge perturbation**

| Module | KOs | Mean co-assignment |
| --- | --- | --- |
| M01 | 57 | 0.965 |
| M02 | 52 | 0.761 |
| M03 | 49 | 0.757 |
| M04 | 27 | 0.737 |
| M05 | 27 | 0.881 |
| M06 | 24 | 0.826 |

Thirty threshold combinations reduced to five distinct edge sets. Identical edge sets reused the same partition, so stochastic variation was not attributed to a threshold change. The mean-strength and sign cutoffs did not distinguish graphs within the tested grid. The following values refer to the best of 50 Leiden starts at resolution 1.

**Supplementary Table S24. Sensitivity to the consensus selection frequency**

| Frequency | Edges | Connected KOs | Major modules | ARI |
| --- | --- | --- | --- | --- |
| 0.70 | 1362 | 350 | 6 | 0.688 |
| 0.75 | 1217 | 337 | 6 | 0.709 |
| 0.80 | 1062 | 327 | 6 | 1.000 |
| 0.85 | 906 | 311 | 8 | 0.664 |
| 0.90 | 761 | 293 | 6 | 0.695 |

The alternative algorithms were applied to the same full candidate graph. Leiden and Louvain used 50 starts; Infomap used a fixed seed and Walktrap was deterministic. ARI was evaluated across reference major-module nodes. Scores are agreement measures, not evidence of biological truth.

**Supplementary Table S25. Sensitivity to the community detection algorithm**

| Algorithm | Major modules | ARI |
| --- | --- | --- |
| louvain | 6 | 0.928 |
| infomap | 4 | 0.491 |
| walktrap | 4 | 0.393 |
| leiden | 6 | 1.000 |

## Convergence and parameter scale

The primary 500-subsample run reproduced the audited consensus network. Edge Jaccard relative to the 500-subsample graph was 0.910 at 100 subsamples, 0.942 at 200, 0.972 at 300 and 0.992 at 400. These values assess Monte Carlo convergence conditional on this cohort. Independent batches of 100 draws had median edge Jaccard 0.881; the two 250-draw halves had Jaccard 0.916. They do not establish reproducibility in an external cohort.

In the prior edge-budget audit, 1,000, 2,000 and 4,000 retained associations per subsample produced 448, 1,060 and 2,435 consensus edges using 100 subsamples. Increasing density changed the connected set and partition stability. Resolution 1 and the 2,000-edge budget are reference choices supported by this audit, not a demonstrated global optimum. Enrichment was not used to optimize either parameter. Pseudocount and candidate-universe alternatives were not evaluated in this revision.

## Sensitivity to metadata study units

Each of 43 BioProject units was removed in turn; the unresolved project label was evaluated separately. Every omission used 500 newly drawn locality-capped subsamples and the same 50-start Leiden policy. The CLR reference and 1,500 candidate KOs were fixed, so these tests evaluate sensitivity to sample removal conditional on candidate selection. Different BioProjects may belong to one publication; this metadata grouping does not establish 43 independent studies.

Excluding the six PRJNA437562 profiles produced the greatest change: 565 edges among 200 connected KOs, four major modules, edge Jaccard 0.111 and partition ARI 0.099. Those profiles came from two recorded Chinese localities, Hainan and Guanxi, and all four realms remained represented. This dependence does not identify its biological or technical cause and prevents a general claim of study-independent modular organization. The canonical cohort was retained; influence alone is insufficient grounds for post hoc exclusion, and any eligibility revision would require a separate protocol audit. Across BioProject omissions, edge Jaccard ranged from 0.111 to 0.992, with median 0.924. ARI across the 236 reference major-module nodes ranged from 0.099 to 1.000, with median 0.843. Networks contained 565–1101 edges and 4–8 major modules. The table shows the six omissions with lowest partition agreement; the complete 44-case results and accession lists accompany the archive.

**Supplementary Table S26. BioProject omissions with the lowest partition agreement**

| BioProject | Removed | Subsample n | Edge Jaccard | ARI |
| --- | --- | --- | --- | --- |
| PRJNA437562 | 6 | 75 | 0.111 | 0.099 |
| PRJEB72453 | 2 | 79 | 0.784 | 0.542 |
| PRJNA1215819 | 3 | 78 | 0.750 | 0.545 |
| PRJNA954358 | 32 | 69 | 0.731 | 0.574 |
| PRJNA1159532 | 6 | 77 | 0.640 | 0.602 |
| PRJEB50284 | 7 | 77 | 0.875 | 0.672 |

Some BioProject omissions also change geographic coverage; removing PRJNA1215819 eliminates the Tropical Eastern Pacific from the retained cohort. Such effects cannot be separated from loss of that realm in this sensitivity analysis. Best-overlap Jaccard compares each reference module with its closest reconstructed module and penalizes both losses and mergers. The minimum summarizes the most influential BioProject omission for each metric, which need not be the same omission across columns.

**Supplementary Table S27. Module sensitivity across BioProject omissions**

| Module | Minimum Jaccard | Median Jaccard | Minimum co-assignment | Minimum KO retention |
| --- | --- | --- | --- | --- |
| M01 | 0.209 | 0.965 | 0.171 | 0.877 |
| M02 | 0.276 | 0.836 | 0.158 | 0.462 |
| M03 | 0.106 | 0.780 | 0.035 | 0.449 |
| M04 | 0.074 | 0.667 | 0.006 | 0.148 |
| M05 | 0.108 | 0.788 | 0.068 | 0.667 |
| M06 | 0.167 | 0.500 | 0.033 | 0.417 |

Removing the unresolved project-label unit (1 profile) gave edge Jaccard 0.980 and partition ARI 0.736. This case is not described as exclusion of a known study.

## Data and reproducibility

The archive contains original inputs, the frozen accession manifest and annotations, scripts, source tables, network objects, the assembled Figure 5 and individual panels, and both Word documents. Input checksums and error logs accompany the analysis. The numerical sensitivity engine reproduced all 1,000,000 primary edge selections across 500 subsamples before exclusions were evaluated. R 4.3.3, igraph 2.2.1 and data.table 1.17.8 define the inferential environment. The code repository contains the current scripts, text results and numerical summaries. Images, binary reports, raw matrices and fitted objects remain in the companion archive.
