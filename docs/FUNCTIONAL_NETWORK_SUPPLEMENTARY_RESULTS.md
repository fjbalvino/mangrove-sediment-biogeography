# Supplementary results for the functional network

The consensus network inferred from 191 metagenomes contained 1,060 associations among 327 connected KOs. Five major modules contained 238 KOs and 989 associations. These modules were dominated by near-global functions, but the predominance was not supported as enrichment after conditioning on KO occurrence, abundance, variance and network connectivity. Functional enrichment and module stability supported different conclusions and were evaluated separately.

## Module composition

Major modules contained at least 10 KOs and were labelled by decreasing size. Their persistence classes comprised 9 absolute-core, 197 near-global and 32 realm-core KOs. No peripheral KO occurred in these five modules. The complete connected network comprised 33 connected components; the largest contained 254 KOs.

| Module | KOs | Absolute core | Near-global | Realm-core |
| --- | --- | --- | --- | --- |
| M01 | 65 | 7 | 55 | 3 |
| M02 | 60 | 1 | 53 | 6 |
| M03 | 59 | 0 | 52 | 7 |
| M04 | 42 | 1 | 29 | 12 |
| M05 | 12 | 0 | 8 | 4 |

The information-processing module M01 contained 29 KOs assigned to that category (2.08-fold enrichment); M04 contained 14 amino acid/cofactor biosynthesis KOs (2.34-fold); and M05 contained 11 transport/membrane KOs (3.97-fold). All three pairs passed both multiplicity criteria (BY q = 0.00968; maxT P = 0.0002, 0.0026 and 0.0002, respectively). No redox/energy module–category pair passed both criteria. Broad categories were assigned using the frozen functional dictionary and the documented text-matching rules. Of the 238 KOs, 15 had tied category scores, resolved by the existing rule order, and 81 were classified as uncharacterized/other.

The supplied KEGG MODULE snapshot covered historical network nodes only. It provided annotations for 100 of the 238 current major-module KOs. Unmapped KOs cannot be interpreted as lacking pathway membership. KEGG MODULE enrichment was therefore not tested; a complete, versioned mapping of the analysis universe would be required for that inference.

## Functional category enrichment

All 35 module–category tests are reported below. O is the observed KO count and E is its mean across 4,999 blocked permutations. Fold enrichment is (O + 0.5)/(E + 0.5). P is the upper-tail permutation probability; q is the Benjamini–Yekutieli adjustment across the 35 tests; maxT P controls the same family using the maximum standardized statistic. A category is supported only when q and maxT P are both ≤ 0.05.

### Amino acid / cofactor biosynthesis

| Module | O | E | Fold | P | BY q | maxT P |
| --- | --- | --- | --- | --- | --- | --- |
| M01 | 2 | 9.61 | 0.247 | 0.9998 | 1.0000 | 1.0000 |
| M02 | 13 | 8.76 | 1.457 | 0.0538 | 0.9773 | 0.7908 |
| M03 | 5 | 9.99 | 0.524 | 0.9958 | 1.0000 | 1.0000 |
| M04 | 14 | 5.69 | 2.344 | 0.0002 | 0.0097 | 0.0026 |
| M05 | 1 | 0.95 | 1.036 | 0.6530 | 1.0000 | 1.0000 |

### Anaerobic redox / energy

| Module | O | E | Fold | P | BY q | maxT P |
| --- | --- | --- | --- | --- | --- | --- |
| M01 | 2 | 6.71 | 0.347 | 0.9964 | 1.0000 | 1.0000 |
| M02 | 7 | 7.62 | 0.924 | 0.6974 | 1.0000 | 1.0000 |
| M03 | 13 | 8.95 | 1.429 | 0.0606 | 0.9773 | 0.8186 |
| M04 | 12 | 8.06 | 1.460 | 0.0562 | 0.9773 | 0.8094 |
| M05 | 0 | 2.67 | 0.158 | 1.0000 | 1.0000 | 1.0000 |

### Carbon / carbohydrate metabolism

| Module | O | E | Fold | P | BY q | maxT P |
| --- | --- | --- | --- | --- | --- | --- |
| M01 | 4 | 4.49 | 0.902 | 0.7288 | 1.0000 | 1.0000 |
| M02 | 2 | 3.91 | 0.567 | 0.9436 | 1.0000 | 1.0000 |
| M03 | 3 | 3.11 | 0.969 | 0.6566 | 1.0000 | 1.0000 |
| M04 | 5 | 1.89 | 2.304 | 0.0226 | 0.6560 | 0.2138 |
| M05 | 0 | 0.60 | 0.454 | 1.0000 | 1.0000 | 1.0000 |

### Information processing

| Module | O | E | Fold | P | BY q | maxT P |
| --- | --- | --- | --- | --- | --- | --- |
| M01 | 29 | 13.70 | 2.077 | 0.0002 | 0.0097 | 0.0002 |
| M02 | 7 | 10.91 | 0.658 | 0.9650 | 1.0000 | 1.0000 |
| M03 | 5 | 9.07 | 0.575 | 0.9822 | 1.0000 | 1.0000 |
| M04 | 0 | 6.12 | 0.076 | 1.0000 | 1.0000 | 1.0000 |
| M05 | 0 | 1.20 | 0.294 | 1.0000 | 1.0000 | 1.0000 |

### Signal / stress / metals

| Module | O | E | Fold | P | BY q | maxT P |
| --- | --- | --- | --- | --- | --- | --- |
| M01 | 1 | 0.88 | 1.088 | 0.6724 | 1.0000 | 1.0000 |
| M02 | 1 | 0.72 | 1.225 | 0.5672 | 1.0000 | 1.0000 |
| M03 | 1 | 0.56 | 1.417 | 0.4734 | 1.0000 | 1.0000 |
| M04 | 0 | 0.59 | 0.457 | 1.0000 | 1.0000 | 1.0000 |
| M05 | 0 | 0.24 | 0.673 | 1.0000 | 1.0000 | 1.0000 |

### Transport / membrane systems

| Module | O | E | Fold | P | BY q | maxT P |
| --- | --- | --- | --- | --- | --- | --- |
| M01 | 6 | 8.09 | 0.756 | 0.8954 | 1.0000 | 1.0000 |
| M02 | 1 | 7.04 | 0.199 | 0.9994 | 1.0000 | 1.0000 |
| M03 | 10 | 6.07 | 1.599 | 0.0492 | 0.9773 | 0.6370 |
| M04 | 2 | 6.41 | 0.362 | 0.9982 | 1.0000 | 1.0000 |
| M05 | 11 | 2.40 | 3.971 | 0.0002 | 0.0097 | 0.0002 |

### Uncharacterized / other

| Module | O | E | Fold | P | BY q | maxT P |
| --- | --- | --- | --- | --- | --- | --- |
| M01 | 21 | 21.52 | 0.976 | 0.6318 | 1.0000 | 1.0000 |
| M02 | 29 | 21.04 | 1.369 | 0.0066 | 0.2395 | 0.2312 |
| M03 | 22 | 21.26 | 1.034 | 0.4672 | 1.0000 | 1.0000 |
| M04 | 9 | 13.24 | 0.691 | 0.9684 | 1.0000 | 1.0000 |
| M05 | 0 | 3.94 | 0.113 | 1.0000 | 1.0000 | 1.0000 |

## Module stability

Mean within-module co-assignment probabilities were highest for M01 and M05. Across 500 perturbations, the adjusted Rand index (ARI), calculated on the 238 KOs belonging to the primary major modules, had a median of 0.530 and 5th–95th percentiles of 0.402–0.666. These percentiles summarize sensitivity of the inferred partition; they are not population confidence intervals.

| Module | KOs | Mean co-assignment | Median node stability | Minimum node stability |
| --- | --- | --- | --- | --- |
| M01 | 65 | 0.866 | 0.918 | 0.433 |
| M02 | 60 | 0.420 | 0.418 | 0.171 |
| M03 | 59 | 0.388 | 0.378 | 0.233 |
| M04 | 42 | 0.516 | 0.562 | 0.087 |
| M05 | 12 | 0.936 | 0.944 | 0.835 |

Each edge was retained with probability equal to its consensus selection frequency. Weighted Louvain clustering was repeated on the complete candidate graph. A node’s stability is its mean co-assignment with other members of its reference module; the module mean averages those node values. Stability concerns reconstruction of the network, not temporal persistence or environmental resilience.

## Community detection sensitivity

Whole-graph alternative community detection also changed module boundaries. The additional Louvain run is a new stochastic optimization and need not reproduce the primary partition. ARI compares assignments for the 238 primary major-module KOs, although communities were inferred on the full candidate graph. Counts below refer to communities among connected nodes; major modules have at least 10 KOs.

| Algorithm | ARI | Modularity | Communities | Major modules |
| --- | --- | --- | --- | --- |
| Louvain | 0.593 | 0.544 | 41 | 6 |
| Infomap | 0.441 | 0.499 | 51 | 3 |
| Walktrap | 0.400 | 0.488 | 60 | 3 |
| Leiden | 0.522 | 0.566 | 40 | 6 |

Together, these analyses support more consistent membership in the information-processing and transport modules, with substantial uncertainty in the remaining module boundaries. They do not support a claim that the complete modular partition was invariant across reconstruction choices.

## Consensus threshold sensitivity

The 30 combinations below varied minimum edge selection frequency, mean absolute proportionality and sign consistency, but yielded only five distinct edge sets. At a fixed frequency threshold, varying the other two thresholds did not change the graph. All six combinations with frequency 0.80 reproduced the primary edge set exactly (edge Jaccard = 1), yet yielded ARI values of 0.517–0.551 because Louvain was optimized again. Variation across the full grid therefore combines changes in edge selection with stochastic community optimization and cannot be attributed wholly to threshold choice. ARI is calculated on the 238 primary major-module KOs; 91.6–100% of those KOs remained connected across the grid. Communities were inferred on the complete candidate graph.

| Frequency | Mean abs ρ | Sign | Edges | KOs | ARI | Retention |
| --- | --- | --- | --- | --- | --- | --- |
| 0.70 | 0.90 | 0.95 | 1388 | 356 | 0.493 | 100.0% |
| 0.70 | 0.90 | 1.00 | 1388 | 356 | 0.448 | 100.0% |
| 0.70 | 0.95 | 0.95 | 1388 | 356 | 0.403 | 100.0% |
| 0.70 | 0.95 | 1.00 | 1388 | 356 | 0.531 | 100.0% |
| 0.70 | 0.97 | 0.95 | 1388 | 356 | 0.504 | 100.0% |
| 0.70 | 0.97 | 1.00 | 1388 | 356 | 0.488 | 100.0% |
| 0.75 | 0.90 | 0.95 | 1233 | 341 | 0.517 | 100.0% |
| 0.75 | 0.90 | 1.00 | 1233 | 341 | 0.484 | 100.0% |
| 0.75 | 0.95 | 0.95 | 1233 | 341 | 0.500 | 100.0% |
| 0.75 | 0.95 | 1.00 | 1233 | 341 | 0.386 | 100.0% |
| 0.75 | 0.97 | 0.95 | 1233 | 341 | 0.550 | 100.0% |
| 0.75 | 0.97 | 1.00 | 1233 | 341 | 0.420 | 100.0% |
| 0.80 | 0.90 | 0.95 | 1060 | 327 | 0.537 | 100.0% |
| 0.80 | 0.90 | 1.00 | 1060 | 327 | 0.545 | 100.0% |
| 0.80 | 0.95 | 0.95 | 1060 | 327 | 0.529 | 100.0% |
| 0.80 | 0.95 | 1.00 | 1060 | 327 | 0.517 | 100.0% |
| 0.80 | 0.97 | 0.95 | 1060 | 327 | 0.551 | 100.0% |
| 0.80 | 0.97 | 1.00 | 1060 | 327 | 0.529 | 100.0% |
| 0.85 | 0.90 | 0.95 | 932 | 312 | 0.387 | 95.8% |
| 0.85 | 0.90 | 1.00 | 932 | 312 | 0.491 | 95.8% |
| 0.85 | 0.95 | 0.95 | 932 | 312 | 0.387 | 95.8% |
| 0.85 | 0.95 | 1.00 | 932 | 312 | 0.367 | 95.8% |
| 0.85 | 0.97 | 0.95 | 932 | 312 | 0.483 | 95.8% |
| 0.85 | 0.97 | 1.00 | 932 | 312 | 0.511 | 95.8% |
| 0.90 | 0.90 | 0.95 | 786 | 295 | 0.453 | 91.6% |
| 0.90 | 0.90 | 1.00 | 786 | 295 | 0.371 | 91.6% |
| 0.90 | 0.95 | 0.95 | 786 | 295 | 0.451 | 91.6% |
| 0.90 | 0.95 | 1.00 | 786 | 295 | 0.446 | 91.6% |
| 0.90 | 0.97 | 0.95 | 786 | 295 | 0.448 | 91.6% |
| 0.90 | 0.97 | 1.00 | 786 | 295 | 0.515 | 91.6% |

## Persistence class selection

Selection tests permuted class-selection labels within blocks of similar KO covariates. Comparisons against the 1,500 candidates conditioned on occurrence, total KO counts and CLR variance. The major-module versus connected-node comparison additionally conditioned on degree and absolute strength. BY q values adjust the upper-tail P values across the five persistence classes within each comparison. The two-sided P values shown are unadjusted. No near-global comparison supported enrichment.

### Connected versus 1,500 candidates

| Persistence class | Observed % | Expected % | Upper P | Two-sided P | BY q |
| --- | --- | --- | --- | --- | --- |
| Absolute core | 3.67 | 2.60 | 0.0422 | 0.0844 | 0.3996 |
| Near-global | 73.09 | 73.14 | 0.5700 | 1.0000 | 1.0000 |
| Realm-core | 20.49 | 22.41 | 0.9760 | 0.0956 | 1.0000 |
| Multi-realm periphery | 2.75 | 1.85 | 0.0700 | 0.1400 | 0.3996 |
| Realm-specific periphery | 0.00 | 0.00 | 1.0000 | 1.0000 | 1.0000 |

### Major modules versus connected nodes

| Persistence class | Observed % | Expected % | Upper P | Two-sided P | BY q |
| --- | --- | --- | --- | --- | --- |
| Absolute core | 3.78 | 4.16 | 0.8954 | 0.7044 | 1.0000 |
| Near-global | 82.77 | 81.97 | 0.2726 | 0.5452 | 1.0000 |
| Realm-core | 13.45 | 12.76 | 0.3140 | 0.6280 | 1.0000 |
| Multi-realm periphery | 0.00 | 1.12 | 1.0000 | 0.0568 | 1.0000 |
| Realm-specific periphery | 0.00 | 0.00 | 1.0000 | 1.0000 | 1.0000 |

### Major modules versus 1,500 candidates

| Persistence class | Observed % | Expected % | Upper P | Two-sided P | BY q |
| --- | --- | --- | --- | --- | --- |
| Absolute core | 3.78 | 3.00 | 0.2194 | 0.4388 | 1.0000 |
| Near-global | 82.77 | 81.86 | 0.3162 | 0.6324 | 1.0000 |
| Realm-core | 13.45 | 14.96 | 0.8980 | 0.3416 | 1.0000 |
| Multi-realm periphery | 0.00 | 0.18 | 1.0000 | 1.0000 | 1.0000 |
| Realm-specific periphery | 0.00 | 0.00 | 1.0000 | 1.0000 | 1.0000 |

Expected percentages are null means. The zero proportion of realm-specific peripheral KOs in the candidate and selected sets provides no evidence of class enrichment; its test probability is 1.0.

## Persistence class mixing

Class mixing was evaluated across the complete connected graph using 999 degree-preserving rewires. Both nominal assortativity and the fraction of edges joining two near-global KOs exceeded their null means. This is a graph-level pattern conditional on the inferred network; it does not establish direct biochemical interaction or independence from environmental and genomic influences.

| Statistic | Observed | Null mean | Upper P | Two-sided P |
| --- | --- | --- | --- | --- |
| Nominal assortativity among five functional persistence classes | 0.2452 | -0.00745 | 0.001 | 0.002 |
| Fraction of edges joining two near-global KOs | 0.7849 | 0.75725 | 0.001 | 0.002 |

## Supplementary methods

The 191-profile accession manifest selected the current cohort from the source KO count matrix. The three air-control profiles SRR31095011, SRR31095534 and SRR31095535 were excluded; TX1_3 was retained. KOs absent from every retained profile were removed, leaving 17,482 KOs. Persistence classes matched those in the functional atlas.

KOs occurring in at least 10 profiles and with at least 10 total assignments formed the association universe (13,308 KOs). A pseudocount of 1 was added before CLR transformation within this universe. Network eligibility required occurrence in at least 20 profiles and at least 50 total assignments. Of 10,614 eligible KOs, the 5,000 with greatest CLR variance were retained, followed by 1,500 network candidates. Variance ties were resolved by KO identifier. Counts represent KO annotation assignments to predicted proteins, not read-mapping abundance or measured functional activity.

One hundred subsamples without replacement retained at most four metagenomes from each marine realm × country × locality stratum. The 26 strata yielded 81 profiles per subsample. This cap limits local replication; it does not equalize realm sample sizes or remove study effects. Proportionality was twice the covariance divided by the sum of variances in CLR space. Each subsample retained the 2,000 strongest absolute associations. Consensus edges required selection frequency ≥ 0.80, complete sign consistency and mean absolute proportionality ≥ 0.95. Weighted Louvain clustering was performed on the full candidate graph; modules containing at least 10 KOs were interpreted.

Enrichment used 4,999 permutations within covariate blocks. Selection tests conditioned on occurrence, total abundance and CLR variance; module tests additionally conditioned on degree and absolute strength. Blocks were constructed by k-means on robustly scaled covariates, targeting 35 KOs per block and merging blocks smaller than 12 with the nearest block. Functional category tests permuted module labels within those blocks among the 238 major-module KOs. Supported categories required both BY q ≤ 0.05 and maxT family-wise P ≤ 0.05. These are conditional post-selection tests on KOs and graph partitions; subsamples and perturbations are not independent biological replicates.

Sensitivity analyses used 999 degree-preserving rewires, 30 threshold combinations, four community algorithms and 500 edge perturbations. The inferential environment was R 4.3.3, data.table 1.17.8 and igraph 2.2.1; the analysis seed was 123. Numerical tables were checked independently against the source inputs and intermediate outputs. Complete machine-readable statistics and code are available in the accompanying repository.
