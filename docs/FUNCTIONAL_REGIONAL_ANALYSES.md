# Figure 5: functional prevalence and regional differentiation

## Scope and inputs

This workflow uses the same original KO matrix, corrected functional metadata and canonical 191-accession manifest as the current functional atlas. Input hashes are frozen in `config/functional_regional_input_checksums.json`. It additionally requires the atlas persistence classification for independent reconciliation and a frozen KO description dictionary for exploratory interpretation. No taxonomic object or coordinate-based complete-case subset enters this analysis. Counts represent protein-cluster KO annotation assignments, not mapped reads, expression or metabolic flux.

Prevalence and persistence use all 191 profiles. A KO is present when its count is positive. Absolute-core KOs occur in every profile; near-global KOs occur in at least 95% after excluding absolute-core KOs. Remaining KOs are realm-core when present in every named locality of at least one realm represented by three or more localities and their equal-locality mean relative annotation abundance is at least 0.0001. Counts are pooled within locality before normalization, and localities receive equal weight. Other KOs are peripheral. The resulting labels must match the supplied current atlas KO by KO. Prevalent, Shell and Cloud are separate prevalence classes, defined as ≥95%, 15% to <95% and <15% of profiles.

## Balanced comparisons

The corrected metadata support three comparison realms: Central Indo-Pacific (10 named localities), Tropical Atlantic (7) and Western Indo-Pacific (8). Tropical Eastern Pacific has one locality and remains descriptive. The primary design selects six localities without replacement per comparison realm, then one metagenome at random per locality. The alternative selects four per realm. Each design uses 499 selections, with identical selections across prevalence groups. Primary persistence-intersection comparisons use the same selections. Groups comprise Prevalent, realm-core Shell, peripheral Shell and peripheral Cloud. Absolute-core and near-global functions are pooled in the Prevalent group.

For each KO group separately, add one to counts, take natural logarithms and centre by the within-profile mean log count. Euclidean distances of these CLR subcompositions define Aitchison distances. Distances are computed once per feature set and exactly subset by selected accessions. Recentring within each feature group means effects describe that subcomposition, not additive shares of all-KO variance.

Unadjusted PERMANOVA fits realm. A complementary model fits log₁₀ total KO annotation counts and realm, reporting the marginal realm sum of squares divided by total sum of squares (unique R²). Partial R², dividing by the residual sum of squares after the count covariate, is archived separately. Unique R² is not a degrees-of-freedom correction. The log-total covariate may reflect both technical and biological differences; this adjustment is a sensitivity analysis, not a causal correction. Each main diagnostic uses 499 permutations; `adonis2(..., by="margin")` supplies the conditional realm test. Permutations are unrestricted among selected locality representatives. Because study and realm are confounded, their P-values are nominal descriptive diagnostics, not evidence that study effects have been controlled.

PERMDISP uses spatial medians with the small-sample correction (`type="median", bias.adjust=TRUE`) and 499 permutations. These settings follow the [vegan betadisper documentation](https://vegandevs.github.io/vegan/reference/betadisper.html), which warns about anti-conservative centroid tests. Regional sums of squares were independently reconciled with [vegan adonis2](https://vegandevs.github.io/vegan/reference/adonis.html); the analysis aborts on diagnostic errors. Repeated selections overlap. Percentiles are sampling-sensitivity ranges, not confidence intervals; nominal P ≤ 0.05 frequencies are not combined tests or independent replication. R² values are compared within each design; their different sample sizes preclude interpreting a larger alternative-design R² as stronger biological differentiation.

## Focused sensitivities

The paired contrast is Shell R² minus max(Prevalent R², Cloud R²), calculated within each selection. Primary and alternative unadjusted and annotation-total-adjusted contrasts are retained. Additional unadjusted checks use 99 selections each: equal-size random feature sets of 1,128 KOs per group; removal of KOs detected in fewer than ten profiles; and alternative prevalent/cloud boundaries of 90%/10% and 97.5%/20%. The occurrence filter changes the retained feature set without redefining the full-cohort classes. Equal-KO sampling redraws both profiles and KOs.

Country and study/BioProject deletions retain the primary three-realm × six-locality design. Each estimable deletion uses 99 selections; its median paired contrast is the summary unit. The supplementary sensitivity figure uses the full range across deletion medians, rather than treating pooled selections as independent exclusions. Six of nine country deletions and 42 of 43 active study/BioProject deletions were estimable. The Tropical Eastern Pacific-only study is outside the regional comparison. Ineligible deletions are retained in explicit skip tables. Class assignments remain those of the full 191-profile cohort; these checks assess sample-selection dependence, not reclassification stability.

## Functional interpretation

The exploratory KO audit decomposes the unique realm sum of squares conditional on log₁₀ total KO counts within the Shell CLR subcomposition. For each primary selection, project the response onto the realm design residualized against the intercept and count covariate. Squared projected coordinates provide per-KO contributions, which are divided by their sum. Rank KOs by mean contribution share across 499 selections. The sum of contributions independently reconciles with the R adjusted realm R² to within 2 × 10⁻¹⁴. The frozen description dictionary covers 17,359 of all 17,482 KOs and 9,815 of 9,886 Shell KOs. Missing descriptions remain explicit; descriptions are not a complete validated pathway mapping. No pathway-enrichment or activity claim is made.

## Selected process scores and locality map

The original six process groups are frozen in `config/functional_process_catalog.tsv`; its provenance records the recovered script 158 hash. The catalog is independent of the new effect sizes. Reconcile every KO with the current classification and retain only Shell KOs. K00441 is Prevalent at 182/191 and is explicitly excluded; 33 KOs remain in six scores. Script 09E averages each group's selected coordinates in the full 9,886-KO Shell CLR (pseudocount one), then z-standardizes the score across all 191 profiles with sample SD. Recentring within a process must not be used, because the mean would then be zero. Reuse the exact 499 draws for each primary and alternative design. Compute per-score univariate between-realm SS / total SS, plus unique realm SS conditional on log10 total KO counts / total SS. These are effect sizes; no process-level significance tests or pathway-enrichment tests are claimed. Heatmap means average primary-draw regional means. Process R2 is not a fraction of total KO compositional variance.

For the map, use all 191 profiles, including Tropical Eastern Pacific, and the same full Shell CLR. Average CLR coordinates within REALM × Country × core_locality. Country does not subdivide any canonical realm-locality unit in this cohort (checked). Give all 26 locality profiles equal weight in the global centroid. Map colour is the RMS Aitchison distance from each locality centroid to that global centroid. Place the 23 located localities at component-wise medians of complete observed coordinate pairs, representing 184 profiles. Guanxi, Hainan and Rabigh Lagoon have no coordinates; retain their centroids but do not place them. The map therefore contains 23 named-locality units, not the 57 coordinate-defined sites used in Figure 4. Resample profiles within each locality 500 times, recomputing the global centroid, to report uncertainty; never replace the direct map estimate with the bootstrap median. Singletons can vary only through the resampled global reference. This bootstrap does not account for study or locality selection. The 5th–95th percentiles are descriptive resampling intervals and can exclude the direct estimate.

Script 09F independently recomputes all 1,146 standardized scores, 11,976 raw/conditional score R2 values and 26 map distances in base R using linear-model residual sums of squares and independently reconstructed CLR profiles. Agreement is required within 1e-12. Objects are saved as `181_182_process_map_objects.npz`; full per-draw tables, catalog inclusion/exclusion and sample-locality manifests are exported. Scripts 158 and 179 supply the historical methods, not historical sample identifiers or effect sizes. Natural Earth basemap geometry is frozen under `config/maps` with source and SHA-256; it does not determine realm assignments.

The main Figure 5 preserves the original image's scientific structure: persistence/prevalence crosswalk, intersection effects, selected-process heatmap/effects, and geographic Shell map. Prevalence distributions, paired sensitivities, annotation-total/dispersion diagnostics and locality-bootstrap intervals appear in the supplement. The older exploratory per-KO SS decomposition remains in reproducibility tables but is not used as a substitute for process-score analysis or emphasized in the manuscript.

## Reproduction and outputs

Install the validated R packages listed in `config/functional_regional_R_packages.tsv` and Python packages in `config/functional_requirements.txt`. Place `functional_KO_counts_original.tsv.gz`, `metadata_functional.csv`, `KO_persistence_classes.tsv` and `KO_function_dictionary.tsv` together. Preview then execute:

```bash
bash workflow/06_RUN_FUNCTIONAL_REGIONAL.sh --inputs /path/to/inputs --outdir /path/to/results
bash workflow/06_RUN_FUNCTIONAL_REGIONAL.sh --inputs /path/to/inputs --outdir /path/to/results --execute
```

When valid balanced outputs already exist, add `--complete-panels` to recompute the process/map analyses and figures while verifying and reusing those outputs.

Set `RSCRIPT` and `PYTHON` if executables have non-default paths. The workflow verifies inputs and selections, computes all retained analyses, records session information and renders the assembled Figure 5, individual panels and three supplementary figures as 600-dpi PNG/TIFF and vector PDF. Filenames and output directories are canonical and replaced on rerun. It does not invoke Git or publish files. Full per-selection statistics, selected accessions, classifications, skipped exclusions and fitted objects accompany the reproducibility package. Compact verified summaries are in `results/functional_regional/`.

Country boundaries use black lines of 0.65 pt in the global map and both insets. To redraw only panel d and the assembled Figure 5 as PNG files from existing validated tables, without rerunning analyses or regenerating the other exports:

```bash
python3 scripts/09_functional_regional_differentiation/09B_PLOT_FUNCTIONAL_REGIONAL_DIFFERENTIATION.py \
  --tables /path/to/results/tables --outdir /path/to/results/plots --map-only
```

The Word report builder uses `python-docx` from the same requirements file. It writes the Results with the assembled figure and the Supplementary Results as separate documents; `--main-only` updates only the main report. Both use plain 11-point text and bold titles.

```bash
python3 tools/build_functional_regional_reports.py \
  --results /path/to/results --outdir /path/to/reports
```

`tests/test_regional_functional_statistics.R` checks the fast R² and F calculations against `vegan`, verifies the conditional realm sum of squares, and rejects missing, negative and non-numeric counts. Run it with the validated R environment:

```bash
Rscript --vanilla tests/test_regional_functional_statistics.R
```
