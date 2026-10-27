# Manuscript-to-code guide

Source: **A global metagenomic atlas reveals shared functional potential across geographically structured mangrove microbiomes**, manuscript and five supplementary documents reviewed on 26 October 2026. This guide records the current manuscript organization. It does not replace the manuscript or redefine its statistical units.

## Results and supplementary material

| Manuscript result | Main figure | Code and execution entry | Supplementary tables |
| --- | --- | --- | --- |
| MAGENTA: a global metagenomic atlas of mangrove-associated microbiomes | 1 | Scripts 01–02; cohort manifests | S1, S33 |
| Taxonomic diversity is structured by environmental and spatial gradients | 2 | Scripts 03–04; workflow 03 | S2–S7; figures S1–S3 |
| Non-selection categories and spatially persistent taxa | 3 | Scripts 05–06; workflow 03 with `--core`; iCAMP separately | S8–S10 |
| A widely shared functional backbone | 4 | Scripts 08; workflow 04 | S11–S13; figure S4 |
| Widely distributed KOs form co-abundance modules | 5 | Scripts 10; workflow 05 | S14–S27 |
| Regional differentiation across functional classes and selected Shell processes | 6 | Scripts 09; workflow 06 | S28–S32; figures S5–S7 |
| Analysis-specific sample definitions | All | Accession manifests and each script's audit outputs | S33 |

The taxonomic supplement covers S1–S10; the functional supplement S11–S13; the network supplement S14–S27; the regional supplement S28–S32. The fifth document is the cohort registry, Table S33.

## Cohort registry

| Analysis | Profiles and analytical units | Interpretation |
| --- | --- | --- |
| Source taxonomic object | 405 profiles; 24,045 input features | Upstream source, not the current analysis cohort |
| Taxonomic reference | 341 profiles; five realms; 24,009 detected features | 18,145 rarefied classification counts per profile |
| Primary realm diversity/composition | 312 profiles at 120 sites in three realms | One site summary; realms need at least five sites |
| Taxonomic geographic comparisons | 323 profiles at 122 coordinate-defined sites | The same profiles are shown at 43 named localities in Figure 1 |
| Shannon GAM | 274 profiles, 33 localities, 14 predictors | Response is log1p(Shannon entropy), not Hill diversity q1 |
| Taxonomic variation partitioning | 277 profiles at 99 sites | Complete covariates; study-conditioned and unconditioned results are distinct |
| Taxonomic core | 339 profiles, 47 realm–locality units | Two unresolved localities excluded; 10,117 eligible features, 2,773 core-associated |
| Core geographic/climatic comparisons | 43 / 41 locality profiles | Pool within named localities; do not substitute sampling sites |
| iCAMP | 400 profiles, six realms, 1,892 genera, 24,692 within-realm pairs | Separate completed analysis; 999 randomizations |
| Functional atlas | 191 paired profiles, 17,482 detected KOs | Select from a 194-profile, 17,508-KO source matrix; remove three air controls and 26 all-zero KOs |
| Functional geography | 184 profiles at 57 sites | Primary replicated-realm tests use 56 sites; TA–CIP uses 45 sites |
| Functional network | 191 profiles; 1,500 candidate KOs | 500 subsamples of 81 profiles; cap four per each of 26 realm–country–locality strata |
| Regional functional comparisons | 191 for classification; 188 eligible in three realms | 499 selections of 18 profiles (3 × 6 localities); alternative: 12 (3 × 4) |
| Shell locality map | All 26 locality centroids; 23 mappable localities, 184 profiles | All localities contribute to the global centroid, including the three without coordinates |

`run_accession` is the canonical profile key. Sampling site, named locality, realm and study are separate metadata fields. Keep missing coordinates and covariates explicit.

## Scientific interpretation to preserve

- Core-associated taxa account for 27.4% of eligible features and 62.6% of counts assigned to those features. These denominators differ from all detected taxa and all classified reads.
- The shared functional fraction comprises 24 absolute-core, 1,104 near-global and 1,415 realm-core KOs: 2,543 features (14.5%) and 86.2% of KO annotation counts. The remaining 14,939 KOs are peripheral.
- Prevalent/Shell/Cloud and absolute-core/near-global/realm-core/periphery are different classifications. The corrected balanced comparisons have the highest median regional differentiation among Prevalent KOs; they do not support a general Shell advantage.
- The network has 327 connected KOs and 1,062 edges; six major modules contain 236 KOs and 999 edges. Information-processing and transport/membrane enrichment persist across annotation treatments; carbon/carbohydrate enrichment is annotation-sensitive.
- Network edges are co-abundance associations. BioProject omission, especially PRJNA437562, changes module organization. They are not direct interactions or evidence of metabolic activity.
- Geography/environment overlap strongly with study identity. Assemblies, KO annotation counts and these observational designs do not establish ecological causation or expression.

## Provenance and remaining reproduction gaps

The source code keeps its established objective-based filenames and independent verification routines. Some plot assemblers retain historical panel organization: script 03G combines four taxonomic panels with an iCAMP panel in its `Figure_2` output, while the current manuscript places iCAMP in Figure 3. The regional plotting/report code retains `Figure_5*` filenames for the current manuscript's Figure 6. Use the manuscript numbering above when assembling a submission; filename compatibility does not certify the current final artwork.

The Figure 1 final atlas layout is not supplied as a dedicated validated plotting entry point. The original validated assembly-to-KO scripts are still missing from block 07. The proposed upstream 2,000-protein exclusion rule has not been independently verified. The original iCAMP reference tree and fitted pairwise objects must be recovered for a full rerun. These limitations are recorded rather than replaced with inferred analyses.

No scientific analysis was rerun as part of the repository reorganization. Numerical results in this guide come from the supplied manuscript and supplements.

When the separate data/figure deposit is published, update the manuscript's code/data availability paragraph: code, manifests and checksums belong in GitHub; numerical tables, matrices, fitted objects and figures belong in the data companion.
