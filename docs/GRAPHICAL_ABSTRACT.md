# Graphical abstract

## Manuscript title

A global metagenomic atlas reveals shared functional potential across geographically structured mangrove microbiomes

## Caption

Graphical abstract | Taxonomic biogeography and shared functional potential in MAGENTA. (a) The taxonomic reference comprises 341 metagenomes from 12 countries or territories and five marine realms; 323 georeferenced profiles represent 43 named localities. (b) Taxonomic composition varied with geography and environmental conditions. Across 339 profiles, realm-supported core-associated taxa comprised 27.4% of eligible taxonomic features but accounted for 62.6% of counts assigned to those features. A separate genus-level iCAMP analysis of 400 metagenomes assigned 86.8% of mean within-realm turnover to categories not attributed to selection, including drift, weak selection or other undominated processes (47.1%) and dispersal limitation (26.7%). (c) Functional profiling of 191 metagenomes detected 17,482 KOs. Absolute-core, near-global and realm-core KOs represented 14.5% of detected KOs and 86.2% of KO annotation counts. Functional similarity declined more slowly with geographic distance than taxonomic similarity. Balanced comparisons showed the highest median regional differentiation among Prevalent KOs, followed by Shell (9,886 KOs) and Cloud KOs. Six major co-abundance modules included 236 KOs; M01 and M06 were enriched in information-processing and transport/membrane functions, respectively. Module organization was sensitive to sample inclusion. Network links represent co-abundance associations, not direct biochemical interactions or metabolic activity.

## Artwork specification

Use three clearly labelled panels (a–c) on a white background, with readable
near-black text and restrained teal, blue and ochre accents. Maintain a clear
visual hierarchy when the artwork is reduced to manuscript size.

- **Panel a — Atlas coverage.** Present the 341-profile taxonomic reference and
  its 12 countries or territories and five marine realms. Identify the
  georeferenced subset separately: 323 profiles from 43 named localities. A world
  map may provide geographic context; sample markers require the actual
  coordinates and must not be invented.
- **Panel b — Taxonomic structure.** Contrast the core-associated fraction's
  27.4% of eligible features with its 62.6% of classification counts across 339
  profiles. Show the separate 400-metagenome iCAMP cohort in a distinct block.
  Label 86.8% as turnover assigned to categories not attributed to selection;
  retain the combined drift, weak-selection or undominated label for 47.1% and
  dispersal limitation for 26.7%. Illustrative microbial communities are
  conceptual and should not be assigned to specific realms.
- **Panel c — Shared functional potential.** Identify the 191-profile cohort and
  17,482 KOs. Distinguish the 14.5% share of detected KOs from the 86.2% share of
  KO annotation counts for the combined absolute-core, near-global and
  realm-core classes. Represent the slower geographic decline in functional
  similarity and the median regional-differentiation ordering
  Prevalent > Shell > Cloud as qualitative summaries, without fabricated data
  points or fitted curves. Label Shell as 9,886 KOs. Show six schematic modules
  containing 236 KOs in total, highlighting supported enrichment in M01
  (information processing) and M06 (transport/membrane functions).
  Make sample sensitivity visible and retain the interpretation of links as
  co-abundance associations.

The artwork is a conceptual summary, not a replacement for the empirical
figures. It does not depict biochemical interactions, measured activity or
experimentally established environmental selection. KO annotation counts
represent assembled protein-cluster members assigned to KOs, not mapped-read
abundances. Schematic network nodes and links do not reproduce the inferred
topology.

## Supporting results

| Component | Canonical evidence |
| --- | --- |
| Taxonomic reference and core-associated fraction | [Taxonomic methods](TAXONOMIC_ANALYSES.md) and [taxonomic summary](../results/taxonomic_summary.json) |
| Separate iCAMP cohort and process assignments | [iCAMP summary](../results/icamp/summary.json) |
| Functional persistence and annotation-count shares | [Functional class summary](../results/functional/functional_class_summary.tsv) |
| Geographic taxonomic–functional contrast | [Distance-decay estimates](../results/functional/geographic_distance_decay.json) |
| Prevalence-class ordering and KO counts | [Balanced comparison statistics](../results/functional_regional/180_verified_main_statistics.tsv) and [persistence–prevalence crosswalk](../results/functional_regional/180_persistence_prevalence_crosswalk.tsv) |
| Module counts, supported enrichment and sample sensitivity | [Figure 5 results](FUNCTIONAL_NETWORK_RESULTS.md) and [network supplementary results](FUNCTIONAL_NETWORK_SUPPLEMENTARY_RESULTS.md) |

## Figure delivery

The artwork is maintained outside the code repository. This document stores its
caption and design specification; figure files are reserved for the final data
and figure compendium. The empirical network is Figure 5 and regional functional
differentiation is Figure 6.
