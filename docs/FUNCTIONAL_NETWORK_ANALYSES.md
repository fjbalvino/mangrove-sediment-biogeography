# Functional co-abundance network

This directory preserves the Figure 5 network workflow and its subsequent
process-context and calibration audits, using the same 191-profile functional
cohort as the current functional atlas. The calibrated chain-coherence procedure
has not been released for manuscript inference; see the decision below.
The network is separate from the Figure 6
Prevalent/Shell/Cloud regional-differentiation analysis. Historical network
counts from 194 profiles are retained only as comparisons, never as targets for
the new inference.

## Required inputs

- `functional_KO_counts_original.tsv.gz`: original matrix with 194 columns; the
  accession manifest selects the current 191 profiles and excludes
  SRR31095011, SRR31095534 and SRR31095535.
- `metadata_functional.csv`: curated metadata, joined by `run_accession`.
- `KO_persistence_classes.tsv`: verified atlas classification. Script 10H
  recomputes it from counts and metadata and requires exact agreement in KO
  identifiers, class, occurrence count and total annotation count.
- `137_KO_function_dictionary_used.tsv`: frozen dictionary with 28,352 KOs.
- `144_KO_to_KEGG_MODULE_long.tsv`: historical mapping with 247 links for 109
  KOs. It is retained for descriptive annotation only. Its restricted coverage
  does not establish that unmapped new nodes lack KEGG MODULE membership.

The companion archive supplies these inputs and the complete scripts. Raw
matrices and fitted R objects are not included in the code repository.

## Analysis sequence

`10H → 10B → 10A → 10C → 10D → 10E → 10F → 10J → 10K → 10G`

1. Verify the cohort and persistence classification. The 14,939 peripheral KOs
   are subdivided into multi-realm non-core and realm-specific classes for
   compatibility with network annotations; this does not alter the atlas's four
   primary persistence classes.
2. Retain 17,482 nonzero KOs. The association filter requires occurrence in at
   least 10 profiles and at least 10 total KO assignments (13,308 KOs). Add a
   pseudocount of 1 and center log counts within this association-filtered
   universe. The network filter requires occurrence in at least 20 profiles and
   at least 50 total assignments (10,614 eligible KOs). Retain the 5,000 highest
   CLR-variance KOs, then the top 1,500 network candidates; resolve ties by KO ID.
3. Generate 500 subsamples without replacement, capped at four profiles per
   marine realm × country × locality stratum. There are 26 strata and 81
   profiles per subsample. This caps local replication; it does not equalize
   sample numbers among realms or remove study effects.
4. Estimate proportionality as `2 cov(x,y) / (var(x)+var(y))` in CLR space.
   Retain the 2,000 strongest absolute associations per subsample. Consensus
   requires selection frequency ≥0.80, complete sign consistency and mean
   absolute proportionality ≥0.95. Discover communities with weighted Leiden clustering at resolution 1, using 50 starts (seeds 2001–2050) and convergence to a stable partition. Retain the largest modularity at this fixed resolution; ties retain the first maximum. Interpret modules with at least 10 KOs. Subsampling seed: 123. The complete 1,500-candidate graph is optimized jointly, including all connected components.
5. Annotate broad functional categories using the frozen dictionary and the
   documented historical text-matching rules. Script 10A records rule scores
   and ties. Module numbers are size ranks and do not carry fixed biological
   labels across runs.
   The driver now also calls `lib/biogeochemical_context.py` after 10E. This
   overlays the atlas's frozen KO-to-process rules on the existing partition,
   retaining many-to-many membership, optional components, alternatives and
   unavailable matrix rows. It does not modify edges, modules or the existing
   category enrichment tests. The context outputs are descriptive.
6. Evaluate enrichment using 4,999 blocked permutations. Connected-node tests
   condition on occurrence, total abundance and CLR variance; module tests
   additionally condition on degree and absolute strength. Functional
   specialization requires both Benjamini–Yekutieli q ≤0.05 and permutation
   maxT family-wise P ≤0.05. These are conditional, post-selection tests on
   KOs and graph partitions, not independent tests using metagenomes as
   replicates.
7. Use KO descriptions alone for the primary broad-category assignment. Test the complete module × seven-category family (42 tests for the current six modules), including categories with zero counts. Repeat the same family after adding the partial historical KEGG MODULE text, using exactly the same permutations. This comparison measures annotation sensitivity; incomplete pathway coverage cannot define pathway absence. Shared keyword rules are defined in `scripts/10_functional_network/lib/functional_category_rules.R`.
8. Evaluate class mixing with 999 degree-preserving rewires; assess 30
   consensus-threshold combinations, four community algorithms, and 500 edge
   perturbations. Reuse partitions for identical edge sets, separating effective threshold changes from algorithm randomness. Leiden and Louvain use 50 starts; Leiden is the primary method in every perturbation. Random edge masks are generated before parallel evaluation, making the result independent of worker count. Perturbations retain each edge with probability equal to its
   subsampling selection frequency. Co-assignment stability is averaged within
   each reference module. Community algorithms operate on the full candidate
   graph, matching primary inference; they are not optimized separately within
   each connected component.

9. Repeat edge inference after removing each metadata study unit, with 500 new locality-capped subsamples per omission. Keep the CLR reference and candidate set fixed, so these are conditional sensitivity tests. The metadata contain 43 BioProject accessions and one unresolved project label; that last omission is reported separately and is not asserted to represent a known study. Different BioProjects may belong to one publication. Script 10J verifies its numerical engine by reproducing all 1,000,000 edge selections from the primary R analysis before inference. Re-estimate each partition using the same 50-start Leiden policy. Report edge Jaccard, reference-node retention, ARI across reference major-module nodes, within-module co-assignment and best-overlap module Jaccard. Co-assignment alone does not penalize module mergers.

KEGG-process permutation tests are disabled for the supplied restricted
snapshot. A complete, versioned mapping covering the analysis universe is
required before explicitly enabling `--kegg_mapping_complete 1` in script 10F.

## Execution

The driver is `workflow/05_RUN_FUNCTIONAL_NETWORK.py`. Its required arguments
are `--counts`, `--metadata`, `--atlas-classes`, `--annotations`, and `--outdir`;
`--rscript` selects the R executable. `--workers 4` controls independent computations. A complete run through G executes the study sensitivity before rendering the figure; `--skip-study-sensitivity` reuses completed study outputs when G is requested; it does not create a main figure without those results. Restarting at G redraws the figure only. `--stop-after G` also creates the main
figure locally. All commands, exit codes and timestamps are saved under the
output directory. `--start-at` permits continuation from existing completed
stages without overwriting them; output directories have timestamps.

The inferential environment is pinned to R 4.3.3, data.table 1.17.8 and igraph
2.2.1, matching the archived runtime report. Auxiliary annotation plots are disabled in the driver (`--make_figures 0`);
the main figure remains a separate stage. Plotting additionally needs
ggplot2 and patchwork; compressed input reading uses R.utils. Python requires
the packages in `config/functional_requirements.txt`.

Use the complete portable command in [the execution guide](../workflow/00_README_HOW_TO_RUN.md). Add `--dry-run` to preview the stage plan without reading data or creating outputs. Source inputs, annotations and the output directory are explicit local arguments; `--accessions` and `--input-checksums` select the cohort and checksum configuration. The manuscript recipe retains its frozen scientific checks.

Execution logs are written beneath the selected output directory. Timestamped computation directories preserve individual run provenance. Figures are generated locally and are not uploaded by the workflow. Archived numerical summaries belong under `data/reference_results/functional_network/`, outside Git.

## Editorial delivery

Keep the main results and Figure 5 caption in `FUNCTIONAL_NETWORK_RESULTS.md`.
Keep complete statistical tables, sensitivity analyses and supplementary
methods in `FUNCTIONAL_NETWORK_SUPPLEMENTARY_RESULTS.md`. Their Word
deliverables are `Functional_network_results.docx` and
`Functional_network_supplementary_results.docx`. Replace the current files when
correcting them; use version history instead of dated or numbered copies.

Retain the accepted narrative structure and the reviewer's text where the
current evidence supports it. Keep decisive limitations in the main text and
full inferential detail in the supplement. Editorial changes must not silently
change cohorts, definitions or analytical results. Update the same companion
archive when the documents change. The companion package includes local figures and illustrated Word reports. GitHub publication remains text-only; images and binary reports are excluded.

## Figure and report regeneration

The archived stage-G figure below has the earlier three-panel organization.
The current four-panel Figure 5 is assembled with `10Q --figure5-from` after the
expanded enrichment has been plotted; see “Current four-panel Figure 5” below.
Its biological annotation replaces the old enrichment and membership-sensitivity
panels. The historical report prose is not automatically updated by this figure
operation and must not supply the new panel legend.

The driver verifies all five source-input SHA256 checksums before execution. Panel a labels every module with its current functional description and KO count. Labels are generated from the tested category table; modules without enrichment are labelled Diverse functions. Panel b contains the supported primary categories; an open point and asterisk flag annotation-sensitive support. Panel c contrasts mean co-assignment across 500 edge perturbations with the median and observed range across 43 BioProject omissions. The unresolved project-label unit is reported separately in the supplement. Stage G validates plotted nodes and edges against the current inference, without imposing historical numbers, and exports the three-panel figure and individual panels at 600 dpi plus vector PDF.

To regenerate only the figure after a completed analysis, supply the same input arguments and add `--start-at G --stop-after G`. Rebuild both Word reports with:

```bash
python tools/build_functional_network_documents.py --repository . --figure /path/to/main_figure/figures/177_reproducible_functional_coabundance_structure.png --outdir /path/to/reports
```

Word text uses Arial 11, black, with bold titles only. The main report embeds Figure 5; the separate supplementary report retains all 42 primary category tests, the informative annotation comparison, and grouped threshold results. Complete sensitivity tables remain in TSV files.

The external archived summaries omit per-omission accession lists and host-specific execution paths. Complete sample manifests remain in the companion archive and are regenerated locally. The network is Figure 5 in the current manuscript; the regional workflow retains its historical filenames for compatibility with archived outputs.

## Biochemical interpretation and redesign pilot

The biological target is reproducible covariation among encoded functions in
their reaction and pathway context. A statistical module and a KEGG biochemical
module are different objects. Do not assign the name of a complete cycle to a
statistical module simply because it contains one shared enzyme.

The annotation layer uses the existing 65-process, 97-rule atlas catalogue. It
reports every KO's process membership, every logical block's representation in
each statistical module, and genes present in the matrix but outside that module.
KOs outside the targeted catalogue are now covered by the repertoire-wide
reference mapping described below; they are not classified as biologically
unknown. Shared KOs remain explicitly
many-to-many and cannot independently confirm all listed pathways.

### Repertoire-wide functional annotation

`10O_EXTEND_FUNCTIONAL_ANNOTATION.py` extends interpretation beyond the selected
biogeochemical cycles. It maps **every original matrix KO**, including zero-only
rows and KOs excluded by network selection, to a frozen KEGG reference. This is
identifier mapping, not sequence reannotation: the historical descriptions,
counts and sample accessions remain unchanged. Current descriptions are added
in a separate field, and changed or no-longer-listed entries are reported.

The reference comprises KO names, PATHWAY and REACTION cross-references, the
complete KO hierarchy, all detailed protein/RNA-family hierarchies listed in its BRITE branch,
and all MODULE definitions listed at retrieval. This captures central metabolism,
amino-acid/cofactor/lipid/nucleotide metabolism, transport, information processing,
repair, protein folding, signaling, motility and other recorded functions. A
downloaded hierarchy is not proof of its biological relevance to these sediments.
Organismal/disease labels remain visible as reference context and cannot establish
a host, pathogenicity, resistance phenotype or clinical risk.

All memberships are retained. Descriptive tables count each KO once per term;
overlapping term percentages need not sum to 100. Detailed BRITE paths preserve
the distinctions among transporter subfamilies and other protein families.
Reaction equations are reference chemistry, not observed flux or direction.
An unclassified KO, a named KO outside a pathway, an explicitly uncharacterized
protein and an identifier missing from the current reference remain distinct.

MODULE definitions use the atlas's tested AND/OR/optional parser, with nested
module expansion. Each definition is audited. Unsupported definitions, including
unassigned `--` steps, remain unevaluated and are exported with the reason; they
are never declared absent. Optional-only blocks are excluded from required-block
coverage. Rules are evaluated separately within each community at thresholds of
one and two annotation records. An unavailable matrix row is unknown, not zero.
Meeting a community rule does not establish a complete pathway in one organism.

The same logical blocks are overlaid on existing statistical modules, recording
genes inside the module, measured genes elsewhere and unavailable rows. Existing
edges receive shared pathway/module context without changing the graph. Broad
global/overview maps and organismal/disease maps are excluded from the shared
biological-pathway edge field. Shared membership does not demonstrate a direct
reaction, adjacency, interaction or causal direction. Statistical modules are
never renamed automatically after the most frequent pathway label. For example,
transporters linked to the quorum-sensing map retain their transport functions;
that map membership alone is not evidence of communication.

This annotation expansion is **descriptive**. It does not reuse the historical
seven-category P values or enable the old KEGG enrichment switch. The separate
conditional tests below use their own declared family and comparison universes.
The frozen atlas rules and the eight-chain calibration remain unchanged.

Retrieve a reference explicitly (academic KEGG API use, no more than two request
starts per second), then run offline annotation:

```bash
python scripts/10_functional_network/10O_EXTEND_FUNCTIONAL_ANNOTATION.py fetch-reference \
  --outdir /path/to/kegg-reference
python scripts/10_functional_network/10O_EXTEND_FUNCTIONAL_ANNOTATION.py annotate \
  --matrix /path/to/functional_KO_counts_original.tsv.gz \
  --samples config/samples/functional_accessions.txt \
  --nodes /path/to/network/final_results/173/tables/173_network_nodes_annotated.tsv \
  --candidates /path/to/network/final_results/171/tables/171_candidate_KO_selection.tsv \
  --edges /path/to/network/final_results/171/tables/171_consensus_edges.tsv.gz \
  --dictionary /path/to/annotations/137_KO_function_dictionary_used.tsv \
  --snapshot /path/to/kegg-reference --outdir /path/to/broad-functional-context
```

For the full network driver, add `--functional-reference /path/to/kegg-reference`.
It invokes this same canonical annotation after stage E, without live downloads
or changes to inference. Omitting the argument preserves historical replay.
Keep the raw reference, manifests and generated tables in the external companion,
never Git. Reuse a verified snapshot for repeated analyses; deliberately fetching
a new reference requires a new external snapshot and an annotation comparison.

Primary documentation: [KEGG API](https://www.kegg.jp/kegg/rest/),
[API manual](https://www.kegg.jp/kegg/rest/keggapi.html),
[BRITE](https://www.kegg.jp/kegg/brite.html), and
[MODULE logic](https://www.kegg.jp/kegg/module.html).
Database citation: Kanehisa et al. (2025), *KEGG: biological systems database as
a model of the real world*, *Nucleic Acids Research* 53:D672–D677,
[doi:10.1093/nar/gkae909](https://doi.org/10.1093/nar/gkae909).

The executed 191-profile audit covered 17,508 matrix KOs, 1,500 candidates and
all 327 connected nodes, while preserving the 1,062 historical edges. Of these
nodes, 306 have hierarchical/family context and 274 classified nodes lie outside
the atlas's selected catalogue. The reference includes 56 detailed BRITE
hierarchies and 583 MODULE definitions; 569 definitions were evaluated and 14
with unsupported unassigned steps were retained without completion scores.
The 217,358 per-sample module records reconcile with all summary counts.
Within historical M06, 11 of 24 KOs belong to transporters, including five
branched-chain amino-acid transport components. Within M01, 17 of 57 KOs map
to aminoacyl-tRNA biosynthesis. These are descriptions of module contents,
not new enrichment or activity claims. Full results and frozen reference bytes
are in the external `Functional_network_pilot.zip` companion under
`broad_functional_context/`.

### Conditional enrichment of the expanded annotation

`10P_ENRICH_NETWORK_FUNCTIONS.py` tests the existing six major modules without
changing nodes, edges or community assignments. The primary comparison universe
is their union (236 KOs); sensitivity includes all 327 connected KOs and assigns
the other 91 nodes to an outside label. The original 1,500 candidates include
disconnected nodes, which have no comparable connected-module assignment and are
not the primary reference. Selection already depends on detection/count filters,
CLR variance and stable associations; the graph is not simply the most abundant
genes or an exhaustive definition of the functional core.

Terms comprise biological pathway categories, nonoverview pathways, BRITE
families and detailed subfamilies, and nonoptional MODULE memberships. Eligibility
requires at least three member and three nonmember KOs in the union of major
modules. Filtering uses that union, never a specific module's enrichment result.
All eligible terms are crossed with all major modules, including zero-count
cells. Overlapping sets and multiple annotations are retained. In this release,
150 sets yield 900 tests in each comparison universe.

KOs are grouped by prevalence, log10(total annotation records + 1), CLR variance,
degree and absolute weighted strength. Median/MAD scaling is followed by
scikit-learn KMeans (Lloyd, 50 starts, seed 20261105, 200 iterations), targeting
35 KOs per block and merging blocks smaller than 12 into the nearest centroid.
Matching does not use functional annotations or module labels. It conditions
on coarse covariate blocks, not exact covariate equality. Block membership and
within-block covariate ranges are exported so this approximation can be audited.

The null exchanges module labels within blocks while preserving block-specific
module sizes and complete annotation vectors. For each term/module pair, the
one-sided count probability is computed by convolution of independent
stratum-specific hypergeometric distributions. Benjamini–Yekutieli adjustment
uses the entire 900-test family. As an additional joint-null safeguard, common
label permutations produce the maximum standardized enrichment across all tests.
The empirical probability is (exceedances + 1)/(permutations + 1). A reported
support flag requires both adjusted probabilities at most 0.05. The maxT
calculation is a guard under the stated joint exchangeability null; no strong
familywise-error guarantee under arbitrary mixtures of alternatives is claimed.

A 9,999-draw computational pilot identified Monte Carlo uncertainty near 0.05.
One documented precision amendment increased both complete families to 99,999
draws while retaining all hypotheses, blocks, exact tests and thresholds. Wilson
intervals for Monte Carlo uncertainty are exported; no further outcome-dependent
extension is used. Exact P values and BY adjustments do not depend on this draw
count. Synthetic exhaustive-allocation tests independently check the conditional
distribution, and empirical null means are checked against its exact expectation.

These are **exploratory competitive annotation tests conditional on the selected
graph**. They are not independent confirmation of network discovery, sample-level
differential abundance, or evidence of reaction direction, activity or flux.
The earlier process-coherence calibration tests a different question and does
not validate these enrichment tests. Original matrix values and annotation
records are never relabeled as transcript abundance or activity.

```bash
python scripts/10_functional_network/10P_ENRICH_NETWORK_FUNCTIONS.py \
  --annotation /path/to/broad-functional-context \
  --nodes /path/to/network/final_results/173/tables/173_network_nodes_annotated.tsv \
  --candidates /path/to/network/final_results/171/tables/171_candidate_KO_selection.tsv \
  --edges /path/to/network/final_results/171/tables/171_consensus_edges.tsv.gz \
  --outdir /path/to/functional-enrichment
python scripts/10_functional_network/10Q_PLOT_FUNCTIONAL_ENRICHMENT.py \
  --enrichment /path/to/functional-enrichment \
  --annotation /path/to/broad-functional-context \
  --matrix /path/to/functional_KO_counts_original.tsv.gz \
  --samples config/samples/functional_accessions.txt \
  --outdir /path/to/functional-enrichment/figures
```

The enrichment output directory must be empty. The plotting matrix must match
the annotation source hash. Generated tables, figures, reference bytes and
permutation draws belong in the external companion, never Git.

The current primary test supports five term/module combinations, all in M01 or
M06; three combinations persist with the connected-node reference. Two of those
three are the same observed KO set: translation and aminoacyl-tRNA biosynthesis
both contain 17 M01 KOs. They are one biological signal, not independent findings.
The other retained combination is M06 membership of KEGG's quorum-sensing map:
all 11 contributing genes are transporters. This does not support a claim of
quorum sensing, and the separately tested broad transporter category does not
meet the joint criterion. M03's seven NADH–quinone oxidoreductase components
remain descriptive. Oxidative phosphorylation in M03 passes BY alone but not
the joint maxT guard; it is not marked as supported.

`10Q` displays KO counts, module proportions and primary support together, keeping
unstarred functions visible. All primary supported terms are shown; identical
memberships can share a row with an exported alias record. Selected respiratory
and branched-chain amino-acid transport components are shown separately from
statistical support. Per-profile rules preserve KEGG alternatives, including
NuoB/C/D fusion options. Community co-detection does not establish a complete
complex within one organism. For a new cohort, rerun annotation, inference and
the declared tests; do not transfer existing support labels.

Implementation sources: [SciPy hypergeometric distribution](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.hypergeom.html),
[SciPy false-discovery adjustment](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html),
[scikit-learn 1.8 KMeans](https://scikit-learn.org/1.8/modules/generated/sklearn.cluster.KMeans.html),
[R multiplicity documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html),
[Phipson and Smyth (2010)](https://doi.org/10.2202/1544-6115.1585), and
[KEGG MODULE M00144](https://www.kegg.jp/entry/M00144).

### Current four-panel Figure 5

The same `10Q_PLOT_FUNCTIONAL_ENRICHMENT.py` now assembles the manuscript figure
from its previously verified plotted tables and the existing stage-G graph
coordinates. It does not rerun inference, enrichment, layout optimization or
community detection. Panel a preserves all 236 major-module nodes, their
coordinates, persistence symbols and the same 999 displayed edges. Module IDs
and KO counts replace historical functional labels so that labels do not imply
support from the superseded category analysis. Panel b is the expanded functional
heatmap; panels c and d show respiratory and transport components. Component
colors match the graph's M03 and M06 colors; the heatmap retains plasma.

```bash
python scripts/10_functional_network/10Q_PLOT_FUNCTIONAL_ENRICHMENT.py \
  --figure5-from /path/to/functional-enrichment/figures \
  --network-nodes /path/to/network/final_results/177/tables/177_network_nodes_used.tsv \
  --network-edges /path/to/network/final_results/177/tables/177_network_edges_used.tsv.gz \
  --outdir /path/to/manuscript-figure-5
```

Outputs are `Figure_5_functional_coabundance.png` (600 dpi), its vector PDF and
`Figure_5_provenance.json`. Module identities, sizes, input identifiers and
heatmap percentages are checked before rendering. The companion archive contains
the current figure, its caption and the exact graph inputs; none are uploaded
to GitHub. The older regional Figure_5 files belong to the current manuscript's
Figure 6 and must not be overwritten with this network figure.

**Figure 5 | Functional organization of the co-abundance network.**
(a) Six major modules from the 191-metagenome analysis, showing 236 KOs and
999 associations; the full network contains 327 connected KOs and 1,062
associations. Colors identify modules; labels give module IDs and KO counts.
Squares, circles and triangles indicate absolute-core, near-global and realm-core
KOs, respectively. (b) Selected functional annotations across modules. Numbers
are distinct KOs; shading gives their percentage of the module, not abundance.
Annotations overlap. Asterisks indicate primary conditional enrichment meeting
both Benjamini–Yekutieli adjusted P ≤ 0.05 and joint-null maxT P ≤ 0.05 across
150 eligible sets and six modules (900 comparisons), using exact blocked count
probabilities and 99,999 common label permutations. Translation and tRNA
aminoacylation share the same observed KO set and are displayed once. The M06
quorum-sensing map membership is supported by transporter genes and is not
evidence of communication. (c) NADH–quinone oxidoreductase components, with
M03 members in green and other detected components in gray. Alternative NuoB/C/D
fusion architectures are indicated. (d) Five branched-chain amino-acid transport
components in M06, grouped by their reference roles. Box fractions give profiles
with at least one annotation record. Co-detection is community-level evidence;
it does not establish same-organism completion, activity or flux. Panels c and d
describe component membership, not inferred reaction order or physical structure.

`config/gene_biogeography/network_reference_chains.json` adds eight selected
reference chains spanning nitrogen, acetate, methane/C1 and sulfur. These are
biochemical context links, never additional inferred edges. Successive chemical
steps are curated explicitly; their order is not guessed from whitespace in a
logical expression. Electron-transfer or recycling machinery is kept separate
from the carbon-transfer chain. All eight atlas panels receive rule/context
coverage, but this initial chain file does not claim exhaustive ordered reaction
coverage of every process. The curated rules remain unchanged.

Source definitions: [KEGG module semantics](https://www.kegg.jp/kegg/module.html),
[nitrogen-oxide reduction](https://www.kegg.jp/entry/M00529),
[ammonium-producing nitrate reduction](https://www.kegg.jp/entry/M00530),
[nitrate assimilation](https://www.kegg.jp/entry/M00531),
[CO2-linked methane pathway](https://www.kegg.jp/entry/M00567),
[acetate-linked methane pathway](https://www.kegg.jp/entry/M00357), and
[sulfate-pathway context](https://www.kegg.jp/entry/M00596).
Reference direction is not measured direction. Nar/Nxr, Mcr, Apr/Dsr and
reversible acetate reactions retain their existing ambiguities.

### Diagnostic execution

Run this before changing the final network inference. Inputs are the accepted
atlas matrix/metadata and the archived network's current 171/173 outputs. The
pilot requires an empty output directory and matches samples only by accession.

```bash
python scripts/10_functional_network/10L_AUDIT_PROCESS_CONTEXT.py \
  --matrix /path/to/atlas/inputs/functional_KO_counts_original.tsv.gz \
  --metadata /path/to/atlas/inputs/metadata_functional.csv \
  --samples config/samples/functional_accessions.txt \
  --reference-nodes /path/to/network/final_results/173/tables/173_network_nodes_annotated.tsv \
  --reference-edges /path/to/network/final_results/171/tables/171_consensus_edges.tsv.gz \
  --reference-candidates /path/to/network/final_results/171/tables/171_candidate_KO_selection.tsv \
  --reference-design /path/to/network/final_results/171/tables/171_bootstrap_design.tsv.gz \
  --outdir /path/to/network-process-pilot
```

The historical candidate identifiers and CLR variances must reproduce before
the pilot proceeds. This guard compares the supplied reference with its source
matrix; it does not impose the old candidate set on Dani's future matrix. New
data first require their own reference run and a new calibration.

The pilot compares run-level global proportionality, equal-study weighted
study-coordinate proportionality, and within-study coordinate deviations. CLR
coordinates are averaged within BioSample labels, then labels within each
study-coordinate group; each study has equal total weight and its coordinate
groups share that weight equally. This is a mean of CLR coordinates, not the
atlas's arithmetic mean of per-million annotations. Coordinates are grouped
at five decimal places for consistency with the atlas, without implying that
the original coordinates have that precision. Studies with one coordinate
group supply no within-study coordinate contrast. These are distinct descriptive
estimands, not alternative tests of the same null.

Sensitivity uses pseudocounts 0.5, 1 and 2 in original annotation-record units;
the original association reference is held fixed for edge-weight comparisons.
A separate comparison uses all detected matrix rows as CLR reference. Feature
selection is independently repeated for each pseudocount, with and without the
influential project. The exact archived 500 subsets also provide unconditional
weights for the same historical edges. Their percentiles are sampling
sensitivity ranges, not confidence intervals. The custom proportionality
implementation is checked against the independent log-ratio-variance identity;
the pilot does not claim execution of the R `propr` package.

### Verified findings from the 191-profile pilot, 2026-11-04

- Reproduced all 1,500 historical candidate identities and their CLR variances.
  All 1,062 historical consensus edges were positive.
- Of 368 queried atlas KOs, 328 have matrix rows, 125 enter the historical
  candidate set and only 19 occur in the six major modules (236 KOs).
  No nitrogen-catalogue KO is connected in the historical network. Sulfur has
  two connected KOs, neither in a major module. Cycle-specific counts overlap
  when KOs have multiple process memberships.
- The selected Nap → NirK/S → cNOR → NosZ reference chain is co-detected in
  153/191 profiles at one record per component and 141/191 at two. This is a
  community-level, architecture-specific statement, not complete denitrification
  in an organism, expression, net nitrogen loss or evidence of a statistical edge.
- The influential project's six profiles have median total assignments 12,255.5
  and median detected KOs 2,018, versus 389,810 and 8,372 in the other 185.
  Those contrasts identify a recovery-related explanation to test; they do not
  demonstrate its cause or justify excluding the project.
- There are 85 study-coordinate groups in 44 provenance labels. Only 10 labels
  have multiple coordinates (51 groups). These labels are not certified
  independent studies. Median proportionality of historical edges is 0.978
  globally, 0.974 under study-coordinate weighting and 0.829 within studies.
  At the historical diagnostic cutoff 0.95, 123/1,062 edges remain above the
  cutoff within studies, and 13/1,062 after omitting the influential project
  (9 labels, 46 groups). These are not significant-edge counts or FDR estimates.
- Omitting the project leaves all historical edges above 0.95 in the run-level
  global calculation. Its effect on the old selected topology therefore cannot
  be read simply as disappearance of strong global covariation: ranking,
  the per-draw edge budget and consensus selection also matter.
- Pseudocount changes alter candidate identity (Jaccard 0.717 for 0.5 and
  0.803 for 2, relative to 1). The conditional-minus-unconditional mean weight
  difference is small in this run (mean 0.000441; maximum 0.003945).

The pilot exports three PNG/PDF figures, complete tabular evidence and checksums
outside Git. Input matrices and metadata are not modified. It does not reinfer
final modules, run null-network simulations, establish edge FDR, validate held-out
studies, test new process enrichment or replace manuscript results. Those stages
must follow the sampling and annotation diagnosis. Preserve the distinction
between the existing six modules and any subsequently inferred partition.

The next calibration must test both a discovery network and a prespecified
process-focused coherence analysis. The latter retains measured members of the
curated chains even when they fail the discovery network's variability ranking;
it must use its own declared hypothesis family and must not add edges merely
because a biochemical reference connects two steps. Module interpretation will
then report which steps are represented together, split between modules,
outside the selected network, or unavailable in the matrix.

## Calibration of process coherence

`10M_CALIBRATE_PROCESS_COHERENCE.py` evaluates the eight reference chains fixed
by the preceding pilot. It is a separate analysis of measured markers and does
not add biochemical edges to the discovery network. It does not calibrate the
complete historical candidate-selection, top-edge-budget and clustering pipeline.

For each adjacent pair of reference steps, calculate proportionality for every
pair of available KOs spanning those steps. Average pairs equally within a
transition and transitions equally within a chain. Alternative enzymes remain
in the same step; their mutual association is not required. Shared KOs across
steps cause an explicit stop rather than a self-edge. A missing step is not
replaced with zero evidence. Two invariant measured markers contribute zero
coherence, not perfect association. This score describes coordinated relative
representation of measured markers, not enzyme abundance, flux or completeness.

Use the original association-filtered CLR reference with one annotation record
as pseudocount. Average runs within BioSample labels, then labels within each
study-coordinate group. Centre coordinates within each provenance label.
For each KO pair, its proportionality is twice the mean within-label covariance
divided by the sum of the mean within-label variances. Labels receive equal
covariance weight; labels with fewer than two coordinates provide no contrast.
The chain score decomposes into one contribution per informative provenance label.

Enumerate all sign flips of these block contributions for a one-sided test of
positive coherence. The exact test assumes independent blocks and joint sign
symmetry under the null, a stronger condition than zero average covariance.
For a chain, reversing the signs of alternating steps reverses all adjacent-step
covariances while preserving variances. This explains the transformation under
the stated symmetry assumptions; the reference arrows do not supply statistical
replicates. Within-block spatial correlation is retained by a whole-block flip.
BioProject labels are provenance labels, not verified independent publications.
Apply Holm correction to all eight chains, retaining unmeasurable chains with
P=1. Scores, block contributions and assumption-conditional P values are saved
regardless of the subsequent release decision.

The fixed plan in `config/network_calibration.json` specifies 1,000 independently
seeded replicates in each of eight scenarios. Three statistical nulls use
independent Gaussian genes, independent spatial Gaussian gene fields, or symmetric
heavy-tailed genes with unequal group precision. Three measurement stress tests
start with independent latent gene intensities and simulate Poisson annotation
recovery at equal yields, the observed per-run yields, or observed yields plus
20% independent coordinate-level dropout. Two positive controls introduce a
shared factor with latent correlation 0.6 among the Nap-branch nitrogen markers.
Shared membership can transmit that signal to another chain; power is reported
for the prespecified target rather than declaring every other chain null.

The count generator retains the observed run/BioSample/study-coordinate hierarchy
and all matrix rows. Study mean relative profiles are mixed with the equal-study
global mean at weights 0.9/0.1. Independent lognormal effects have log SD 0.75 at
coordinate level and 0.35 at BioSample level, with mean correction before Poisson
sampling. These are declared scenarios, not fitted explanations of real gene
recovery or a reconstruction of assembly. Recovery totals are annotation records,
not read counts. The Gaussian generator uses observed study means and pooled
within-study marker scales, with a 0.05 scale floor; the heavy-tailed scenario
divides t3 innovations by sqrt(3) and by sqrt(number of BioSample labels per group).
The spatial kernel is 0.8 times an exponential chord-distance kernel plus 0.2
identity, with range equal to each study's median positive distance.

Statistical-null rejection rates measure type-I error under those models.
Rejection rates in the independent-intensity count scenarios instead diagnose
apparent biological coupling caused by recovery and compositional representation;
they are not automatically type-I errors for a relative-data null. Report Wilson
95% intervals for simulation rates. The declared diagnostic gate requires the
upper bound of each family rejection rate to be at most 0.075; this tolerance
addresses simulation uncertainty and does not change the testing alpha of 0.05.
Failure of either statistical calibration or the recovery stress check prevents
promotion of observed P values to manuscript evidence. Passing is conditional,
not proof that the real independence and recovery assumptions hold.

Sensitivity retains the same eight-chain family and changes pseudocounts to
0.5 or 2, uses all detected KOs as reference, omits each informative provenance
label, and requires at least three coordinates per label. These are diagnostics,
not additional opportunities to choose a significant result. The same 1,062
historical edges are also screened at rho >=0.95 in simulations; those fixed-edge
counts are not a graph-wide FDR or a new network.

```bash
python scripts/10_functional_network/10M_CALIBRATE_PROCESS_COHERENCE.py \
  --matrix /path/to/atlas/inputs/functional_KO_counts_original.tsv.gz \
  --metadata /path/to/atlas/inputs/metadata_functional.csv \
  --samples config/samples/functional_accessions.txt \
  --reference-edges /path/to/network/final_results/171/tables/171_consensus_edges.tsv.gz \
  --reference-candidates /path/to/network/final_results/171/tables/171_candidate_KO_selection.tsv \
  --outdir /path/to/network-calibration --workers 4

python scripts/10_functional_network/lib/plot_calibration.py \
  --run /path/to/network-calibration
```

`--resume` checks input, source, curation and plan hashes before continuing.
Independent seeds are indexed by scenario and replicate; changing worker count
does not change results. `--replicates` permits a smoke test but a campaign below
the configured count cannot pass the release gate. Exact enumeration is supported
for up to 20 informative labels; a larger new cohort requires an explicit update
to the computational strategy, not silent approximation.

Method sources: [Lovell et al. (2015)](https://doi.org/10.1371/journal.pcbi.1004075),
[Erb and Notredame (2016)](https://doi.org/10.1007/s12064-015-0220-8),
[Winkler et al. (2014)](https://doi.org/10.1016/j.neuroimage.2014.01.060),
[Holm (1979)](https://www.jstor.org/stable/4615733), and
[Weiss et al. (2016)](https://doi.org/10.1038/ismej.2015.235).
The custom score and exact sign enumeration are tested against the independent
proportionality identity and SciPy's exact paired-sample permutation calculation.
SciPy provides the binomial Wilson intervals. The implementation uses NumPy,
pandas and SciPy; figures use Matplotlib.

### Complete edge-selector stress controls

`10N_AUDIT_NETWORK_SELECTOR_STRESS.py` reuses `10J.infer_case` rather than
introducing another consensus-edge implementation. Before controls, it replays
all 500 archived subsets and compares every selected edge with the original R
outputs. It then generates ten independent datasets with common expected
annotation yield and ten with the observed yields. All provenance labels receive
the same global baseline profile: there are neither latent cross-gene links nor
latent study shifts in these controls. The coordinate/BioSample hierarchy and
independent gene-specific variation remain. Parameters and seeds are fixed in
`config/network_selector_stress.json`.

Each simulated dataset receives a newly estimated CLR reference and candidate
set, followed by all 500 archived subsets, the 2,000-edge per-subset budget, sign
consistency and original consensus thresholds. Clustering is not needed to
evaluate whether edges survive these controls and is not executed. These controls
test susceptibility to the stated recovery model. Their null concerns latent
genes, not necessarily marginal relative-data association, and twenty controls
do not estimate graph-wide FDR or certify all possible measurement models.

```bash
python scripts/10_functional_network/10N_AUDIT_NETWORK_SELECTOR_STRESS.py \
  --matrix /path/to/atlas/inputs/functional_KO_counts_original.tsv.gz \
  --metadata /path/to/atlas/inputs/metadata_functional.csv \
  --samples config/samples/functional_accessions.txt \
  --reference-edges /path/to/network/final_results/171/tables/171_consensus_edges.tsv.gz \
  --reference-candidates /path/to/network/final_results/171/tables/171_candidate_KO_selection.tsv \
  --reference-design /path/to/network/final_results/171/tables/171_bootstrap_design.tsv.gz \
  --reference-membership /path/to/network/final_results/171/tables/171_bootstrap_edge_membership.tsv.gz \
  --outdir /path/to/network-selector-stress --workers 4
```

### Calibration result and interpretation

The 8,000-dataset chain campaign completed with 1,000 replicates per scenario.
Family rejection rates in the three statistical nulls were 3.9%, 4.3% and 5.1%,
with all Wilson upper 95% bounds below the declared 7.5% diagnostic limit.
Recovery-stress rejection rates were 4.3% at equal expected yields, 12.5% at
observed yields (95% interval 10.6–14.7%), and 9.7% with additional dropout
(8.0–11.7%). The measurement-stress gate therefore failed. The prespecified
nitrogen positive control was detected in 90.5% of Gaussian and 72.2% of recovered
count simulations; these detection rates do not offset the failed stress check.

In the observed data, none of the eight primary chain tests passed Holm 0.05:
adjusted assumption-conditional P values ranged from 0.0546875 to 0.2001953.
Some sensitivity variants crossed the threshold, which is not a reason to
replace the primary specification. The scores remain descriptive and their
tests are not promoted to biological evidence. Within-study coherence answers
a different question from geographic differences between studies; lack of
support here does not establish absence of a pathway or absence of biogeography.

The full historical selector reproduced all 1,000,000 archived edge selections
and all 1,062 consensus edges; the maximum mean-weight discrepancy was
`1.9984e-15`. None of the twenty complete-selector controls produced a consensus
edge (0/10 under each recovery condition). This limited result provides no basis
for declaring the historical graph an artifact of recovery alone. It also does
not resolve its provenance sensitivity, validate its communities, or make its
modules equivalent to biogeochemical chains.

**Decision:** retain the biochemical annotation and reaction-context layer.
Keep co-detection, relative representation and statistical co-variation as
separate evidence. Do not add cycle-coupling claims or a newly validated network
to the manuscript from this campaign. The historical graph can be examined as
an exploratory global co-abundance pattern; its biological labels must follow
the actual annotated members. Re-run the declared diagnostics and calibration
for the next input release instead of carrying forward a significance label.
