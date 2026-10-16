# Functional co-abundance network

This workflow updates the Figure 5 network using the same 191-profile functional
cohort as the current functional atlas. It is separate from the Figure 6
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
2.2.1, matching the supplied Landa runtime report. Auxiliary annotation plots are disabled in the driver (`--make_figures 0`);
the main figure remains a separate stage. Plotting additionally needs
ggplot2 and patchwork; compressed input reading uses R.utils. Python requires
the packages in `config/functional_requirements.txt`.

For the standalone companion package, extract or refresh the same package directory:

```bash
unzip -o MAGENTA_functional_network.zip
bash MAGENTA_red_191/ejecutar_en_landa.sh
```

The wrapper uses `Rscript` unless `RSCRIPT` is explicitly set and writes its
execution record beneath `MAGENTA_red_191/reanalisis/`. Timestamped computation
directories preserve individual run provenance; they are not parallel versions
of the published manuscript or scripts. Figures are generated locally and are
not uploaded by this workflow.

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

The driver verifies all five source-input SHA256 checksums before execution. Panel a labels every module with its current functional description and KO count. Labels are generated from the tested category table; modules without enrichment are labelled Diverse functions. Panel b contains the supported primary categories; an open point and asterisk flag annotation-sensitive support. Panel c contrasts mean co-assignment across 500 edge perturbations with the median and observed range across 43 BioProject omissions. The unresolved project-label unit is reported separately in the supplement. Stage G validates plotted nodes and edges against the current inference, without imposing historical numbers, and exports the three-panel figure and individual panels at 600 dpi plus vector PDF.

To regenerate only the figure after a completed analysis, supply the same input arguments and add `--start-at G --stop-after G`. Rebuild both Word reports with:

```bash
python tools/build_functional_network_documents.py --repository . --figure /path/to/main_figure/figures/177_reproducible_functional_coabundance_structure.png --outdir /path/to/reports
```

Word text uses Arial 11, black, with bold titles only. The main report embeds Figure 5; the separate supplementary report retains all 42 primary category tests, the informative annotation comparison, and grouped threshold results. Complete sensitivity tables remain in TSV files.

Repository summaries omit per-omission accession lists and host-specific execution paths. Complete sample manifests remain in the companion archive and are regenerated locally. The network is Figure 5 in the current manuscript; the regional workflow retains its historical filenames in this network-only update.
