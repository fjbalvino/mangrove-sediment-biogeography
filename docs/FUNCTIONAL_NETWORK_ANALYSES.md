# Functional co-abundance network

This workflow updates the Figure 6 network using the same 191-profile functional
cohort as the current functional atlas. It is separate from the Figure 5
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

`10H → 10B → 10A → 10C → 10D → 10E → 10F → 10G`

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
3. Generate 100 subsamples without replacement, capped at four profiles per
   marine realm × country × locality stratum. There are 26 strata and 81
   profiles per subsample. This caps local replication; it does not equalize
   sample numbers among realms or remove study effects.
4. Estimate proportionality as `2 cov(x,y) / (var(x)+var(y))` in CLR space.
   Retain the 2,000 strongest absolute associations per subsample. Consensus
   requires selection frequency ≥0.80, complete sign consistency and mean
   absolute proportionality ≥0.95. Discover communities with weighted Louvain
   clustering; interpret modules with at least 10 KOs. Seed: 123.
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
7. Evaluate class mixing with 999 degree-preserving rewires; assess 30
   consensus-threshold combinations, four community algorithms, and 500 edge
   perturbations. Perturbations retain each edge with probability equal to its
   subsampling selection frequency. Co-assignment stability is averaged within
   each reference module. Community algorithms operate on the full candidate
   graph, matching primary inference; they are not optimized separately within
   each connected component.

KEGG-process permutation tests are disabled for the supplied restricted
snapshot. A complete, versioned mapping covering the analysis universe is
required before explicitly enabling `--kegg_mapping_complete 1` in script 10F.

## Execution

The driver is `workflow/05_RUN_FUNCTIONAL_NETWORK.py`. Its required arguments
are `--counts`, `--metadata`, `--atlas-classes`, `--annotations`, and `--outdir`;
`--rscript` selects the R executable. `--stop-after G` also creates the main
figure locally. All commands, exit codes and timestamps are saved under the
output directory. `--start-at` permits continuation from existing completed
stages without overwriting them; output directories have timestamps.

The inferential environment is pinned to R 4.3.3, data.table 1.17.8 and igraph
2.2.1, matching the supplied Landa runtime report. Auxiliary annotation plots are disabled in the driver (`--make_figures 0`);
the main figure remains a separate stage. Plotting additionally needs
ggplot2 and patchwork; compressed input reading uses R.utils. Python requires
the packages in `config/functional_requirements.txt`.

For the standalone companion package on Landa:

```bash
cd /home/fjbalvino/Biogeografia_2026
unzip -n MAGENTA_red_191_20261015.zip
bash MAGENTA_red_191/ejecutar_en_landa.sh
```

The wrapper uses `/data/env/dastool/bin/Rscript` unless `RSCRIPT` is explicitly
set, and writes a new timestamped run beneath
`/home/fjbalvino/Biogeografia_2026/MAGENTA_red_191/reanalisis/`. Figures are
generated locally by the scripts and are not uploaded by this workflow.
