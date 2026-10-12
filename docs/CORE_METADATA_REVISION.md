# Core analysis with corrected sample metadata

The corrected core analysis uses 343 metagenomes from 48 sampling localities
in five marine realms. These are the 345 profiles selected for taxonomic
diversity, excluding `SRR31131347` and `SRR31131364` because their marine realm
could not be resolved from the available source metadata. Both profiles remain
in the full corrected 405-profile object and the selected 345-profile object.
Missing realms are not imputed or treated as an additional realm.

The metadata revision corrects source-verified sample locations and regional
assignments. Counts, taxonomic annotations and within-sample read totals are
unchanged. Each supplied profile already contains 18,145 retained reads;
these inputs are standardized taxonomic counts, not raw sequencing libraries.
The fixed sample identifiers in
`config/samples/taxonomic_core_assigned_realm_identifiers.txt` match the
`run_accession` field exactly, including project sample identifiers.

## External inputs

Count objects and metadata are distributed separately from this code. The
runner requires the corrected 343-profile object and the corrected metadata;
the latter may include all 405 profiles because samples are joined by
`run_accession`. Do not combine a corrected manifest with the older metadata.

| Input | Default path | SHA-256 of the validated file |
|---|---|---|
| Core count object, 343 profiles | `data/taxonomic_core_counts.rds` | `f06995e70b5d7988ba5d67345021aedb491915f1214654248a003ee4d3c5a08b` |
| Corrected metadata used by the analyses | `data/environmental_metadata.csv` | `18d289f79dd55db18885957de9b1ca6d631baf2ff6a92300016f5f939f8016f1` |
| Full corrected count object, 405 profiles | Supplied separately | `edc6f2e631733854144236dabc7b9e4c081283196e56e59c5aa3c7f36adccd31` |
| Selected corrected count object, 345 profiles | Supplied separately | `577086ca6a0d45d4f343b09e67ecefaa1a53a8a8381ae7ddee837e425295329e` |

Checksums identify file bytes; reserialization can change a checksum without
changing its count matrix. The count objects preserve the original feature
space. Stage 05A separately harmonizes abbreviated taxon labels only when an
unambiguous full-lineage match exists and records every merge.

## Run stages 05A–D

Install the R dependencies listed in `config/R_packages.tsv`, supply the two
required input files above, and run from the repository root:

```bash
bash scripts/05_taxonomic_core_biogeography/run_taxonomic_core.sh
```

`MAGENTA_CORE_PHYLOSEQ`, `MAGENTA_METADATA_CSV` and `CORE_SAMPLE_MANIFEST`
override the count, metadata and manifest paths. `BIOGEOGRAPHY_BASE` sets the
project directory, `BIOGEOGRAPHY_RESULTS` sets the output directory, and
`RSCRIPT_BIN` selects the Rscript executable. The default output directory is
`results/taxonomic_core/`. The runner stops on failure and passes the completed
05A classification to the subsequent stages through its output pointer.
The older `MAGENTA_PHYLOSEQ_401` variable remains a compatibility alias in
05A; the new runner uses `MAGENTA_CORE_PHYLOSEQ` explicitly.

| Stage | Output |
|---|---|
| 05A | Core membership, eligibility audit, realm/locality counts and core read/richness proportions |
| 05B | Taxa detected in only one realm, with separate all-detected and eligible denominators |
| 05C | Phylum composition of the core, including pooled and locality-balanced proportions |
| 05D | Core/non-core similarity against geographic and climatic distance; 9,999 matrix permutations per gradient |

The 343-profile manifest is validated against the input object before analysis.
Each stage records its inputs, parameters and session information. Intermediate
tables retain their established filenames to preserve downstream compatibility.

## Definitions and denominators

Counts are pooled within each sampling locality. The code uses a joint
`REALM` and `Locality` key so identically named localities in different realms
cannot be merged silently. Global eligibility requires at least 50 retained
reads and detection in at least five locality units.

An eligible taxon is core-associated if it occurs in every sampled locality
within at least one realm represented by three or more localities and has
mean relative abundance at least 0.01% within that realm. The relative abundance
at each locality uses all retained reads at that locality. The core is the
union of the taxa meeting these criteria; it does not require presence in
every locality worldwide or imply temporal persistence.

Figure 3a uses only globally eligible taxa for both denominators: eligible
richness and reads assigned to eligible taxa. In the corrected primary analysis,
2,773 of 10,189 eligible features are core-associated; their 3,830,379 reads
represent 62.107% of 6,167,390 eligible reads. The same core reads represent
61.545% of all 6,223,735 retained reads. Features excluded by the eligibility
thresholds are reported separately and are not labeled eligible non-core.

Locality-balanced phylum composition is calculated by pooling core counts
within each locality, expressing each phylum as a proportion of that locality's
core reads, and averaging those proportions with equal weight across the 48
localities. Thus, localities with more metagenomes do not receive more weight.
This differs from the pooled eligible-read proportion in Figure 3a.

Realm restriction is assessed among all detected features independently of core
eligibility. Detection in only one sampled realm is not evidence of endemism.
The exported tables identify which restricted taxa pass global eligibility.

## Analysis status

Stages 05A–D were rerun with the corrected metadata and the 343-profile
manifest. A broader 399-profile analysis, obtained by excluding the same two
unresolved realms from the previous 401-profile cohort, serves only as a
sensitivity reference. Its results must not be substituted for the primary
343-profile analysis.

The earlier Figure 2 diversity, environmental GAM, ordination, distance-decay
and regional-comparison results require recalculation with the corrected
metadata. This branch does not validate their previous numerical results,
realm labels, sample counts or figure captions. The iCAMP analysis also requires
recalculation with the corrected cohort and realm assignments; its previous
process percentages are not updated or endorsed here. Its external GTDB
reference phylogeny and null-model assumptions still require explicit reporting.
