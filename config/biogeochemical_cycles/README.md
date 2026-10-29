# Frozen biogeochemical analysis definitions

These small text files configure an analysis. They contain no sample counts,
sequences, reference database or numerical results.

| File | Role |
| --- | --- |
| `reglas_operativas.tsv` | 97 stable rule IDs, parent process, logical expression, evidence type, limitation, source and representative (`primary`) flag. |
| `procesos_curados.tsv` | 65 parent categories, candidate genes and interpretation limits. Candidate genes are broader than the operational KO rules. |
| `procesos_no_resueltos.json` | Nine categories requiring evidence beyond the available KO matrix. |
| `modulos_KEGG_definiciones.tsv` | 32 selected reference definitions, frozen expressions and source URLs. |
| `module_dependencies.json` | M00597 and M00598 expressions used to expand nested references. |
| `provenance.json` | Configuration version, units, thresholds and scope. |

Evidence types retain their original Spanish keys: `conjunto` (component set),
`ambiguo` (ambiguous assignment/direction), `contexto` (contextual component), and
`modulo` (selected module definition). A category has exactly one representative
rule unless explicitly unresolved; all alternative rules remain in the long tables.

Sources are resolved in [the references](../../docs/BIOGEOCHEMICAL_REFERENCES.md).
Broad category sources do not mean that every gene listed is diagnostic or that
the cited tools/databases were executed. CCycDB is retained as a preprint source
in this research snapshot. The packaged configuration date is not the date or
version of the historical matrix annotation database.

To revise a rule, document the scientific reason and source, preserve its meaning
or assign a new ID, update the configuration version and tests, and compare outputs
with the preceding configuration. Do not silently fetch current KEGG definitions
inside a reproduction run. External sources retain their own terms; these selected
analysis expressions do not constitute a redistribution of those databases.
