# Relationship to MAGENTA

MAGENTA is the broader project and historical code base. This repository is
the manuscript-specific reproducibility companion for the global biogeography
analysis.

The repositories are linked but serve different purposes:

| Repository | Scope |
|---|---|
| [`fjbalvino/magenta`](https://github.com/fjbalvino/magenta) | Broader MAGENTA project, exploration and related analyses |
| `fjbalvino/mangrove-sediment-biogeography` | Frozen, ordered and audited workflow supporting the global biogeography manuscript |

To preserve reproducibility, this repository does not source scripts directly
from the changing MAGENTA default branch. Shared scripts are copied at a known
version, renamed by scientific objective and recorded in the provenance tables.
Future fixes should be propagated explicitly through reviewed commits rather
than through hidden cross-repository dependencies.

