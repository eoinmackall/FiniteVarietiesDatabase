This directory contains scripts for collecting polynomial representatives for hypersurfaces in projective space under projective equivalence. The corresponding database will be available by a GitHub release in a future version.

## hypersurfaces.db Database

The following tables show what equivalence classes of hypersurfaces are included in the database. Rows are labeled by the dimension of the hypersurface and columns are labeled by degree.

---

**GF(2)**
| dim/deg| 2 | 3 | 4 | 5 | 6 | 7 |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | X | X | X | X | X | K64 |
| 2 | X | X | X | K64 |  |  |
| 3 | X | X | K64 |  |  |  |
| 4 | X | K64 |  |  |  |  |
| 5 |   |   |   |   |   |   |
---

**GF(3)**
| dim/deg| 2 | 3 | 4 | 5 | 6 | 7 |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | X | X | X | X |  |  |
| 2 | X  |  |  |  |  |  |
| 3 |  |  |  |  |  |  |
| 4 |  |  |  |  |  |  |

---

**GF(5)**
| dim/deg| 2 | 3 | 4 | 5 | 6 | 7 |
| --- | --- | --- | --- | --- | --- | --- |
| 1 |  |  |  |  |  |  |
| 2 |  |  |  |  |  |  |
| 3 |  |  |  |  |  |  |
| 4 |  |  |  |  |  |  |

---
- An `X` in any box above means that the script to collect representatives for equivalence classes of hypersurfaces, of the given dimension and degree over the specified field, has completed. 
- Two exclamation points `!!` means that their is a data-collection script running in this dimension and degree. 
- A `K##` entry indicates that a data-collection script was run, and killed, on a machine with `##` G of RAM and no script has been run on a machine with more available memory.
- An empty box indicates that no script has been ran for this equivalence class.

