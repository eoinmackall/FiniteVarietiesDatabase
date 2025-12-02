This directory contains scripts for collecting polynomial representatives for hypersurfaces in projective space under projective equivalence. The corresponding database will be available by a GitHub release in a future version.

## hypersurfaces.db Database

The following tables show what equivalence classes of hypersurfaces are included in the database.

---

**GF(2)**
| dim/deg| 2 | 3 | 4 | 5 | 6 | 7 |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | X | X | X | X | X | !! |
| 2 | X | X | X |  |  |  |
| 3 |  |  |  |  |  |  |
| 4 |  |  |  |  |  |  |

---

**GF(3)**
| dim/deg| 2 | 3 | 4 | 5 | 6 | 7 |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | X | X | X | X |  |  |
| 2 |  |  |  |  |  |  |
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
- An `X` in any box above means that the script to collect hypersurfaces, of the given dimension and degree over the specified field, has completed. 
- Two exclamation points `!!` means that the script is running. 
- An empty box indicates that no script has been ran for this equivalence class.

