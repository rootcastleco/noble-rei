# Verification record

Date: 2026-08-18

Environment used for the packaged build:

- Python 3.13.5
- SymPy 1.14.0
- NumPy 2.3.5
- `wolframscript`: not installed in this environment

## Point-group closure

Generated orders:

| Group | Order |
|---|---:|
| 332 | 12 |
| *332 | 24 |
| 3*2 | 24 |
| 432 | 24 |
| *432 | 48 |
| 532 | 60 |
| *532 | 120 |

## 0-DOF independent enumeration

All seven zero-degree-of-freedom orbit types were executed through the final CLI.

| Orbit | Generated | Published invariant |
|---|---:|---:|
| T | 1 | 1 |
| O | 1 | 1 |
| CO | 0 | 0 |
| C | 1 | 1 |
| I | 4 | 4 |
| ID | 6 | 6 |
| D | 7 | 7 |
| **Total** | **20** | **20** |

Notable per-symmetry observations from the run:

- `ID`: `*532 -> 6`, `532 -> 3`, canonical union `6`.
- `D`: `*532 -> 6`, `532 -> 8`, canonical union `7`.

These differences are expected: a polyhedron may be found under multiple symmetry
descriptions, so the final count is after full-group canonicalization.

## 1-DOF independent enumeration — verified subset

The final package was executed for the bounded SymPy subset:

| Orbit | Generated | Generic class facetings | Published invariant |
|---|---:|---:|---:|
| tT | 0 | 0 | 0 |
| rT | 0 | 0 | 0 |
| rP | 0 | 0 | 0 |
| tO | 1 | 0 | 1 |
| tC | 1 | 0 | 1 |
| rC | 1 | 0 | 1 |

`tO` is generated at the positive real root approximately
`a = 2.1478990357047874` of `a^3 - a^2 - 2a - 1`.

A regression was found during development: a naive critical-plane implementation
counted two `rP` candidates at `a = phi`. Exact normalized-Gram comparison showed that
this vertex orbit is a scaled/rotated copy of the 0-DOF icosahedral `I` orbit. The
specialization filter removes it, restoring the published `rP = 0` result.

## Explicitly not verified independently in this environment

- `tI = 17`
- `tD = 6`
- `rD = 19`
- `sT/gT/gP/sC/gC/sD/gD` 2D enumeration

The published reference distribution for these cases is present only as a post-generation
assertion table. Full reproduction is delegated to `upstream-run`, which requires the
published repository and Wolfram Language. Because `wolframscript` is absent here, that
full CAS pipeline was not executed for this build.
