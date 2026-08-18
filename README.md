# noble-rei

Exact finite-group/faceting kernel plus a reproducibility gateway for Connor Hill's
2026 classification of noble polyhedra.

Primary paper: **Connor Hill, “The complete set of noble polyhedra”**, arXiv:2607.28711.
Published reference code: **Plasmath/noble-tools-revised** (GPL-3.0).

## Scientific status

This package deliberately separates **independently reproduced results** from the
**canonical upstream proof pipeline**.

| Scope | Independent backend | Verified here | Full-proof route |
|---|---:|---:|---|
| 0-DOF: T/O/CO/C/I/ID/D | implemented | **20/20** | independent Python |
| 1-DOF: tT/rT/rP/tO/tC/rC | implemented | **0/0/0/1/1/1** | independent Python |
| 1-DOF: tI/tD/rD | algebraic primitives implemented | **not claimed** | upstream Wolfram pipeline |
| 2-DOF: sT/gT/gP/sC/gC/sD/gD | resultant primitive only | **not claimed** | upstream Wolfram pipeline |
| Published nonprismatic total | reference invariant | **146** in paper/upstream | upstream full run |

**Do not interpret the independent backend alone as a new proof of all 146 cases.**
The package refuses to run the unverified heavy icosahedral 1D path as if it were a
validated result.

## What is independently implemented

- Exact Coxeter point groups `332`, `*332`, `3*2`, `432`, `*432`, `532`, `*532`.
- Group closure until no new element exists; no fixed “N iterations” assumption.
- Exact point orbits and permutation actions.
- Coplanarity/volume configuration polynomials.
- A fast exact `Q(sqrt(2), sqrt(5))` coefficient path for 1D configurations.
- Adjacency graph construction and simple face-cycle enumeration.
- Abstract-polyhedron checks: two faces per edge, connectivity, vertex figures,
  injective face sets, and no adjacent coplanar faces.
- Exact normalized-Gram orbit similarity, used to detect lower-dimensional
  specializations such as `rP(phi)` collapsing to the icosahedral `I` orbit.
- Deterministic OFF export.
- A gateway that executes Hill's published Python/Wolfram scripts without modifying
  the independent core.

## Install

Python 3.12+ is required.

```bash
python -m pip install -e .
```

If the environment has no package-index access but already contains NumPy/SymPy:

```bash
python -m pip install -e . --no-build-isolation
```

## Commands

```bash
noble-rei doctor
noble-rei groups
noble-rei reference
noble-rei enumerate-0d-orbit ID --output out/id
noble-rei enumerate-0d-orbit D --output out/d
noble-rei enumerate-1d-orbit tO --output out/to
noble-rei enumerate-1d-orbit rP
```

For the complete canonical enumeration, clone the published upstream repository,
install Wolfram Language / `wolframscript`, then run:

```bash
noble-rei upstream-run /path/to/noble-tools-revised --phase all
```

## Tests

```bash
pytest -q
```

See `VERIFICATION.md` for the verified scope and runtime checks.

## License

GPL-3.0-or-later. See `LICENSE` and `NOTICE.md`.
