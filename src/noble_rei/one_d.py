"""One-parameter critical-class enumeration."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

from .algebra import permutation_action
from .critical import A, CriticalFactor1D, critical_factors_1d, merge_planes, volume_configuration
from .faceting import Faceting, PolyhedronSignature, canonical_polyhedron_under_group, find_facetings_in_plane, unique_plane_representatives, validate_faceting
from .orbits import OrbitSpec, ZERO_D_ORBITS
from .point_groups import matrix_group
from .geometry import are_similar_orbits


@dataclass(frozen=True, slots=True)
class OneDRealization:
    polynomial: str
    root: float
    faces: PolyhedronSignature
    symmetry_sources: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class OneDEnumeration:
    orbit_name: str
    vertex_count: int
    critical_factor_count: int
    critical_root_count: int
    generic_faceting_count: int
    realizations: tuple[OneDRealization, ...]

    @property
    def count(self) -> int:
        return len(self.realizations)


def _generic_facetings(symbolic_orbit: tuple, configuration: dict, action: tuple[tuple[int, ...], ...]) -> tuple[Faceting, ...]:
    found: dict[PolyhedronSignature, Faceting] = {}
    vertex_count = len(symbolic_orbit)
    for i, j in combinations(range(1, vertex_count), 2):
        face = (0, i, j)
        plane = frozenset(face)
        result = validate_faceting(symbolic_orbit, face, action, seed_plane=plane)
        if result is not None:
            found.setdefault(result.faces, result)
    shared_triples = configuration.get(0, ())
    if shared_triples:
        shared_planes = merge_planes(list(shared_triples))
        for plane in unique_plane_representatives(shared_planes, action):
            for result in find_facetings_in_plane(symbolic_orbit, action, plane, min_cycle_length=4):
                found.setdefault(result.faces, result)
    return tuple(found.values())


def _critical_candidates_for_group(symbolic_orbit: tuple, action: tuple[tuple[int, ...], ...], factors: tuple[CriticalFactor1D, ...]) -> dict[tuple[str, int, PolyhedronSignature], Faceting]:
    candidates: dict[tuple[str, int, PolyhedronSignature], Faceting] = {}
    for factor in factors:
        representatives = unique_plane_representatives(factor.planes, action)
        factor_facetings: dict[PolyhedronSignature, Faceting] = {}
        for plane in representatives:
            for faceting in find_facetings_in_plane(symbolic_orbit, action, plane, min_cycle_length=4):
                factor_facetings.setdefault(faceting.faces, faceting)
        for root_index, _root in enumerate(factor.positive_roots):
            for signature, faceting in factor_facetings.items():
                candidates[(str(factor.polynomial), root_index, signature)] = faceting
    return candidates


def _is_zero_d_specialization(spec: OrbitSpec, polynomial: str, root_index: int) -> bool:
    sympy = __import__("sympy")
    expr = sympy.sympify(polynomial, locals={"a": A})
    poly = sympy.Poly(expr, A, extension=True)
    if poly.degree() != 1 or root_index != 0:
        return False
    exact_root = sympy.solve(poly.as_expr(), A)[0]
    specialized = tuple(tuple(sympy.simplify(value.subs(A, exact_root)) for value in point) for point in spec.vertices())
    for zero_spec in ZERO_D_ORBITS:
        if len(zero_spec.vertices()) != len(specialized):
            continue
        if are_similar_orbits(specialized, zero_spec.vertices()):
            return True
    return False


def enumerate_one_d_orbit(spec: OrbitSpec) -> OneDEnumeration:
    symbolic_orbit = spec.vertices()
    full_action = permutation_action(symbolic_orbit, matrix_group(spec.full_vertex_group))
    configuration = volume_configuration(symbolic_orbit, A, symmetry_action=full_action)
    factors = critical_factors_1d(configuration, A)
    generic: dict[PolyhedronSignature, Faceting] = {}
    realization_sources: dict[tuple[str, int, PolyhedronSignature], set[str]] = {}
    roots_by_factor = {str(f.polynomial): f.positive_roots for f in factors}
    for group_name in spec.enumeration_groups:
        action = permutation_action(symbolic_orbit, matrix_group(group_name))
        for faceting in _generic_facetings(symbolic_orbit, configuration, action):
            key = canonical_polyhedron_under_group(faceting.faces, full_action)
            generic.setdefault(key, faceting)
        group_candidates = _critical_candidates_for_group(symbolic_orbit, action, factors)
        for (polynomial, root_index, signature), faceting in group_candidates.items():
            canonical = canonical_polyhedron_under_group(signature, full_action)
            key = (polynomial, root_index, canonical)
            realization_sources.setdefault(key, set()).add(group_name)
    realizations: list[OneDRealization] = []
    for (polynomial, root_index, faces), sources in realization_sources.items():
        roots = roots_by_factor[polynomial]
        if not 0 <= root_index < len(roots):
            raise RuntimeError("critical-root index invariant violated")
        if _is_zero_d_specialization(spec, polynomial, root_index):
            continue
        realizations.append(OneDRealization(polynomial=polynomial, root=roots[root_index], faces=faces, symmetry_sources=tuple(sorted(sources))))
    realizations.sort(key=lambda r: (r.polynomial, r.root, r.faces))
    return OneDEnumeration(orbit_name=spec.name, vertex_count=len(symbolic_orbit), critical_factor_count=len(factors), critical_root_count=sum(len(f.positive_roots) for f in factors), generic_faceting_count=len(generic), realizations=tuple(realizations))
