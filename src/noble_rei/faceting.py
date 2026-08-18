"""Exact faceting search for a finite point-group orbit."""

from __future__ import annotations

from collections import Counter, deque
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from typing import TypeAlias

import numpy as np
import sympy as sp

from .algebra import Permutation, Vec3, apply_permutation

Face: TypeAlias = tuple[int, ...]
Plane: TypeAlias = frozenset[int]
Edge: TypeAlias = tuple[int, int]
PolyhedronSignature: TypeAlias = tuple[Face, ...]


@dataclass(frozen=True, slots=True)
class Faceting:
    seed_face: Face
    faces: tuple[Face, ...]
    edges: tuple[Edge, ...]


def canonical_face(face: Sequence[int]) -> Face:
    if len(face) < 3:
        raise ValueError("a face must have at least three vertices")
    if len(set(face)) != len(face):
        raise ValueError("a face cannot revisit a vertex")
    cycle = tuple(face)
    variants: list[Face] = []
    for oriented in (cycle, tuple(reversed(cycle))):
        for offset in range(len(oriented)):
            variants.append(oriented[offset:] + oriented[:offset])
    return min(variants)


def face_edges(face: Sequence[int]) -> tuple[Edge, ...]:
    if len(face) < 3:
        raise ValueError("a face must have at least three vertices")
    return tuple(sorted((min(face[i], face[(i + 1) % len(face)]), max(face[i], face[(i + 1) % len(face)])) for i in range(len(face))))


def _vector(point: Vec3) -> sp.ImmutableMatrix:
    return sp.ImmutableMatrix(3, 1, point)


@lru_cache(maxsize=None)
def _numeric_orbit(orbit: tuple[Vec3, ...]) -> np.ndarray:
    return np.asarray([[float(sp.N(value, 20)) for value in point] for point in orbit], dtype=np.float64)


@lru_cache(maxsize=None)
def _plane_through_cached(orbit: tuple[Vec3, ...], i: int, j: int, k: int) -> Plane | None:
    if len({i, j, k}) != 3:
        raise ValueError("plane seed indices must be distinct")
    xyz = _numeric_orbit(orbit)
    p0, p1, p2 = xyz[i], xyz[j], xyz[k]
    normal = np.cross(p1 - p0, p2 - p0)
    normal_norm = float(np.linalg.norm(normal))
    if normal_norm <= 1e-13:
        e0, e1, e2 = (_vector(orbit[x]) for x in (i, j, k))
        exact_normal = (e1 - e0).cross(e2 - e0)
        if all(sp.simplify(component) == 0 for component in exact_normal):
            return None
    scale = max(1.0, float(np.max(np.linalg.norm(xyz - p0, axis=1))))
    tolerance = 2e-11 * max(1.0, normal_norm * scale)
    signed = (xyz - p0) @ normal
    e0 = _vector(orbit[i]); e1 = _vector(orbit[j]); e2 = _vector(orbit[k])
    exact_normal = (e1 - e0).cross(e2 - e0)
    members: set[int] = set()
    for index, approximate in enumerate(signed):
        if abs(float(approximate)) > tolerance:
            continue
        delta = _vector(orbit[index]) - e0
        if sp.simplify(exact_normal.dot(delta)) == 0:
            members.add(index)
    if len(members) < 3:
        raise RuntimeError("a non-collinear plane seed produced fewer than three members")
    return frozenset(members)


def plane_through(orbit: Sequence[Vec3], i: int, j: int, k: int) -> Plane | None:
    orbit_tuple = tuple(orbit)
    a, b, c = sorted((i, j, k))
    return _plane_through_cached(orbit_tuple, a, b, c)


def planes_through_anchor(orbit: Sequence[Vec3], anchor: int = 0) -> tuple[Plane, ...]:
    if not 0 <= anchor < len(orbit):
        raise ValueError("anchor is outside the orbit")
    planes: set[Plane] = set()
    others = [i for i in range(len(orbit)) if i != anchor]
    for j, k in combinations(others, 2):
        plane = plane_through(orbit, anchor, j, k)
        if plane is not None:
            planes.add(plane)
    return tuple(sorted(planes, key=lambda p: (len(p), tuple(sorted(p)))))


def apply_perm_to_plane(plane: Plane, permutation: Permutation) -> Plane:
    return frozenset(permutation[i] for i in plane)


def unique_plane_representatives(planes: Iterable[Plane], group: Sequence[Permutation]) -> tuple[Plane, ...]:
    if not group:
        raise ValueError("permutation group is empty")
    seen_orbits: set[tuple[int, ...]] = set()
    representatives: list[Plane] = []
    for plane in planes:
        images = [tuple(sorted(apply_perm_to_plane(plane, g))) for g in group]
        orbit_key = min(images)
        if orbit_key in seen_orbits:
            continue
        seen_orbits.add(orbit_key)
        representatives.append(plane)
    return tuple(representatives)


def adjacency_graph(plane: Plane, group: Sequence[Permutation]) -> dict[int, frozenset[int]]:
    edges: set[Edge] = set()
    for permutation in group:
        image = apply_perm_to_plane(plane, permutation)
        intersection = plane.intersection(image)
        if len(intersection) == 2:
            a, b = sorted(intersection)
            edges.add((a, b))
    adjacency: dict[int, set[int]] = {vertex: set() for vertex in plane}
    for a, b in edges:
        adjacency[a].add(b); adjacency[b].add(a)
    return {vertex: frozenset(neighbors) for vertex, neighbors in adjacency.items()}


def cycles_through_anchor(adjacency: dict[int, frozenset[int]], anchor: int = 0, min_length: int = 3) -> tuple[Face, ...]:
    if min_length < 3:
        raise ValueError("minimum cycle length must be at least three")
    if anchor not in adjacency or len(adjacency[anchor]) < 2:
        return ()
    found: set[Face] = set(); max_length = len(adjacency)
    def visit(path: list[int], visited: set[int]) -> None:
        if len(path) > max_length:
            return
        current = path[-1]
        for nxt in adjacency[current]:
            if nxt == anchor:
                if len(path) >= min_length:
                    found.add(canonical_face(path))
                continue
            if nxt in visited:
                continue
            visited.add(nxt); path.append(nxt); visit(path, visited); path.pop(); visited.remove(nxt)
    visit([anchor], {anchor})
    return tuple(sorted(found))


def face_orbit(seed_face: Face, group: Sequence[Permutation]) -> tuple[Face, ...]:
    return tuple(sorted({canonical_face(apply_permutation(seed_face, permutation)) for permutation in group}))


def _connected_graph(adjacency: dict[int, set[int]]) -> bool:
    if not adjacency:
        return False
    start = next(iter(adjacency)); seen = {start}; queue: deque[int] = deque([start])
    while queue:
        current = queue.popleft()
        for nxt in adjacency[current]:
            if nxt not in seen:
                seen.add(nxt); queue.append(nxt)
    return len(seen) == len(adjacency)


def _face_adjacency_connected(faces: Sequence[Face]) -> bool:
    edge_sets = [set(face_edges(face)) for face in faces]
    adjacency: dict[int, set[int]] = {i: set() for i in range(len(faces))}
    for i, j in combinations(range(len(faces)), 2):
        if edge_sets[i].intersection(edge_sets[j]):
            adjacency[i].add(j); adjacency[j].add(i)
    return _connected_graph(adjacency)


def _vertex_figures_connected(vertex_count: int, faces: Sequence[Face]) -> bool:
    for vertex in range(vertex_count):
        incident_edges: set[Edge] = set(); adjacency: dict[Edge, set[Edge]] = {}
        for face in faces:
            if vertex not in face:
                continue
            position = face.index(vertex); prev_vertex = face[position - 1]; next_vertex = face[(position + 1) % len(face)]
            e1 = (min(vertex, prev_vertex), max(vertex, prev_vertex)); e2 = (min(vertex, next_vertex), max(vertex, next_vertex))
            incident_edges.update((e1, e2)); adjacency.setdefault(e1, set()).add(e2); adjacency.setdefault(e2, set()).add(e1)
        if not incident_edges:
            return False
        for edge in incident_edges:
            adjacency.setdefault(edge, set())
        if not _connected_graph(adjacency):
            return False
    return True


def validate_faceting(orbit: Sequence[Vec3], seed_face: Face, group: Sequence[Permutation], seed_plane: Plane | None = None) -> Faceting | None:
    if len(orbit) < 4:
        return None
    try:
        canonical_seed = canonical_face(seed_face)
    except ValueError:
        return None
    if any(index < 0 or index >= len(orbit) for index in canonical_seed):
        return None
    faces = face_orbit(canonical_seed, group)
    if not faces:
        return None
    face_vertex_sets = [frozenset(face) for face in faces]
    if len(set(face_vertex_sets)) != len(face_vertex_sets):
        return None
    edge_incidence: Counter[Edge] = Counter()
    for face in faces:
        edge_incidence.update(face_edges(face))
    if not edge_incidence or any(count != 2 for count in edge_incidence.values()):
        return None
    if not _face_adjacency_connected(faces) or not _vertex_figures_connected(len(orbit), faces):
        return None
    used_vertices = {vertex for face in faces for vertex in face}
    if used_vertices != set(range(len(orbit))):
        return None
    face_plane: dict[Face, Plane] = {}
    if seed_plane is not None:
        if not set(canonical_seed).issubset(seed_plane):
            return None
        for permutation in group:
            transformed_face = canonical_face(apply_permutation(canonical_seed, permutation))
            transformed_plane = apply_perm_to_plane(seed_plane, permutation)
            prior = face_plane.setdefault(transformed_face, transformed_plane)
            if prior != transformed_plane:
                return None
    else:
        for face in faces:
            plane = plane_through(orbit, face[0], face[1], face[2])
            if plane is None or not set(face).issubset(plane):
                return None
            face_plane[face] = plane
    edge_sets = [set(face_edges(face)) for face in faces]
    for i, j in combinations(range(len(faces)), 2):
        if edge_sets[i].intersection(edge_sets[j]) and face_plane[faces[i]] == face_plane[faces[j]]:
            return None
    return Faceting(canonical_seed, faces, tuple(sorted(edge_incidence)))


def find_facetings_in_plane(orbit: Sequence[Vec3], group: Sequence[Permutation], plane: Plane, min_cycle_length: int = 3) -> tuple[Faceting, ...]:
    adjacency = adjacency_graph(plane, group)
    cycles = cycles_through_anchor(adjacency, anchor=0, min_length=min_cycle_length)
    valid: dict[PolyhedronSignature, Faceting] = {}
    for cycle in cycles:
        result = validate_faceting(orbit, cycle, group, seed_plane=plane)
        if result is not None:
            valid.setdefault(result.faces, result)
    return tuple(valid.values())


def facet_all(orbit: Sequence[Vec3], group: Sequence[Permutation]) -> tuple[Faceting, ...]:
    planes = planes_through_anchor(orbit)
    representatives = unique_plane_representatives(planes, group)
    found: dict[PolyhedronSignature, Faceting] = {}
    for plane in representatives:
        for faceting in find_facetings_in_plane(orbit, group, plane):
            found.setdefault(faceting.faces, faceting)
    return tuple(found.values())


def canonical_polyhedron_under_group(faces: Sequence[Face], full_vertex_group: Sequence[Permutation]) -> PolyhedronSignature:
    if not faces:
        raise ValueError("polyhedron must contain faces")
    variants: list[PolyhedronSignature] = []
    for permutation in full_vertex_group:
        transformed = tuple(sorted(canonical_face(apply_permutation(face, permutation)) for face in faces))
        variants.append(transformed)
    return min(variants)
