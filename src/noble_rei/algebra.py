"""Exact algebra helpers for finite 3D point groups."""

from __future__ import annotations

from collections import deque
from collections.abc import Iterable, Sequence
from typing import TypeAlias

import sympy as sp

Expr: TypeAlias = sp.Expr
Vec3: TypeAlias = tuple[Expr, Expr, Expr]
Mat3: TypeAlias = tuple[Expr, Expr, Expr, Expr, Expr, Expr, Expr, Expr, Expr]
Permutation: TypeAlias = tuple[int, ...]


def canon_expr(value: sp.Expr) -> sp.Expr:
    return sp.factor(sp.cancel(sp.expand(value)))


def vec_key(vector: Sequence[sp.Expr]) -> Vec3:
    if len(vector) != 3:
        raise ValueError("expected a 3-vector")
    return tuple(canon_expr(sp.sympify(x)) for x in vector)  # type: ignore[return-value]


def matrix_key(matrix: sp.MatrixBase) -> Mat3:
    if matrix.shape != (3, 3):
        raise ValueError("expected a 3x3 matrix")
    return tuple(canon_expr(sp.sympify(matrix[i, j])) for i in range(3) for j in range(3))  # type: ignore[return-value]


def matrix_from_key(key: Mat3) -> sp.ImmutableMatrix:
    if len(key) != 9:
        raise ValueError("invalid Mat3 key")
    return sp.ImmutableMatrix(3, 3, key)


def compose(left: Permutation, right: Permutation) -> Permutation:
    if len(left) != len(right):
        raise ValueError("permutations must have equal degree")
    return tuple(left[right[i]] for i in range(len(right)))


def matrix_group_closure(generators: Sequence[sp.MatrixBase]) -> tuple[sp.ImmutableMatrix, ...]:
    if not generators:
        raise ValueError("at least one generator is required")
    if any(g.shape != (3, 3) for g in generators):
        raise ValueError("all generators must be 3x3")

    identity = sp.ImmutableMatrix(sp.eye(3))
    seen: dict[Mat3, sp.ImmutableMatrix] = {matrix_key(identity): identity}
    queue: deque[sp.ImmutableMatrix] = deque([identity])
    canonical_generators = tuple(matrix_from_key(matrix_key(g)) for g in generators)

    while queue:
        current = queue.popleft()
        for generator in canonical_generators:
            candidate = matrix_from_key(matrix_key(generator * current))
            key = matrix_key(candidate)
            if key in seen:
                continue
            seen[key] = candidate
            queue.append(candidate)

    return tuple(seen.values())


def orbit_from_group(initial: Sequence[sp.Expr], group: Sequence[sp.MatrixBase]) -> tuple[Vec3, ...]:
    if len(initial) != 3:
        raise ValueError("initial point must have three coordinates")
    if not group:
        raise ValueError("group must be non-empty")

    p = sp.ImmutableMatrix(3, 1, [sp.sympify(x) for x in initial])
    initial_key = vec_key(tuple(p))
    points: dict[Vec3, Vec3] = {}
    for matrix in group:
        transformed = vec_key(tuple(matrix * p))
        points.setdefault(transformed, transformed)

    if initial_key not in points:
        raise RuntimeError("identity image missing from generated orbit")

    rest = sorted((key for key in points if key != initial_key), key=lambda v: tuple(map(str, v)))
    return (initial_key, *rest)


def permutation_action(orbit: Sequence[Vec3], group: Sequence[sp.MatrixBase]) -> tuple[Permutation, ...]:
    if not orbit:
        raise ValueError("orbit must be non-empty")
    index = {point: i for i, point in enumerate(orbit)}
    if len(index) != len(orbit):
        raise ValueError("orbit contains duplicate points")

    permutations: set[Permutation] = set()
    for matrix in group:
        image: list[int] = []
        for point in orbit:
            vector = sp.ImmutableMatrix(3, 1, point)
            key = vec_key(tuple(matrix * vector))
            try:
                image.append(index[key])
            except KeyError as exc:
                raise ValueError("group does not preserve supplied orbit") from exc
        permutations.add(tuple(image))

    identity = tuple(range(len(orbit)))
    if identity not in permutations:
        raise RuntimeError("permutation action is missing identity")
    return tuple(sorted(permutations))


def apply_permutation(items: Iterable[int], permutation: Permutation) -> tuple[int, ...]:
    return tuple(permutation[i] for i in items)
