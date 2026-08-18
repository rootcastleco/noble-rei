"""Exact similarity checks between finite spherical vertex orbits."""

from __future__ import annotations

from collections import Counter

import sympy as sp

from .algebra import Vec3


def _gram_canon(value: sp.Expr) -> sp.Expr:
    return sp.factor(sp.radsimp(sp.cancel(sp.together(sp.expand(value)))))


def _dot(left: Vec3, right: Vec3) -> sp.Expr:
    return _gram_canon(sum((left[i] * right[i] for i in range(3)), sp.Integer(0)))


def normalized_gram(orbit: tuple[Vec3, ...]) -> tuple[tuple[sp.Expr, ...], ...]:
    if not orbit:
        raise ValueError("orbit is empty")
    radius2 = _dot(orbit[0], orbit[0])
    if radius2 == 0:
        raise ValueError("orbit contains the origin")
    for point in orbit[1:]:
        if _gram_canon(_dot(point, point) - radius2) != 0:
            raise ValueError("orbit is not spherical about the origin")
    return tuple(tuple(_gram_canon(_dot(p, q) / radius2) for q in orbit) for p in orbit)


def _row_signature(row: tuple[sp.Expr, ...]) -> tuple[str, ...]:
    return tuple(sorted(str(value) for value in row))


def are_similar_orbits(left: tuple[Vec3, ...], right: tuple[Vec3, ...]) -> bool:
    if len(left) != len(right):
        return False
    if not left:
        return False
    gl = normalized_gram(left)
    gr = normalized_gram(right)
    n = len(left)

    left_sig = [_row_signature(row) for row in gl]
    right_sig = [_row_signature(row) for row in gr]
    if Counter(left_sig) != Counter(right_sig):
        return False

    candidates: dict[int, tuple[int, ...]] = {
        i: tuple(j for j in range(n) if right_sig[j] == left_sig[i]) for i in range(n)
    }
    mapping: dict[int, int] = {}
    used: set[int] = set()

    def search() -> bool:
        if len(mapping) == n:
            return True
        best_i = -1
        best_options: list[int] | None = None
        for i in range(n):
            if i in mapping:
                continue
            options: list[int] = []
            for j in candidates[i]:
                if j in used:
                    continue
                if all(gl[i][mi] == gr[j][mj] for mi, mj in mapping.items()):
                    options.append(j)
            if not options:
                return False
            if best_options is None or len(options) < len(best_options):
                best_i = i
                best_options = options
                if len(options) == 1:
                    break
        if best_options is None:
            return False
        for j in best_options:
            mapping[best_i] = j
            used.add(j)
            if search():
                return True
            used.remove(j)
            del mapping[best_i]
        return False

    return search()
