"""Coxeter point groups and orbit definitions used by the noble-polyhedra search."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Final

import sympy as sp

from .algebra import Vec3, matrix_group_closure, orbit_from_group

SQRT2: Final = sp.sqrt(2)
PHI: Final = (sp.Integer(1) + sp.sqrt(5)) / 2
HALF: Final = sp.Rational(1, 2)

R1_TET = sp.ImmutableMatrix([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
R2_TET = sp.ImmutableMatrix([[0, -1, 0], [-1, 0, 0], [0, 0, 1]])
R3_TET = sp.ImmutableMatrix([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
R1_OCT = sp.ImmutableMatrix([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
R2_OCT = R1_TET
R3_OCT = R3_TET
R1_ICO = R1_OCT
R2_ICO = sp.ImmutableMatrix([[1, 0, 0], [0, -1, 0], [0, 0, 1]])
R3_ICO = HALF * sp.ImmutableMatrix([[1 - PHI, -PHI, 1], [-PHI, 1, PHI - 1], [1, PHI - 1, PHI]])

BASIS_TET: Final[tuple[Vec3, Vec3, Vec3]] = (
    (-SQRT2 / 2, SQRT2 / 2, SQRT2 / 2),
    (SQRT2 / 2, SQRT2 / 2, SQRT2 / 2),
    (sp.Integer(0), sp.Integer(0), SQRT2),
)
BASIS_OCT: Final[tuple[Vec3, Vec3, Vec3]] = (
    (sp.Integer(1), sp.Integer(1), sp.Integer(1)),
    (sp.Integer(0), SQRT2, SQRT2),
    (sp.Integer(0), sp.Integer(0), SQRT2),
)
BASIS_ICO: Final[tuple[Vec3, Vec3, Vec3]] = (
    (sp.Integer(1), sp.Integer(0), PHI + 1),
    (sp.Integer(0), sp.Integer(1), PHI),
    (sp.Integer(0), sp.Integer(0), 2 * PHI),
)


@dataclass(frozen=True, slots=True)
class GroupDefinition:
    name: str
    generators: tuple[sp.ImmutableMatrix, ...]


GROUP_DEFINITIONS: Final[dict[str, GroupDefinition]] = {
    "*332": GroupDefinition("*332", (R1_TET, R2_TET, R3_TET)),
    "332": GroupDefinition("332", (R1_TET * R2_TET, R1_TET * R3_TET, R2_TET * R3_TET)),
    "*432": GroupDefinition("*432", (R1_OCT, R2_OCT, R3_OCT)),
    "432": GroupDefinition("432", (R1_OCT * R2_OCT, R1_OCT * R3_OCT, R2_OCT * R3_OCT)),
    "*532": GroupDefinition("*532", (R1_ICO, R2_ICO, R3_ICO)),
    "532": GroupDefinition("532", (R1_ICO * R2_ICO, R1_ICO * R3_ICO, R2_ICO * R3_ICO)),
    "3*2": GroupDefinition("3*2", (R1_OCT, R2_OCT * R1_OCT * R2_OCT, R3_OCT * R2_OCT)),
}


@lru_cache(maxsize=None)
def matrix_group(name: str) -> tuple[sp.ImmutableMatrix, ...]:
    try:
        definition = GROUP_DEFINITIONS[name]
    except KeyError as exc:
        raise ValueError(f"unknown point group: {name}") from exc
    return matrix_group_closure(definition.generators)


def _linear_combination(basis: tuple[Vec3, Vec3, Vec3], a: sp.Expr, b: sp.Expr, c: sp.Expr) -> Vec3:
    values = tuple(sp.expand(a * basis[0][i] + b * basis[1][i] + c * basis[2][i]) for i in range(3))
    return values  # type: ignore[return-value]


def orbit(group_name: str, a: sp.Expr, b: sp.Expr, c: sp.Expr) -> tuple[Vec3, ...]:
    if group_name in {"*332", "332"}:
        basis = BASIS_TET
    elif group_name in {"*432", "432", "3*2"}:
        basis = BASIS_OCT
    elif group_name in {"*532", "532"}:
        basis = BASIS_ICO
    else:
        raise ValueError(f"unsupported orbit group: {group_name}")
    initial = _linear_combination(basis, sp.sympify(a), sp.sympify(b), sp.sympify(c))
    return orbit_from_group(initial, matrix_group(group_name))
