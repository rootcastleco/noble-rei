"""Orbit-type catalog from the non-prismatic cases in the paper."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

import sympy as sp

from .algebra import Vec3
from .point_groups import orbit
from .critical import A


@dataclass(frozen=True, slots=True)
class OrbitSpec:
    name: str
    generator_group: str
    parameters: tuple[sp.Expr, sp.Expr, sp.Expr]
    enumeration_groups: tuple[str, ...]
    full_vertex_group: str

    def vertices(self) -> tuple[Vec3, ...]:
        return orbit(self.generator_group, *self.parameters)


ZERO_D_ORBITS: Final[tuple[OrbitSpec, ...]] = (
    OrbitSpec("T", "*332", (sp.Integer(1), sp.Integer(0), sp.Integer(0)), ("*332", "332"), "*332"),
    OrbitSpec("O", "*432", (sp.Integer(0), sp.Integer(0), sp.Integer(1)), ("*332", "332", "3*2", "*432", "432"), "*432"),
    OrbitSpec("CO", "*432", (sp.Integer(0), sp.Integer(1), sp.Integer(0)), ("3*2", "*432", "432"), "*432"),
    OrbitSpec("C", "*432", (sp.Integer(1), sp.Integer(0), sp.Integer(0)), ("3*2", "*432", "432"), "*432"),
    OrbitSpec("I", "*532", (sp.Integer(0), sp.Integer(1), sp.Integer(0)), ("*532", "532"), "*532"),
    OrbitSpec("ID", "*532", (sp.Integer(0), sp.Integer(0), sp.Integer(1)), ("*532", "532"), "*532"),
    OrbitSpec("D", "*532", (sp.Integer(1), sp.Integer(0), sp.Integer(0)), ("*532", "532"), "*532"),
)

EXPECTED_ZERO_D_COUNTS: Final[dict[str, int]] = {"T": 1, "O": 1, "CO": 0, "C": 1, "I": 4, "ID": 6, "D": 7}

ONE_D_ORBITS: Final[tuple[OrbitSpec, ...]] = (
    OrbitSpec("tT", "*332", (sp.Integer(0), A, sp.Integer(1)), ("*332", "332"), "*332"),
    OrbitSpec("rT", "*332", (A, sp.Integer(1), sp.Integer(0)), ("*332", "332"), "*332"),
    OrbitSpec("rP", "3*2", (sp.Integer(0), A, sp.Integer(1)), ("3*2",), "3*2"),
    OrbitSpec("tO", "*432", (sp.Integer(0), A, sp.Integer(1)), ("*432", "432"), "*432"),
    OrbitSpec("tC", "*432", (A, sp.Integer(1), sp.Integer(0)), ("*432", "432", "3*2"), "*432"),
    OrbitSpec("rC", "*432", (A, sp.Integer(0), sp.Integer(1)), ("*432", "432", "3*2"), "*432"),
    OrbitSpec("tI", "*532", (sp.Integer(0), A, sp.Integer(1)), ("*532", "532"), "*532"),
    OrbitSpec("tD", "*532", (A, sp.Integer(0), sp.Integer(1)), ("*532", "532"), "*532"),
    OrbitSpec("rD", "*532", (A, sp.Integer(1), sp.Integer(0)), ("*532", "532"), "*532"),
)

EXPECTED_ONE_D_COUNTS: Final[dict[str, int]] = {"tT": 0, "rT": 0, "rP": 0, "tO": 1, "tC": 1, "rC": 1, "tI": 17, "tD": 6, "rD": 19}
