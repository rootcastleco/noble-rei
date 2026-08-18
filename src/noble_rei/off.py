"""Minimal deterministic OFF exporter."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import sympy as sp

from .algebra import Vec3
from .faceting import Face


def write_off(
    path: Path,
    vertices: Sequence[Vec3],
    faces: Sequence[Face],
    *,
    substitutions: dict[sp.Symbol, float] | None = None,
    precision: int = 16,
) -> None:
    if not vertices:
        raise ValueError("cannot export an empty vertex set")
    if not faces:
        raise ValueError("cannot export a polyhedron without faces")
    if precision < 6:
        raise ValueError("precision must be at least 6 digits")
    substitutions = substitutions or {}
    path.parent.mkdir(parents=True, exist_ok=True)

    lines = ["OFF", f"{len(vertices)} {len(faces)} 0"]
    for point in vertices:
        if len(point) != 3:
            raise ValueError("OFF exporter expects 3D vertices")
        values = [float(sp.N(coord.subs(substitutions), precision + 4)) for coord in point]
        if not all(sp.Float(value).is_finite for value in values):
            raise ValueError("non-finite coordinate generated during OFF export")
        lines.append(" ".join(f"{value:.{precision}g}" for value in values))

    vertex_count = len(vertices)
    for face in faces:
        if len(face) < 3:
            raise ValueError("OFF face has fewer than three vertices")
        if any(index < 0 or index >= vertex_count for index in face):
            raise ValueError("OFF face references an invalid vertex")
        lines.append(f"{len(face)} " + " ".join(str(index) for index in face))

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
