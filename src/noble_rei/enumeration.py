"""Enumeration orchestration."""

from __future__ import annotations

from dataclasses import dataclass
import multiprocessing as mp
from multiprocessing.connection import Connection

from .algebra import permutation_action
from .faceting import Faceting, PolyhedronSignature, canonical_polyhedron_under_group, facet_all
from .orbits import EXPECTED_ZERO_D_COUNTS, ZERO_D_ORBITS, OrbitSpec
from .point_groups import matrix_group


@dataclass(frozen=True, slots=True)
class OrbitEnumeration:
    orbit_name: str
    vertex_count: int
    symmetry_counts: tuple[tuple[str, int], ...]
    unique_polyhedra: tuple[PolyhedronSignature, ...]

    @property
    def count(self) -> int:
        return len(self.unique_polyhedra)


def enumerate_zero_d_orbit(spec: OrbitSpec) -> OrbitEnumeration:
    vertices = spec.vertices()
    full_action = permutation_action(vertices, matrix_group(spec.full_vertex_group))

    unique: dict[PolyhedronSignature, Faceting] = {}
    symmetry_counts: list[tuple[str, int]] = []
    for group_name in spec.enumeration_groups:
        action = permutation_action(vertices, matrix_group(group_name))
        facetings = facet_all(vertices, action)
        symmetry_counts.append((group_name, len(facetings)))
        for faceting in facetings:
            key = canonical_polyhedron_under_group(faceting.faces, full_action)
            unique.setdefault(key, faceting)

    return OrbitEnumeration(
        orbit_name=spec.name,
        vertex_count=len(vertices),
        symmetry_counts=tuple(symmetry_counts),
        unique_polyhedra=tuple(sorted(unique)),
    )


def _zero_d_worker(spec: OrbitSpec, sender: Connection) -> None:
    try:
        sender.send((True, enumerate_zero_d_orbit(spec)))
    except BaseException as exc:
        sender.send((False, (type(exc).__name__, str(exc))))
    finally:
        sender.close()


def enumerate_zero_d(*, isolated: bool = True) -> tuple[OrbitEnumeration, ...]:
    if not isolated:
        return tuple(enumerate_zero_d_orbit(spec) for spec in ZERO_D_ORBITS)

    context = mp.get_context("spawn")
    results: list[OrbitEnumeration] = []
    for spec in ZERO_D_ORBITS:
        receiver, sender = context.Pipe(duplex=False)
        process = context.Process(target=_zero_d_worker, args=(spec, sender), daemon=False)
        process.start()
        sender.close()
        try:
            ok, payload = receiver.recv()
        except EOFError as exc:
            process.join()
            raise RuntimeError(
                f"worker for {spec.name} exited without returning a result; exit={process.exitcode}"
            ) from exc
        finally:
            receiver.close()
        process.join()
        if process.exitcode != 0:
            raise RuntimeError(f"worker for {spec.name} failed with exit code {process.exitcode}")
        if not ok:
            error_type, message = payload
            raise RuntimeError(f"worker for {spec.name} failed: {error_type}: {message}")
        if not isinstance(payload, OrbitEnumeration):
            raise RuntimeError(f"worker for {spec.name} returned an invalid payload")
        results.append(payload)
    return tuple(results)


def verify_zero_d(results: tuple[OrbitEnumeration, ...]) -> None:
    actual = {result.orbit_name: result.count for result in results}
    if actual != EXPECTED_ZERO_D_COUNTS:
        raise AssertionError(f"zero-dimensional enumeration mismatch: {actual}")
    if sum(actual.values()) != 20:
        raise AssertionError(f"expected 20 zero-dimensional polyhedra, got {sum(actual.values())}")
