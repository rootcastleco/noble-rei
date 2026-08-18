"""Command-line interface."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sys

from .critical import A
from .enumeration import enumerate_zero_d_orbit
from .off import write_off
from .one_d import enumerate_one_d_orbit
from .orbits import ONE_D_ORBITS, ZERO_D_ORBITS
from .point_groups import GROUP_DEFINITIONS, matrix_group
from .reference import PUBLISHED_COUNTS, PUBLISHED_TOTAL, verify_reference_table
from .upstream import run_upstream


def _orbit_spec(name: str, catalog: tuple):
    try:
        return next(spec for spec in catalog if spec.name == name)
    except StopIteration as exc:
        raise SystemExit(f"unknown orbit type: {name}") from exc


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="noble-rei")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("groups", help="print generated finite point-group orders")
    sub.add_parser("reference", help="print the published 146-result distribution")
    sub.add_parser("doctor", help="report local runtime/CAS capabilities")
    p0 = sub.add_parser("enumerate-0d-orbit", help="enumerate one 0-DOF orbit type")
    p0.add_argument("orbit", choices=[spec.name for spec in ZERO_D_ORBITS])
    p0.add_argument("--output", type=Path)
    p1 = sub.add_parser("enumerate-1d-orbit", help="enumerate one 1-DOF orbit type")
    p1.add_argument("orbit", choices=[spec.name for spec in ONE_D_ORBITS])
    p1.add_argument("--output", type=Path)
    up = sub.add_parser("upstream-run", help="run the published upstream GPL pipeline")
    up.add_argument("repo", type=Path, help="path to Plasmath/noble-tools-revised checkout")
    up.add_argument("--phase", choices=("0d", "1d", "2d", "all"), default="all")
    return parser


def _emit_zero_d(name: str, output: Path | None) -> None:
    spec = _orbit_spec(name, ZERO_D_ORBITS)
    result = enumerate_zero_d_orbit(spec)
    expected = PUBLISHED_COUNTS[name]
    if result.count != expected:
        raise RuntimeError(f"{name}: generated {result.count}, expected {expected}")
    if output is not None:
        vertices = spec.vertices()
        for index, faces in enumerate(result.unique_polyhedra, start=1):
            write_off(output / f"{name}-{index}.off", vertices, faces)
    print(json.dumps({"orbit": name, "vertices": result.vertex_count, "count": result.count, "expected": expected, "symmetry_counts": dict(result.symmetry_counts)}, indent=2))


def _emit_one_d(name: str, output: Path | None) -> None:
    if name in {"tI", "tD", "rD"}:
        raise SystemExit(f"{name} is intentionally disabled in the independent SymPy backend: use upstream-run with Wolfram Language for the canonical proof pipeline.")
    spec = _orbit_spec(name, ONE_D_ORBITS)
    result = enumerate_one_d_orbit(spec)
    expected = PUBLISHED_COUNTS[name]
    if result.count != expected:
        raise RuntimeError(f"{name}: generated {result.count}, expected {expected}")
    if result.generic_faceting_count != 0:
        raise RuntimeError(f"{name}: unexpected generic-class faceting(s)")
    if output is not None:
        symbolic_vertices = spec.vertices()
        for index, realization in enumerate(result.realizations, start=1):
            write_off(output / f"{name}-{index}.off", symbolic_vertices, realization.faces, substitutions={A: realization.root})
    print(json.dumps({"orbit": name, "vertices": result.vertex_count, "count": result.count, "expected": expected, "critical_factors": result.critical_factor_count, "critical_roots": result.critical_root_count, "generic_facetings": result.generic_faceting_count, "roots": [r.root for r in result.realizations]}, indent=2))


def main() -> None:
    args = _build_parser().parse_args()
    if args.command == "groups":
        print(json.dumps({name: len(matrix_group(name)) for name in GROUP_DEFINITIONS}, indent=2, sort_keys=True)); return
    if args.command == "reference":
        verify_reference_table(); print(json.dumps({"total": PUBLISHED_TOTAL, "counts": PUBLISHED_COUNTS}, indent=2, sort_keys=True)); return
    if args.command == "doctor":
        print(json.dumps({"python": sys.version.split()[0], "wolframscript": shutil.which("wolframscript"), "full_upstream_1d_2d_ready": shutil.which("wolframscript") is not None}, indent=2)); return
    if args.command == "enumerate-0d-orbit":
        _emit_zero_d(args.orbit, args.output); return
    if args.command == "enumerate-1d-orbit":
        _emit_one_d(args.orbit, args.output); return
    if args.command == "upstream-run":
        results = run_upstream(args.repo, args.phase)
        print(json.dumps([{"script": r.script, "returncode": r.returncode, "seconds": r.seconds} for r in results], indent=2)); return
    raise RuntimeError(f"unhandled command: {args.command}")


if __name__ == "__main__":
    main()
