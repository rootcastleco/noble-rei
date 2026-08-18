"""Adapter for executing Connor Hill's published GPL enumeration checkout.

This is deliberately a gateway, not a dependency of the independent enumeration core.
It lets a researcher reproduce the canonical pipeline when Wolfram Language is present.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import shutil
import subprocess
import sys
import time


@dataclass(frozen=True, slots=True)
class StepResult:
    script: str
    returncode: int
    seconds: float


PHASES: dict[str, tuple[str, ...]] = {
    "0d": ("enumerate0D.py",),
    "1d": ("enumerate1D-A.py", "enumerate1D-B.wls", "enumerate1D-C.py"),
    "2d": (
        "initialize2D-A.py", "initialize2D-B.wls", "initialize2D-C.py",
        "enumerate2D-A.py", "enumerate2D-B.py", "enumerate2D-C.wls",
        "enumerate2D-D.wls", "enumerate2D-E.py",
    ),
}


def _command_for(script: Path) -> list[str]:
    if script.suffix == ".py":
        return [sys.executable, script.name]
    if script.suffix == ".wls":
        executable = shutil.which("wolframscript")
        if executable is None:
            raise RuntimeError("wolframscript is required for the upstream 1D/2D pipeline")
        return [executable, "-file", script.name]
    raise ValueError(f"unsupported upstream script type: {script}")


def validate_checkout(repo: Path, phase: str) -> tuple[Path, ...]:
    if phase not in PHASES and phase != "all":
        raise ValueError(f"unknown phase: {phase}")
    if not repo.is_dir():
        raise FileNotFoundError(f"upstream repository not found: {repo}")
    phases = ("0d", "1d", "2d") if phase == "all" else (phase,)
    scripts = tuple(repo / name for p in phases for name in PHASES[p])
    missing = [path.name for path in scripts if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"upstream checkout is missing: {', '.join(missing)}")
    return scripts


def run_upstream(repo: Path, phase: str) -> tuple[StepResult, ...]:
    scripts = validate_checkout(repo, phase)
    results: list[StepResult] = []
    for script in scripts:
        command = _command_for(script)
        started = time.perf_counter()
        completed = subprocess.run(command, cwd=repo, check=False)
        elapsed = time.perf_counter() - started
        results.append(StepResult(script.name, completed.returncode, elapsed))
        if completed.returncode != 0:
            raise RuntimeError(f"upstream step failed: {script.name}, exit={completed.returncode}")
    return tuple(results)
