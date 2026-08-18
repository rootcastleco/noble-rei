"""Published result invariants used only for verification, never generation."""

from __future__ import annotations

from typing import Final

PUBLISHED_COUNTS: Final[dict[str, int]] = {
    "T": 1, "O": 1, "CO": 0, "C": 1, "I": 4, "ID": 6, "D": 7,
    "tT": 0, "rT": 0, "rP": 0, "tO": 1, "tC": 1, "rC": 1,
    "tI": 17, "tD": 6, "rD": 19,
    "sT": 0, "gT": 0, "gP": 0, "sC": 7, "gC": 3, "sD": 33, "gD": 38,
}

PUBLISHED_TOTAL: Final[int] = 146


def verify_reference_table() -> None:
    total = sum(PUBLISHED_COUNTS.values())
    if total != PUBLISHED_TOTAL:
        raise AssertionError(f"reference table invariant failed: {total} != {PUBLISHED_TOTAL}")
