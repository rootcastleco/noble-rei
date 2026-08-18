import sympy as sp

from noble_rei.critical import A, B, configuration_entry, merge_planes, resultant_candidates_2d


def test_configuration_entry_detects_coplanarity() -> None:
    p = (sp.Integer(0), sp.Integer(0), sp.Integer(0)); q = (sp.Integer(1), sp.Integer(0), sp.Integer(0)); r = (sp.Integer(0), sp.Integer(1), sp.Integer(0)); s = (sp.Integer(1), sp.Integer(1), sp.Integer(0))
    assert configuration_entry(p, q, r, s) == 0


def test_plane_merge_uses_shared_vertex_pairs() -> None:
    merged = merge_planes([(1, 2, 3), (1, 2, 4), (7, 8, 9)])
    assert frozenset({0, 1, 2, 3, 4}) in merged
    assert frozenset({0, 7, 8, 9}) in merged


def test_2d_resultant_adapter_eliminates_requested_variable() -> None:
    f = A + B - 1; g = A - B
    eliminate_b, eliminate_a = resultant_candidates_2d(f, g)
    assert sp.simplify(eliminate_b - (2 * A - 1)) == 0
    assert B not in eliminate_b.free_symbols
    assert A not in eliminate_a.free_symbols
