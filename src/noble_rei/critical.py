"""Symbolic critical-configuration machinery for 1D noble-orbit families."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from fractions import Fraction
from typing import TypeAlias

import sympy as sp

from .algebra import Permutation, Vec3

Triple: TypeAlias = tuple[int, int, int]
CriticalPlane: TypeAlias = frozenset[int]
VolumeConfiguration: TypeAlias = dict[sp.Expr, tuple[Triple, ...]]
A = sp.Symbol("a", positive=True, real=True)
B = sp.Symbol("b", positive=True, real=True)


@dataclass(frozen=True, slots=True)
class CriticalFactor1D:
    polynomial: sp.Expr
    positive_roots: tuple[float, ...]
    planes: tuple[CriticalPlane, ...]


def configuration_entry(p: Vec3, q: Vec3, r: Vec3, s: Vec3) -> sp.Expr:
    p0 = sp.ImmutableMatrix(3, 1, p); q0 = sp.ImmutableMatrix(3, 1, q); r0 = sp.ImmutableMatrix(3, 1, r); s0 = sp.ImmutableMatrix(3, 1, s)
    return sp.expand((q0 - p0).dot((r0 - p0).cross(s0 - p0)))

Q25: TypeAlias = tuple[Fraction, Fraction, Fraction, Fraction]
Q25Poly: TypeAlias = tuple[Q25, ...]
_Q0: Q25 = (Fraction(0), Fraction(0), Fraction(0), Fraction(0))
_SQRT2 = sp.sqrt(2); _SQRT5 = sp.sqrt(5); _SQRT10 = sp.sqrt(10)


def _fraction(value: sp.Expr) -> Fraction:
    value = sp.cancel(value)
    if value.is_Rational is not True:
        raise ValueError(f"coefficient is not rational after Q(sqrt2,sqrt5) expansion: {value}")
    rational = sp.Rational(value)
    return Fraction(int(rational.p), int(rational.q))


def _q25_from_expr(expr: sp.Expr) -> Q25:
    expanded = sp.expand(sp.radsimp(expr)); collected = sp.collect(expanded, [_SQRT2, _SQRT5, _SQRT10], evaluate=False)
    allowed = {sp.Integer(1), _SQRT2, _SQRT5, _SQRT10}
    if any(key not in allowed for key in collected):
        raise ValueError(f"expression escaped Q(sqrt2,sqrt5): {expr} -> {collected}")
    return (_fraction(collected.get(sp.Integer(1), sp.Integer(0))), _fraction(collected.get(_SQRT2, sp.Integer(0))), _fraction(collected.get(_SQRT5, sp.Integer(0))), _fraction(collected.get(_SQRT10, sp.Integer(0))))


def _q25_to_expr(value: Q25) -> sp.Expr:
    def r(x: Fraction) -> sp.Rational: return sp.Rational(x.numerator, x.denominator)
    return r(value[0]) + r(value[1]) * _SQRT2 + r(value[2]) * _SQRT5 + r(value[3]) * _SQRT10

def _q25_add(x: Q25, y: Q25) -> Q25: return (x[0]+y[0], x[1]+y[1], x[2]+y[2], x[3]+y[3])
def _q25_sub(x: Q25, y: Q25) -> Q25: return (x[0]-y[0], x[1]-y[1], x[2]-y[2], x[3]-y[3])
def _q25_mul(x: Q25, y: Q25) -> Q25:
    a,b,c,d=x; e,f,g,h=y
    return (a*e + 2*b*f + 5*c*g + 10*d*h, a*f + b*e + 5*(c*h + d*g), a*g + c*e + 2*(b*h + d*f), a*h + d*e + b*g + c*f)
def _q25_conj2(x: Q25) -> Q25: a,b,c,d=x; return (a,-b,c,-d)
def _q25_conj5(x: Q25) -> Q25: a,b,c,d=x; return (a,b,-c,-d)
def _q25_conj25(x: Q25) -> Q25: a,b,c,d=x; return (a,-b,-c,d)
def _q25_inv(x: Q25) -> Q25:
    if x == _Q0: raise ZeroDivisionError("cannot invert zero in Q(sqrt2,sqrt5)")
    numerator = _q25_mul(_q25_mul(_q25_conj2(x), _q25_conj5(x)), _q25_conj25(x)); norm = _q25_mul(x, numerator)
    if norm[1:] != (Fraction(0), Fraction(0), Fraction(0)) or norm[0] == 0: raise ArithmeticError(f"number-field norm invariant failed: {x} -> {norm}")
    n = norm[0]; return tuple(component / n for component in numerator)  # type: ignore[return-value]

def _q25_linear_coefficients(expr: sp.Expr, variable: sp.Symbol) -> tuple[Q25, Q25]:
    poly = sp.Poly(expr, variable)
    if poly.degree() > 1: raise ValueError(f"expected coordinates linear in {variable}, got {expr}")
    return _q25_from_expr(poly.nth(0)), _q25_from_expr(poly.nth(1))
def _qpoly_add(left: Q25Poly, right: Q25Poly) -> Q25Poly:
    size=max(len(left),len(right)); return tuple(_q25_add(left[i] if i<len(left) else _Q0, right[i] if i<len(right) else _Q0) for i in range(size))
def _qpoly_sub(left: Q25Poly, right: Q25Poly) -> Q25Poly:
    size=max(len(left),len(right)); return tuple(_q25_sub(left[i] if i<len(left) else _Q0, right[i] if i<len(right) else _Q0) for i in range(size))
def _qpoly_mul(left: Q25Poly, right: Q25Poly) -> Q25Poly:
    output=[_Q0 for _ in range(len(left)+len(right)-1)]
    for i,a_coeff in enumerate(left):
        if a_coeff == _Q0: continue
        for j,b_coeff in enumerate(right):
            if b_coeff != _Q0: output[i+j]=_q25_add(output[i+j], _q25_mul(a_coeff,b_coeff))
    return tuple(output)
def _qpoly_monic_key(coefficients: Q25Poly) -> Q25Poly:
    degree=-1
    for index in range(len(coefficients)-1,-1,-1):
        if coefficients[index] != _Q0: degree=index; break
    if degree<0: return ()
    inv=_q25_inv(coefficients[degree]); return tuple(_q25_mul(coefficients[i],inv) for i in range(degree+1))


def volume_configuration(orbit: tuple[Vec3, ...], variable: sp.Symbol = A, symmetry_action: tuple[Permutation, ...] | None = None) -> VolumeConfiguration:
    if len(orbit) < 4: raise ValueError("orbit must contain at least four vertices")
    point_coeffs = tuple(tuple(_q25_linear_coefficients(coord, variable) for coord in point) for point in orbit); anchor = point_coeffs[0]
    differences = [tuple((_q25_sub(point[axis][0], anchor[axis][0]), _q25_sub(point[axis][1], anchor[axis][1])) for axis in range(3)) for point in point_coeffs]
    all_triples = tuple(combinations(range(1, len(orbit)), 3)); classes: list[tuple[Triple, tuple[Triple, ...]]] = []
    if symmetry_action:
        seen: set[Triple] = set()
        for triple in all_triples:
            if triple in seen: continue
            subset = (0, *triple); members: set[Triple] = set()
            for permutation in symmetry_action:
                image = {permutation[index] for index in subset}
                if 0 not in image: continue
                image.remove(0)
                if len(image) != 3: raise RuntimeError("permutation collapsed a 4-point configuration")
                members.add(tuple(sorted(image)))
            if not members or triple not in members: raise RuntimeError("configuration orbit failed to contain its representative")
            seen.update(members); classes.append((triple, tuple(sorted(members))))
        if len(seen) != len(all_triples): raise RuntimeError("symmetry quotient did not cover all anchored configurations")
    else: classes = [(triple, (triple,)) for triple in all_triples]
    pair_cross: dict[tuple[int, int], tuple[tuple, tuple, tuple]] = {}
    for j, k in combinations(range(1, len(orbit)), 2):
        r, t = differences[j], differences[k]
        pair_cross[(j, k)] = (_qpoly_sub(_qpoly_mul(r[1], t[2]), _qpoly_mul(r[2], t[1])), _qpoly_sub(_qpoly_mul(r[2], t[0]), _qpoly_mul(r[0], t[2])), _qpoly_sub(_qpoly_mul(r[0], t[1]), _qpoly_mul(r[1], t[0])))
    raw_groups: dict[tuple, list[Triple]] = defaultdict(list)
    for (i, j, k), members in classes:
        cross0, cross1, cross2 = pair_cross[(j, k)]; q = differences[i]
        coefficients = _qpoly_add(_qpoly_mul(q[0], cross0), _qpoly_mul(q[1], cross1)); coefficients = _qpoly_add(coefficients, _qpoly_mul(q[2], cross2)); raw_groups[coefficients].extend(members)
    grouped_keys: dict[tuple, list[Triple]] = defaultdict(list)
    for raw_key, triples in raw_groups.items(): grouped_keys[_qpoly_monic_key(raw_key)].extend(triples)
    grouped: VolumeConfiguration = {}
    for key, triples in grouped_keys.items():
        expr = sp.Integer(0) if not key else sp.expand(sp.Add(*(_q25_to_expr(coefficient) * variable**degree for degree, coefficient in enumerate(key))))
        grouped[expr] = tuple(triples)
    return grouped


def merge_planes(triples: tuple[Triple, ...] | list[Triple]) -> tuple[CriticalPlane, ...]:
    unique = tuple(sorted(set(triples)))
    if not unique: return ()
    parent = list(range(len(unique))); rank = [0] * len(unique)
    def find(x: int) -> int:
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(x: int, y: int) -> None:
        rx, ry = find(x), find(y)
        if rx == ry: return
        if rank[rx] < rank[ry]: rx, ry = ry, rx
        parent[ry] = rx
        if rank[rx] == rank[ry]: rank[rx] += 1
    pair_owner: dict[tuple[int, int], int] = {}
    for index, triple in enumerate(unique):
        for pair in combinations(triple, 2): union(index, pair_owner.setdefault(pair, index))
    vertices_by_root: dict[int, set[int]] = defaultdict(set)
    for index, triple in enumerate(unique): vertices_by_root[find(index)].update(triple)
    return tuple(sorted((frozenset({0, *vertices}) for vertices in vertices_by_root.values()), key=lambda plane: (len(plane), tuple(sorted(plane)))))


def _number_field() -> sp.polys.domains.AlgebraicField: return sp.QQ.algebraic_field(sp.sqrt(2), sp.sqrt(5))
def _normalize_factor(expr: sp.Expr, variable: sp.Symbol) -> sp.Expr: return sp.expand(sp.Poly(expr, variable, domain=_number_field()).monic().as_expr())
def _positive_numeric_roots(expr: sp.Expr, variable: sp.Symbol) -> tuple[float, ...]:
    poly = sp.Poly(expr, variable, domain=_number_field()); roots = sp.nroots(poly.as_expr(), n=50, maxsteps=200); real_values: list[float] = []
    for root in roots:
        real_part, imag_part = root.as_real_imag(); rv = float(real_part); iv = abs(float(imag_part))
        if iv <= 1e-30 and rv > 0.0 and not any(abs(rv - prior) <= 1e-12 * max(1.0, abs(rv)) for prior in real_values): real_values.append(rv)
    return tuple(sorted(real_values))


def critical_factors_1d(configuration: VolumeConfiguration, variable: sp.Symbol = A) -> tuple[CriticalFactor1D, ...]:
    field = _number_field(); shared: list[Triple] = list(configuration.get(sp.Integer(0), ())); factor_triples: dict[sp.Expr, set[Triple]] = defaultdict(set)
    for expr, triples in configuration.items():
        if expr == 0: continue
        poly = sp.Poly(expr, variable, domain=field)
        if poly.is_zero: shared.extend(triples); continue
        _, factors = sp.factor_list(poly)
        for factor_poly, _multiplicity in factors: factor_triples[_normalize_factor(factor_poly.as_expr(), variable)].update(triples)
    shared_tuple = tuple(shared); results: list[CriticalFactor1D] = []
    for factor, triples in factor_triples.items():
        roots = _positive_numeric_roots(factor, variable)
        if not roots: continue
        all_triples = tuple(sorted(set(triples).union(shared_tuple))); results.append(CriticalFactor1D(factor, roots, merge_planes(all_triples)))
    results.sort(key=lambda item: (sp.degree(item.polynomial, variable), str(item.polynomial)))
    for left, right in combinations(results, 2):
        gcd = sp.gcd(sp.Poly(left.polynomial, variable, domain=field), sp.Poly(right.polynomial, variable, domain=field))
        if gcd.degree() > 0: raise RuntimeError(f"critical factors are not coprime: {left.polynomial}, {right.polynomial}")
    return tuple(results)


def resultant_candidates_2d(f: sp.Expr, g: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    if f == 0 or g == 0: raise ValueError("critical curves must be nonzero")
    return sp.factor(sp.resultant(f, g, B)), sp.factor(sp.resultant(f, g, A))
