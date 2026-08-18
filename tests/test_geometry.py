import sympy as sp

from noble_rei.critical import A
from noble_rei.geometry import are_similar_orbits
from noble_rei.orbits import ONE_D_ORBITS, ZERO_D_ORBITS


def test_rp_phi_is_lower_dimensional_icosahedral_specialization() -> None:
    phi = (1 + sp.sqrt(5)) / 2
    rp = next(spec for spec in ONE_D_ORBITS if spec.name == "rP").vertices()
    specialized = tuple(tuple(sp.simplify(x.subs(A, phi)) for x in point) for point in rp)
    ico = next(spec for spec in ZERO_D_ORBITS if spec.name == "I").vertices()
    assert are_similar_orbits(specialized, ico)
