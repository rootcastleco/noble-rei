from noble_rei.enumeration import enumerate_zero_d_orbit
from noble_rei.orbits import ZERO_D_ORBITS


def test_tetrahedral_zero_d_enumeration() -> None:
    spec = next(spec for spec in ZERO_D_ORBITS if spec.name == "T")
    result = enumerate_zero_d_orbit(spec)
    assert result.vertex_count == 4
    assert result.count == 1
