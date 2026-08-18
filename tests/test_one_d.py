from noble_rei.one_d import enumerate_one_d_orbit
from noble_rei.orbits import ONE_D_ORBITS


def test_rp_specialization_is_not_counted_as_new_polyhedron() -> None:
    spec = next(spec for spec in ONE_D_ORBITS if spec.name == "rP")
    result = enumerate_one_d_orbit(spec)
    assert result.generic_faceting_count == 0
    assert result.count == 0
