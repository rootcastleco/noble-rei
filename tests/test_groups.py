from noble_rei.point_groups import matrix_group


def test_expected_point_group_orders() -> None:
    assert len(matrix_group("332")) == 12
    assert len(matrix_group("*332")) == 24
    assert len(matrix_group("3*2")) == 24
    assert len(matrix_group("432")) == 24
    assert len(matrix_group("*432")) == 48
    assert len(matrix_group("532")) == 60
    assert len(matrix_group("*532")) == 120
