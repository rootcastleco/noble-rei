from noble_rei.reference import PUBLISHED_COUNTS, PUBLISHED_TOTAL, verify_reference_table


def test_published_distribution_sums_to_146() -> None:
    verify_reference_table()
    assert sum(PUBLISHED_COUNTS.values()) == PUBLISHED_TOTAL == 146
