from __future__ import annotations

from experiments.scalar_lift.paley_coefficient_forms import normalized_coefficient_forms


def test_complete_coefficient_enumeration_recovers_the_four_norm_triples():
    report = normalized_coefficient_forms()
    assert len(report["normalized_first_rows"]) == 18
    expected = []
    for a, b, c in ((1, 4, 22), (1, 10, 20), (7, 14, 16), (17, 4, 14)):
        for sign in (1, -1):
            expected.append(
                [
                    [a, 0, b, c],
                    [0, a, sign * c, -sign * b],
                    [b, sign * c, -a, 0],
                    [c, -sign * b, 0, -a],
                ]
            )
    assert report["matrices"] == expected
    # All entries obey the bounds used in the exhaustive argument; no assumed
    # zero pattern was used to obtain the two zero off-diagonal pairs.
    assert all(
        sum(matrix[i][j] == 0 for i in range(4) for j in range(i + 1, 4)) == 2
        for matrix in report["matrices"]
    )
