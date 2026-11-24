"""Search proofs, exact terminal recipes, deduplication and crash recovery."""

import json

import numpy as np
import pytest

from experiments.scalar_lift.forward_descent_search import (
    ForwardSearch,
    affine_fit,
    alphabet_after,
    endpoint_entry,
    signed_triangle,
    weighted_metric_record,
)
from experiments.scalar_lift.forward_source_registry import sources


def catalogue(identifier):
    return [row for row in sources() if row["id"] == identifier]


def test_affine_domains_and_exact_uniform_signed_obstruction():
    fit = affine_fit((0, 1), 3, (1, 2))
    assert fit["found"] and fit["transformed_alphabet"] == [1, 2]
    assert not affine_fit(range(3), 3, (1, 2))["found"]
    witness = signed_triangle(range(3), 3)
    assert witness and witness["uniform_linear_raw_alphabet_sign_frame_excluded"]
    assert signed_triangle((0, 1), 3) is None
    assert set(alphabet_after(range(8), 8, 2, "primary")) == set(range(4))


def test_weighted_nonflat_actions_are_exact_and_retained():
    for phase in (3, 5, 24, 120):
        record = weighted_metric_record(phase)
        gram = np.asarray(record["trace_metric"])
        assert np.array_equal(gram, gram.T)
        assert np.linalg.eigvalsh(gram).min() > 0
        assert not record["flat_frame_constructed_by_this_record"]


def test_registered_forward_search_resume_and_actual_sign_recipe(tmp_path):
    path = tmp_path / "forward.sqlite3"
    source_catalogue = catalogue("fifth-5-q20")
    search = ForwardSearch(path, catalogue=source_catalogue)
    try:
        first = search.run(seconds=10, max_steps=1)
        assert first["progress"]["verified_sources"] == 1
        assert first["states"]["pending"] == 1
        with pytest.raises(RuntimeError, match="another forward"):
            ForwardSearch(path, catalogue=source_catalogue)
    finally:
        search.close()
    resumed = ForwardSearch(path, catalogue=source_catalogue)
    try:
        result = resumed.run(seconds=60)
        assert result["progress"]["registered_branches_finished"]
        assert result["verified_endpoint_recipes"] > 0
        identifier, raw = resumed.db.execute("SELECT id,record FROM endpoints ORDER BY rowid LIMIT 1").fetchone()
        record = json.loads(raw)
        assert record["order"] == 40
        assert result["states"]["open"] > 0
        before = result["verified_endpoint_recipes"]
        assert resumed.run(seconds=1)["verified_endpoint_recipes"] == before
    finally:
        resumed.close()
    signs = np.asarray([[endpoint_entry(path, identifier, i, j) for j in range(40)] for i in range(40)])
    assert set(np.unique(signs)) == {-1, 1}
    assert np.array_equal(signs @ signs.T, 40 * np.eye(40, dtype=np.int64))


def test_interruption_rolls_back_source_proof_and_children(tmp_path):
    search = ForwardSearch(tmp_path / "rollback.sqlite3", catalogue=catalogue("fifth-5-q20"))
    evaluate = search.evaluate

    def interrupt(recipe, identifier):
        evaluate(recipe, identifier)
        raise KeyboardInterrupt

    search.evaluate = interrupt
    try:
        with pytest.raises(KeyboardInterrupt):
            search.run(seconds=20)
        assert search._get("progress")["verified_sources"] == 0
        assert search.db.execute("SELECT COUNT(*) FROM states").fetchone()[0] == 1
        assert search.db.execute("SELECT proof FROM sources").fetchone()[0] is None
    finally:
        search.close()
