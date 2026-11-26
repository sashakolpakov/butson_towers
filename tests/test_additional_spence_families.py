import json
import hashlib
import re
from pathlib import Path

import numpy as np
import pytest

from experiments.scalar_lift import additional_spence_families as family
from experiments.scalar_lift import williamson_spence_family as original
from experiments.scalar_lift.williamson_scalar_endpoints import (
    ImplicitScalarEndpoint, verify_orthogonal_design,
)


def test_complete_additional_sources_and_finite_members_replay():
    before = dict(original.PARAMETERS)
    path = Path(__file__).resolve().parents[1] / "certificates/scalar-lift/additional-spence-families.json"
    report = family.certificate()
    assert report == json.loads(path.read_text())
    assert [x["source_order"] for x in report["first_expanded_members"]] == [15, 1431, 4753, 5995, 6441, 9453]
    assert [x["endpoint_order"] for x in report["higher_original_family_members"]] == [2985892, 38816596, 291145924, 6564924]
    assert original.PARAMETERS == before == {13: (13, 1), 17: (17, 1), 61: (61, 1), 101: (101, 1), 121: (11, 2)}
    with pytest.raises(ValueError):
        original.SpenceSource(5, 1)


def test_additional_literal_source_rejects_a_false_group_norm():
    source = family.source(137, 1)
    bad = source.coefficient_rows()
    bad[0, 0] *= -1
    with pytest.raises(ValueError, match="group-ring norm"):
        original.exact_group_norms(bad, source.group_dimensions, 137, 1)


def test_additional_source_domain_and_implicit_entry_access():
    source = family.source(5, 2)
    assert source.order == 75 and source.group_dimensions == (3, 5, 5)
    assert source.entry(2, 74, 12) == int(source.row(2, 74)[12])
    assert family.source(137, 12).entry(1, 0, 100) in (-1, 1)
    for q, r in ((True, 0), (5, False), (5, -1), (73, 0)):
        with pytest.raises(ValueError):
            family.source(q, r)
    for channel, i, j in ((True, 0, 0), (4, 0, 0), (0, -1, 0), (0, 0, source.order)):
        with pytest.raises(ValueError):
            source.entry(channel, i, j)
    with pytest.raises(ValueError, match="NTT replay bound"):
        family.finite_member(137, 12)


def test_additional_source_substitution_has_full_exact_H420_Gram():
    t = np.zeros((4, 7), dtype=np.int64)
    t[0, 0] = 1
    t[0, [1, 2, 4]] = -1
    t[1, 3] = t[2, 6] = t[3, 5] = 1
    x, _ = verify_orthogonal_design(t)
    h = ImplicitScalarEndpoint(x, family.source(5, 1))
    assert h.order == 420
    packed = []
    for i in range(h.order):
        row = h.row(i)
        assert np.all(abs(row) == 1)
        assert all(h.entry(i, j) == int(row[j]) for j in (0, 14, 15, h.order-1))
        packed.append(int.from_bytes(np.packbits(row < 0, bitorder="little").tobytes(), "little"))
    assert all((row ^ other).bit_count() == h.order//2
               for i, row in enumerate(packed) for other in packed[:i])


def test_completed_inventory_counts_only_positive_distinct_orders():
    root = Path(__file__).resolve().parents[1]
    ledger = json.loads((root / "certificates/scalar-lift/completed-order-inventory.json").read_text())
    inputs = {}
    for path, digest in ledger["source_certificate_sha256"].items():
        data = (root / path).read_bytes()
        assert hashlib.sha256(data).hexdigest() == digest
        inputs[path] = json.loads(data)

    def selected(data, selector):
        for key, index in re.findall(r"([^.[\]]+)|\[(\d+)\]", selector):
            data = data[int(index)] if index else data[key]
        return data

    orders = ledger["orders"]
    assert len({record["order"] for record in orders}) == len(orders) == 65
    new = [x["order"] for x in orders if x["status"] == "new_to_best_of_authors_knowledge"]
    known = [x["order"] for x in orders if x["status"] == "known"]
    assert len(new) == 34 and len(known) == 31
    assert new == ledger["completed_new_real_orders"]
    assert known == ledger["completed_known_real_orders"]
    for record in orders:
        assert record["construction_complete"] and record["positive_certificates"]
        for ref in record["positive_certificates"]:
            proof = selected(inputs[ref["path"]], ref["selector"])
            assert proof.get("endpoint_order", proof.get("output_order")) == record["order"]
            assert (any(proof.get(key) is True for key in (
                "all_coefficient_pair_cross_and_endpoint_integer_Grams_verified",
                "complete_algebraic_endpoint_Gram_verified", "endpoint_Gram_proved_exactly",
                "full_endpoint_Gram_verified_algebraically")) or proof.get("full_endpoint_Gram_proof"))
        for ref in record["status_sources"]:
            status = selected(inputs[ref["path"]], ref["selector"])["status"]
            assert status == "known" or status.startswith("new_to_best")
            assert ("known" if status == "known" else "new_to_best_of_authors_knowledge") == record["status"]
