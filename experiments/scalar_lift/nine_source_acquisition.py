"""Acquire a previously saved literal BH(9,6), without constructing a source.

The source is copied from an existing exact certificate and bound to its
literal record hash. The degree-191 arithmetic is a field preflight, not
a flatness assertion or a real Hadamard construction.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from functools import cache
from pathlib import Path

import numpy as np
from sympy import isprime, n_order, totient

from experiments.scalar_lift.forward_butson_descent import verify_source
from experiments.scalar_lift.forward_source_registry import VerifiedSource
from experiments.scalar_lift.scalar_restriction_obstructions import alphabet_obstruction

ROOT = Path(__file__).resolve().parents[2]
PRIOR_PATH = ROOT / "certificates/scalar-lift/inverse-dense-lifts.json"
PRIOR_RECORD_HASH = "b78a2e434048f4ce2e197dff4389b55ded30605ffb1d56bc4e67cb0345198e4f"
PRIOR_FILE_HASH_AT_ACQUISITION = "0bc5859820c53953f3a9de19cd8e70d425a782e11da16b2c3632f30b8f2b05f4"
SOURCE_HASH = "ef412f4ec5d8416f900c1e3d6059cb3d9771e1ae6bcd9f9fa378aed78abda127"
IDENTIFIER = "saved-unreal9-q6"
PERIOD_CANDIDATES = (
    (1431737, 27861, 788, 25),
    (1519597, 2669, 402, 5),
    (1819849, 2671, 220, 3),
)


def sources():
    return [{"id": IDENTIFIER, "order": 9, "phase": 6,
             "family": "previously-saved-literal-input"}]


@cache
def load_source(identifier=IDENTIFIER):
    if identifier != IDENTIFIER:
        raise ValueError("unknown saved order-nine input")
    encoded = PRIOR_PATH.read_bytes()
    records = json.loads(encoded)["additional_verified_unreal_sixth_root_sources"]
    record = next(row for row in records if row["name"] == "unreal-order-nine-by-tripling")
    digest = hashlib.sha256(json.dumps(record, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if digest != PRIOR_RECORD_HASH:
        raise ValueError("the preexisting literal source record changed")
    array = np.asarray(record["exponents"], dtype=np.int64)
    if array.shape != (9, 9) or record["phase"] != 6:
        raise ValueError("the saved source has unexpected dimensions or phase")
    proof = verify_source(array, 6)
    if proof["exponent_sha256"] != SOURCE_HASH:
        raise ValueError("the literal source hash changed")
    used = tuple(map(int, np.unique(array)))
    if used != (1, 2, 4, 5):
        raise ValueError("the saved source is not in the recorded unreal alphabet")
    gate = alphabet_obstruction(6, used)
    if gate["excluded"] or not gate["complete_for_odd_integer_relations"]:
        raise ArithmeticError("the complete odd-relation control changed")
    proof = {**proof, "acquisition": "read literal entries from an existing certificate",
             "prior_certificate": str(PRIOR_PATH.relative_to(ROOT)),
             "prior_record_sha256": digest,
             "prior_file_sha256_at_acquisition": PRIOR_FILE_HASH_AT_ACQUISITION,
             "prior_file_sha256_now": hashlib.sha256(encoded).hexdigest(),
             "prior_record_name": record["name"],
             "family_reference_recorded_by_prior_artifact": record["known_family_reference"],
             "source_construction_performed": False,
             "odd_integer_relation_gate": gate}
    array = array.astype(np.uint16)
    array.setflags(write=False)
    return VerifiedSource(identifier, 9, 6, SOURCE_HASH, used, proof, array,
                          lambda i, j: int(array[i, j]))


def period_preflight():
    r = 191
    result = []
    for ell, a, b, t in PERIOD_CANDIDATES:
        h, remainder = divmod(ell - 1, r)
        if (remainder or h % 2 or not isprime(ell) or pow(2, h, ell) != 1
                or a*a + r*b*b != ell*t*t or a % 2 != 1 or b % 2 or t % 2 != 1):
            raise ArithmeticError("the literal period-field arithmetic failed")
        result.append({"real_period_degree": r, "prime_conductor": ell,
                       "subgroup_order": h, "order_of_two": int(n_order(2, ell)),
                       "two_splits_completely": True,
                       "norm_witness_a_b_t": [a, b, t],
                       "norm_equality_value": a*a + r*b*b,
                       "trace_basis_parity": "a,t odd; b even; resulting S is identity modulo2",
                       "compositum": "Q(zeta12) times the real degree191 period field",
                       "compositum_degree": 4*r,
                       "degree_reason": "coprime field degrees4 and191",
                       "rational_source_matrix_order": 9*4*r,
                       "actual_period_multiplication_not_built_by_this_module": True,
                       "sign_flattenability_proved": False})
    return result


def certificate():
    source = load_source()
    array = source.exponents.astype(np.int64)
    cycles = (array + array[0, 0] - array[:, 0, None] - array[None, 0, :]) % 6
    projective_order = 6 // math.gcd(6, *map(int, cycles.flat))
    return {"schema": "saved-order-nine-source-acquisition-v1",
            "source": source.metadata(), "literal_exponents": source.exponents.tolist(),
            "projective_root_order": projective_order,
            "projective_field_degree": int(totient(projective_order)),
            "phase_field_degree_divides764": 764 % int(totient(6)) == 0,
            "target_arithmetic": "6876=9*764=9*4*191",
            "period_field_preflight": period_preflight(),
            "source_construction_tree_edge_added": False,
            "real_Hadamard_of_order6876_constructed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument("--report", type=Path)
    actions.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = certificate()
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    elif json.loads(args.check.read_text()) != report:
        raise ValueError("order-nine acquisition certificate failed replay")
    print(json.dumps({"source_id": IDENTIFIER, "source_Gram_verified": True,
                      "literal_source_sha256": SOURCE_HASH,
                      "degree191_conductors": [r["prime_conductor"] for r in report["period_field_preflight"]],
                      "real_Hadamard_constructed": False}))


if __name__ == "__main__":
    main()
