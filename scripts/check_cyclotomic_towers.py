#!/usr/bin/env python3
"""Replay exact two-stage lifts and audit target factorizations and coefficient fields."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sympy import I, Matrix, eye, simplify, totient

from experiments.scalar_lift.cyclotomic_towers import (
    auxiliary_degree_167_subfield,
    cyclotomic_orders_with_degree_dividing,
    target_factor_audit,
)
from experiments.scalar_lift.inverse_blocks import (
    known_order_thirteen_source,
    triple_unreal_source,
)
from hadamard.cgw import is_cgw_exponent_matrix
from hadamard.dense_lift import (
    TWELFTH_POWER_QUATERNARY_BLOCKS,
    lift_quaternary_hadamard,
    unreal_tenth_to_quaternary,
    unreal_twelfth_to_quaternary,
)
from hadamard.inverse_lift import is_sign_hadamard
from scripts.check_inverse_dense_lifts import SNAPSHOT, _snapshot_candidates

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = ROOT / "certificates/scalar-lift/cyclotomic-towers.json"


def _control(name, source, phase, transfer):
    n = len(source)
    if not is_cgw_exponent_matrix(source, phase, n):
        raise ValueError("a tower source failed its exact original-field Gram identity")
    middle = transfer(source, phase)
    if not is_cgw_exponent_matrix(middle, 4, len(middle)):
        raise ValueError("the intermediate quaternary matrix is not Hadamard")
    output = lift_quaternary_hadamard(middle)
    h = Matrix(output)
    if not is_sign_hadamard(output) or h * h.T != len(output) * eye(len(output)):
        raise ValueError("the iterated transfer failed its exact sign or integer Gram check")
    return {
        "name": name,
        "source_phase": phase,
        "source_order": n,
        "stage_degrees": [len(middle) // n, 2],
        "intermediate_phase": 4,
        "intermediate_order": len(middle),
        "intermediate_exponents": middle,
        "real_output_order": len(output),
        "source_and_intermediate_cyclotomic_grams_exact": True,
        "real_signs_and_integer_gram_exact": True,
        "known_order_control": True,
    }


def _summary(audit):
    return {
        "target_order": audit["target_order"],
        "prime_factorization": audit["target_prime_factorization"],
        "arithmetic_survivors": [
            {
                "source_order": branch["source_order"],
                "block_order": branch["block_order"],
                "phases": branch["surviving_phases"],
                "field_degrees": [int(totient(k)) for k in branch["surviving_phases"]],
                "gaussian_terminal_degree_compatible_phases": [
                    field["phase"]
                    for field in branch["coefficient_fields"]
                    if not field["obstructions"]
                    and field["gaussian_terminal_degree_test"]["passes_necessary_degree_test"]
                ],
            }
            for branch in audit["branches"]
            if branch["surviving_phases"]
        ],
        "necessary_tests_only": True,
    }


def build_report():
    s = Matrix([[I, 1 + I], [-1 + I, -I]])
    t = I * (eye(2) - s) / 2

    def clean(x):
        return x.applyfunc(lambda z: simplify(z.expand()))

    if not (
        clean(t * t.conjugate().T) == eye(2)
        and clean(t**4 - t**2 + eye(2)) == Matrix.zeros(2)
        and clean(t**3) == I * eye(2)
    ):
        raise ValueError("the relative cyclotomic representation is invalid")
    for e, exponents in TWELFTH_POWER_QUATERNARY_BLOCKS.items():
        block = Matrix([[I**a for a in row] for row in exponents])
        if clean((1 + I) * t**e) != block:
            raise ValueError("a stored twelfth-root scalar block differs from the field image")

    inverse_path = ROOT / "certificates/scalar-lift/inverse-dense-lifts.json"
    inverse = json.loads(inverse_path.read_text())
    controls = [
        _control(p["name"], p["exponents"], 12, unreal_twelfth_to_quaternary)
        for p in inverse["compatible_torsion_points"]
        if p["phase"] == 12
    ]
    source_three = tuple(tuple(4 if i == j else 2 for j in range(3)) for i in range(3))
    controls.extend(
        _control(name, source, 6, unreal_twelfth_to_quaternary)
        for name, source in (
            ("unreal source order nine", triple_unreal_source(source_three)),
            ("catalogue source M13A", known_order_thirteen_source()),
        )
    )
    source_five = tuple(tuple((3 * (i + j) ** 2 + 1) % 5 for j in range(5)) for i in range(5))
    controls.append(_control("chirped Fourier source", source_five, 5, unreal_tenth_to_quaternary))

    _, candidates = _snapshot_candidates()
    audits = [target_factor_audit(n) for n in candidates]
    spotlight = {2396, 4676, 8684, 9308, 6876, 4792, 5752}
    audit6876 = next(a for a in audits if a["target_order"] == 6876)
    branch36 = next(b for b in audit6876["branches"] if b["block_order"] == 36)
    higher_fields = [f for f in branch36["coefficient_fields"] if f["phase"] in (38, 76, 114)]
    if len(higher_fields) != 3 or any(f["survives_necessary_tests"] for f in higher_fields):
        raise ValueError("the higher-degree branch escaped the proved flat-alphabet obstruction")
    return {
        "scope": "Prime-factor and cyclotomic-field filters for fixed-frame scalar towers",
        "snapshot_sha256": hashlib.sha256(SNAPSHOT.read_bytes()).hexdigest(),
        "inverse_source_report_sha256": hashlib.sha256(inverse_path.read_bytes()).hexdigest(),
        "candidate_status": "dated corrected snapshot; no current priority claim",
        "candidate_orders_audited": len(audits),
        "composition_degree": "product of the stage degrees",
        "composition_frame": "rho_low(Z_high) (I tensor Z_low)",
        "composition_preserves_gram_when_each_frame_gram_is_scalar": True,
        "each_stage_entry_domain_must_be_checked": True,
        "known_two_stage_controls": controls,
        "block_668_all_root_orders_passing_degree_test": cyclotomic_orders_with_degree_dividing(
            668
        ),
        "auxiliary_cyclotomic_subfield_degree_witness": auxiliary_degree_167_subfield(),
        "all_candidate_factor_summaries": [_summary(audit) for audit in audits],
        "spotlight_factor_audits": [a for a in audits if a["target_order"] in spotlight],
        "gaussian_final_stage_for_block_668": {
            "stage_degrees": [334, 2],
            "first_stage_coefficient_field": "Q(i)",
            "required_quaternary_pair_order": 334,
            "required_identity": "UU*=VV*=334I and UV*+VU*=334I",
            "sufficient_for_real_668_pair": True,
            "pair_constructed": False,
            "target_orders_if_constructed": [4676, 8684],
        },
        "dyadic_flat_alphabet_theorem": {
            "block_order_congruence": "d=4 mod 8",
            "source_order": "odd n>1",
            "used_source_alphabets_after_global_phase": [6, 12],
            "scope": "uniform faithful field scalar lifts with one common frame",
            "candidate_targets_to_which_this_applies": sum(n % 8 == 4 for n in candidates),
            "unused_auxiliary_extensions_remain_allowed": True,
        },
        "higher_degree_branch_audit": {
            "target_order": 6876,
            "source_order": 191,
            "block_order": 36,
            "degree_norm_candidate_phases": [38, 76, 114],
            "field_degrees": [18, 36, 36],
            "remaining_phases_after_flat_source_reduction": [],
            "excluded_for_uniform_faithful_field_scalar_lifts": True,
            "field_tests": higher_fields,
            "source_existence_decided": False,
        },
        "new_real_hadamard_order_claimed": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    if args.output is not None and args.check is not None:
        parser.error("use either --output or --check")
    report = build_report()
    encoded = json.dumps(report, indent=2) + "\n"
    if args.check is not None:
        if json.loads(encoded) != json.loads(args.check.read_text()):
            raise SystemExit("cyclotomic tower report replay differs")
    else:
        path = args.output or DEFAULT_REPORT
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(encoded)
    print(
        json.dumps(
            {
                "candidate_orders_audited": report["candidate_orders_audited"],
                "known_control_real_orders": [
                    c["real_output_order"] for c in report["known_two_stage_controls"]
                ],
                "higher_degree_branch_excluded": report["higher_degree_branch_audit"][
                    "excluded_for_uniform_faithful_field_scalar_lifts"
                ],
                "new_real_order_claimed": False,
            }
        )
    )


if __name__ == "__main__":
    main()
