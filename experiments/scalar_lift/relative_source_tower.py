"""Exact sparse source layers in relative cyclotomic extension towers.

All controls have known orders8 or40. Splitting a source preserves its
order and yields coupled sparse layers, not independent smaller Butsons.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from experiments.scalar_lift.fifth_source_assembly import (
    complete_tensor_source,
    exact_column_gram,
    extended_eight_seed,
    root_vectors,
)

ZERO = -1


def extension_data(lower_phase, upper_phase):
    if upper_phase % lower_phase:
        raise ValueError("the lower conductor must divide the upper conductor")
    ratio = upper_phase // lower_phase
    lower = root_vectors(lower_phase)
    upper = root_vectors(upper_phase)
    if ratio < 2 or upper.shape[1] != ratio * lower.shape[1]:
        raise ValueError("the phase ratio must equal the relative field degree")
    embedded = np.zeros((lower_phase, upper.shape[1]), dtype=np.int64)
    embedded[:, ::ratio] = lower
    if not np.array_equal(embedded, upper[::ratio]):
        raise ValueError("the power-basis scalar inclusion failed")
    return ratio, lower, upper


def split_layers(exponents, lower_phase, upper_phase):
    """Use -1 for zero entries; all other entries are phase exponents."""
    ratio, _, _ = extension_data(lower_phase, upper_phase)
    source = np.asarray(exponents, dtype=np.int64)
    if source.ndim != 2 or np.any(source < ZERO) or np.any(source >= upper_phase):
        raise ValueError("canonical exponents or the zero sentinel are required")
    layers = np.full((ratio, *source.shape), ZERO, dtype=np.int64)
    present = source != ZERO
    for index in range(ratio):
        mask = present & (source % ratio == index)
        layers[index][mask] = source[mask] // ratio
    if not np.array_equal((layers != ZERO).sum(axis=0), present.astype(np.int64)):
        raise ValueError("the relative layers do not have disjoint exhaustive supports")
    return layers


def reconstruct_source(layers, lower_phase, upper_phase):
    ratio, _, _ = extension_data(lower_phase, upper_phase)
    layers = np.asarray(layers, dtype=np.int64)
    if layers.ndim != 3 or len(layers) != ratio or np.any(layers < ZERO) or np.any(layers >= lower_phase):
        raise ValueError("invalid sparse relative layers")
    multiplicity = (layers != ZERO).sum(axis=0)
    if np.any(multiplicity > 1):
        raise ValueError("the source layers overlap")
    result = np.full(layers.shape[1:], ZERO, dtype=np.int64)
    for index, layer in enumerate(layers):
        mask = layer != ZERO
        result[mask] = ratio * layer[mask] + index
    return result


def sparse_row_product(first, second, phase):
    """Compute first*second^dagger in the exact lower power basis."""
    first, second = np.asarray(first, dtype=np.int64), np.asarray(second, dtype=np.int64)
    if first.ndim != 2 or second.ndim != 2 or first.shape[1] != second.shape[1]:
        raise ValueError("compatible sparse row matrices are required")
    roots = root_vectors(phase)
    result = np.zeros((len(first), len(second), roots.shape[1]), dtype=np.int64)
    for i, row in enumerate(first):
        for j, other in enumerate(second):
            mask = (row != ZERO) & (other != ZERO)
            result[i, j] = roots[(row[mask] - other[mask]) % phase].sum(axis=0)
    return result


def relative_gram_coefficients(layers, lower_phase, upper_phase):
    ratio, roots, _ = extension_data(lower_phase, upper_phase)
    layers = np.asarray(layers, dtype=np.int64)
    if len(layers) != ratio:
        raise ValueError("the number of layers differs from the relative degree")
    products = [[sparse_row_product(first, second, lower_phase) for second in layers] for first in layers]
    result = np.zeros((ratio, *products[0][0].shape), dtype=np.int64)
    result[0] = sum(products[index][index] for index in range(ratio))
    # Rows of this matrix are the images of the lower power basis under
    # multiplication by zeta_lower^-1; this is the negative-exponent carry.
    negative_carry = roots[(np.arange(roots.shape[1]) - 1) % lower_phase]
    for k in range(1, ratio):
        direct = sum(products[b + k][b] for b in range(ratio - k))
        wrapped = sum(products[b + k - ratio][b] for b in range(ratio - k, ratio))
        result[k] = direct + np.matmul(wrapped, negative_carry)
    return result


def reconstruct_gram(coefficients, lower_phase, upper_phase):
    ratio, lower, upper = extension_data(lower_phase, upper_phase)
    coefficients = np.asarray(coefficients, dtype=np.int64)
    if len(coefficients) != ratio or coefficients.shape[-1] != lower.shape[1]:
        raise ValueError("invalid relative coefficient array")
    result = np.zeros((*coefficients.shape[1:-1], upper.shape[1]), dtype=np.int64)
    for k in range(ratio):
        result[..., k::ratio] = coefficients[k]
    return result


def tower_compatibility(exponents):
    """Compare20->40->80 against the direct relative degree-four split."""
    source = np.asarray(exponents, dtype=np.int64)
    direct = split_layers(source, 20, 80)
    upper_binary = split_layers(source, 40, 80)
    lower_binary = [split_layers(layer, 20, 40) for layer in upper_binary]
    assembled = np.array([lower_binary[j % 2][j // 2] for j in range(4)])
    if not np.array_equal(assembled, direct):
        raise ValueError("the binary tower did not produce the same four sparse layers")
    direct_gram = relative_gram_coefficients(direct, 20, 80)
    upper_gram = relative_gram_coefficients(upper_binary, 40, 80)
    for index in range(4):
        if not np.array_equal(direct_gram[index], upper_gram[index % 2][..., index // 2::2]):
            raise ValueError("the binary tower and direct four-layer Gram coefficients disagree")
    if not np.array_equal(reconstruct_source(direct, 20, 80), source):
        raise ValueError("the four-layer tower did not reconstruct the same source")
    return {"lower_phase": 20, "middle_phase": 40, "upper_phase": 80,
            "direct_layer_index": "j=j80+2*j40",
            "all_four_sparse_layers_identical": True,
            "all_four_Gram_coefficients_identical": True,
            "same_source_reconstructed": True}


def cross_ratio(exponents, phase):
    source = np.asarray(exponents, dtype=np.int64)
    value = int((source[0, 0] + source[1, 4] - source[0, 4] - source[1, 0]) % phase)
    if value != 5 or value % 2 != 1:
        raise ValueError("the control lost its primitive odd cross-ratio exponent")
    return {"rows": [0, 1], "columns": [0, 4], "exponent": value,
            "root_order": phase // 5, "odd_exponent": True,
            "excludes_row_column_phase_gauge_descent_to_half_conductor": True}


def relative_control(source, lower_phase, upper_phase):
    source = np.asarray(source, dtype=np.int64)
    if source.shape[0] != source.shape[1] or np.any(source % 5 == 0):
        raise ValueError("a square source in the admitted fifth alphabet is required")
    layers = split_layers(source, lower_phase, upper_phase)
    ratio = upper_phase // lower_phase
    channel_rules = []
    for index, layer in enumerate(layers):
        present = layer != ZERO
        if np.any((ratio * layer[present] + index) % 5 == 0):
            raise ValueError("a sparse layer violates its shifted fifth-color channel")
        channel_rules.append({"layer_index": index,
                              "allowed_lower_exponent_residues_mod5": [a for a in range(5)
                                                                       if (ratio * a + index) % 5 != 0]})
    if not np.array_equal(reconstruct_source(layers, lower_phase, upper_phase), source):
        raise ValueError("the sparse layers changed the source")
    coefficients = relative_gram_coefficients(layers, lower_phase, upper_phase)
    expected = np.zeros_like(coefficients)
    expected[0, np.arange(len(source)), np.arange(len(source)), 0] = len(source)
    if not np.array_equal(coefficients, expected):
        raise ValueError("the complete relative Gram equations failed")
    upper_row_gram = sparse_row_product(source, source, upper_phase)
    upper_column_gram = exact_column_gram(source, upper_phase)
    expected_upper = np.zeros_like(upper_row_gram)
    expected_upper[np.arange(len(source)), np.arange(len(source)), 0] = len(source)
    if (not np.array_equal(reconstruct_gram(coefficients, lower_phase, upper_phase), upper_row_gram)
            or not np.array_equal(upper_row_gram, expected_upper)
            or not np.array_equal(upper_column_gram, expected_upper)):
        raise ValueError("the independent full upper-field Gram check failed")
    digest = hashlib.sha256(json.dumps(source.tolist(), separators=(",", ":")).encode()).hexdigest()
    report = {"source_order": len(source), "source_exponent_sha256": digest,
            "lower_phase": lower_phase, "upper_phase": upper_phase,
            "relative_degree": upper_phase // lower_phase,
            "zero_sentinel": ZERO, "sparse_layer_exponents": layers.tolist(),
            "layer_support_sizes": [int(np.count_nonzero(layer != ZERO)) for layer in layers],
            "supports_disjoint_and_exhaustive": True,
            "source_reconstruction_verified": True,
            "C0_equals_order_times_identity": True,
            "all_nonconstant_Ck_zero": True,
            "independent_upper_row_and_column_Grams_verified": True,
            "all_source_exponents_nonzero_mod5": True,
            "shifted_fifth_color_rule": "(r*a+j)%5 != 0 for every occupied A_j exponent a",
            "fifth_color_channels": channel_rules,
            "generator_is_primitive_full_conductor_root": True}
    if ratio == 2:
        cross_nonzero = bool(np.any(sparse_row_product(layers[0], layers[1], lower_phase)))
        collapsed = layers.max(axis=0)
        collapsed_gram = sparse_row_product(collapsed, collapsed, lower_phase)
        expected_lower = np.zeros_like(collapsed_gram)
        expected_lower[np.arange(len(source)), np.arange(len(source)), 0] = len(source)
        report["quadratic_cross_coupling"] = {
            "A0_A1_dagger_nonzero": cross_nonzero,
            "collapsed_A0_plus_A1_is_lower_phase_Butson": bool(np.array_equal(collapsed_gram, expected_lower)),
        }
    return report


def target_metadata():
    reports = []
    for order, upper, lower in ((1499, 40, 20), (2999, 80, 40)):
        if (any(order % p == 0 for p in range(2, int(order**0.5) + 1))
                or order % lower != lower - 1):
            raise ValueError("the prime/self-conjugate lower-field obstruction failed")
        reports.append({"source_order": order, "upper_phase": upper, "excluded_lower_phase": lower,
                        "prime_is_minus_one_mod_lower_phase": True,
                        "lower_Butson_source_excluded_by_odd_determinant_norm_valuation": True,
                        "any_upper_source_must_be_nongauge_descending": True,
                        "quadratic_A_B_dagger_must_be_nonzero": True,
                        "reason": "If AB^dagger=0, disjoint full support makes A+B a lower-phase Butson matrix of the same forbidden order.",
                        "upper_source_constructed": False})
    if 1499 % 40 == 39:
        raise ValueError("the order1499 phase40 case was incorrectly excluded")
    reports.append({"source_order": 1499, "upper_phase": 80,
                    "descent_to_phase40_not_excluded_by_this_test": True,
                    "existence_at_phase40_or80_asserted": False})
    return reports


def certificate():
    controls = []
    for phase in (40, 80):
        for label, source in (("known_order40_tensor_source", complete_tensor_source(phase)),
                              ("known_order8_seed", extended_eight_seed(phase))):
            lower_phases = (20,) if phase == 40 else (40, 20)
            reports = [relative_control(source, lower, phase) for lower in lower_phases]
            controls.append({"name": label, "source_order": len(source), "phase": phase,
                             "source_exponents": source.tolist(), "cross_ratio": cross_ratio(source, phase),
                             "relative_checks": reports,
                             "binary_tower_check": None if phase == 40 else tower_compatibility(source)})
    return {"schema": "relative-source-tower-v1",
            "coefficient_formula": {
                "source": "H=sum_j eta^j A_j, eta=zeta_Q, eta^r=zeta_q, r=Q/q=[K_Q:K_q]",
                "C0": "sum_j A_j A_j^dagger",
                "Ck": "sum_(b=0)^(r-1-k) A_(b+k) A_b^dagger + zeta_q^-1 sum_(b=r-k)^(r-1) A_(b+k-r) A_b^dagger",
                "Hadamard_condition": "C0=nI and Ck=0 for1<=k<r",
                "individual_sparse_layers_are_Not_claimed_Butson": True,
                "matrix_order_preserved": True},
            "controls": controls, "prime_target_field_requirements": target_metadata(),
            "new_unknown_source_constructed": False,
            "new_real_Hadamard_order_constructed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--output", type=Path)
    group.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = certificate()
    if args.check:
        if json.loads(args.check.read_text()) != report:
            raise ValueError("the relative-source tower certificate differs from replay")
    else:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"known_control_orders": [8, 40], "relative_tower": [20, 40, 80],
                      "all_exact_relative_Gram_checks_passed": True,
                      "new_unknown_source_constructed": False}))


if __name__ == "__main__":
    main()
