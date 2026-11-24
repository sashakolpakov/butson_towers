"""Concrete smaller-source rectangular branches toward BH(1499,120).

Known restricted sources of orders567,585 and729 provide actual boundaries.
Root-valued corner and simultaneous cross-row completion are still open.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from functools import cache
from math import gcd
from pathlib import Path

import numpy as np

from experiments.scalar_lift.fifth_source_assembly import exact_column_gram, root_vectors
from hadamard.coefficients import CyclotomicField

ROOT = Path(__file__).resolve().parents[2]
CYCLE = np.array([[0, 1, -1], [-1, 0, 1], [1, -1, 0]], dtype=np.int64)


def branch_parameters(source_order, input_order):
    n, d = int(source_order), int(input_order)
    if n != source_order or d != input_order or (n - d) % 2:
        raise ValueError("the source and input orders must have the same parity")
    s, delta = (n - d) // 2, n - 2 * d
    if not 0 < s <= d or delta <= 0:
        raise ValueError("the rectangular branch requires positive metric and at most d selected rows")
    return {"source_order": n, "input_order": d, "selected_rows": s,
            "omitted_rows": d - s, "selected_metric_norm": delta,
            "complement_metric_norm": n, "phase": 120,
            "known_input_constructor_available": d in (567, 585, 729),
            "known_input_route": ({567: "restricted BH7 followed by four exact triplings",
                                   585: "restricted BH13 followed by two triplings, tensor the admitted fifth-root chirp",
                                   729: "restricted BH3 followed by five exact triplings"}.get(d)),
            "characteristic_zero_source_constructed": False}


def _components(source):
    source = np.asarray(source, dtype=np.int64) % 6
    if np.any(source % 3 == 0):
        raise ValueError("the saved input seed must use only unreal sixth roots")
    a = np.array([0, 1, -1, 0, -1, 1], dtype=np.int64)[source]
    b = np.array([0, 1, 1, 0, -1, -1], dtype=np.int64)[source]
    return a, b


def _triple(source):
    a, b = _components(source)
    aa = np.kron(b, np.ones((3, 3), dtype=np.int64))
    bb = np.kron(a, np.eye(3, dtype=np.int64)) + np.kron(b, CYCLE)
    return np.array([[4, 2], [5, 1]], dtype=np.int64)[(aa + 1) // 2, (bb + 1) // 2]


@cache
def _sixth_source(input_order):
    seed_order, repetitions = {117: (13, 2), 567: (7, 4), 729: (3, 5)}[int(input_order)]
    if seed_order == 13:
        from experiments.scalar_lift.inverse_blocks import known_order_thirteen_source

        source = np.asarray(known_order_thirteen_source(), dtype=np.int64)
    else:
        data = json.loads((ROOT / "certificates/scalar-lift/large-joint-targets.json").read_text())
        source = np.asarray(next(row["exponents"] for row in data["verified_sources"]
                                 if row["order"] == seed_order), dtype=np.int64)
    _components(source)
    for _ in range(repetitions):
        source = _triple(source)
    if source.shape != (input_order, input_order):
        raise ValueError("the smaller-source tripling order is wrong")
    return source


@cache
def known_source(input_order):
    """Actual complete BH(d,120), with raw admitted cubic/fifth colors."""
    if int(input_order) == 585:
        cubic = (20 * _sixth_source(117) + 80) % 120
        fifth = _fifth_chirp()
        source = (cubic[:, None, :, None] + fifth[None, :, None, :]).reshape(585, 585) % 120
    else:
        source = (20 * _sixth_source(int(input_order)) + 104) % 120
    if np.any(source % 3 == 2) or np.any(source % 5 == 0):
        raise ValueError("the constructed smaller source left the permitted alphabet")
    return source


def _fifth_chirp():
    return np.array([[24 * ((3 * (i + j) ** 2 + 1) % 5) for j in range(5)] for i in range(5)], dtype=np.int64)


@cache
def _verify_sixth_gram(input_order):
    source = _sixth_source(int(input_order))
    a, b = _components(source)
    aw = [sum(1 << j for j, value in enumerate(row) if value == 1) for row in a]
    bw = [sum(1 << j for j, value in enumerate(row) if value == 1) for row in b]
    n = len(source)
    for i in range(n):
        for j in range(i):
            aa, bb = n - 2 * (aw[i] ^ aw[j]).bit_count(), n - 2 * (bw[i] ^ bw[j]).bit_count()
            ab, ba = n - 2 * (aw[i] ^ bw[j]).bit_count(), n - 2 * (bw[i] ^ aw[j]).bit_count()
            if aa + 3 * bb or ab != ba:
                raise ValueError("an exact smaller-source Gram entry failed")
    return {"order": n, "all_row_pairs_checked": n * (n - 1) // 2,
            "exact_unreal_component_Gram_verified": True}


@cache
def source_certificate(input_order):
    if int(input_order) == 585:
        factor = _verify_sixth_gram(117)
        fifth = _fifth_chirp()
        gram = exact_column_gram(fifth, 120)
        expected = np.zeros_like(gram)
        expected[np.arange(5), np.arange(5), 0] = 5
        if not np.array_equal(gram, expected):
            raise ValueError("the fifth-root tensor factor has an incorrect Gram")
        upper = known_source(input_order)
        return {"order": 585, "phase": 120,
                "construction": "BH117 tensor chirped BH5, with independent cubic and fifth colors",
                "factor117_exact_Gram": factor, "factor5_complete_exact_Gram_verified": True,
                "complete_585_Gram_verified_by_tensor_identity": True,
                "all_raw_cubic_fifth_alphabet_entries_verified": True,
                "cubic_colors": sorted(set((upper % 3).ravel().tolist())),
                "fifth_colors": sorted(set((4 * upper % 5).ravel().tolist())),
                "source_exponent_sha256": hashlib.sha256(upper.astype("<i8").tobytes()).hexdigest()}
    checked = _verify_sixth_gram(input_order)
    n = checked["order"]
    upper = known_source(input_order)
    return {"order": n, "phase": 120, "all_row_pairs_checked": n * (n - 1) // 2,
            "exact_unreal_component_Gram_verified": True,
            "all_raw_cubic_fifth_alphabet_entries_verified": True,
            "cubic_colors": sorted(set((upper % 3).ravel().tolist())),
            "fifth_colors": sorted(set((4 * upper % 5).ravel().tolist())),
            "source_exponent_sha256": hashlib.sha256(upper.astype("<i8").tobytes()).hexdigest()}


def row_selection(input_order, descriptor=None, source_order=1499):
    data = branch_parameters(source_order, input_order)
    descriptor = {"kind": "stride", "step": 1, "shift": 0} if descriptor is None else descriptor
    if descriptor.get("kind", "stride") not in ("stride", "prefix"):
        raise ValueError("this finite branch uses exact stride row subsets")
    step, shift = int(descriptor.get("step", 1)), int(descriptor.get("shift", 0))
    if gcd(step, input_order) != 1:
        raise ValueError("the stride must permute the input rows")
    selected = [(step * i + shift) % input_order for i in range(data["selected_rows"])]
    if len(set(selected)) != len(selected):
        raise ValueError("a boundary row was repeated")
    return selected


def boundary_problem(input_order, descriptor=None, source_order=1499):
    """Independent frames: T_L from G, T_R=conjugate(T_L).

    Thus T_L and T_R* both have the admitted alphabet. Their Gram
    projectors need not agree; the general corner is not assumed normal.
    """
    parameters = branch_parameters(source_order, input_order)
    source = known_source(input_order)
    rows = row_selection(input_order, descriptor, source_order)
    selected = set(rows)
    omitted = [i for i in range(input_order) if i not in selected]
    return {**parameters, "selected_row_indices": rows, "omitted_row_indices": omitted,
            "T_L": source[rows], "T_R": -source[rows] % 120,
            "V_L": source[omitted], "V_R": -source[omitted] % 120,
            "right_frame_is_conjugate_of_left": True}


def metric_certificate(input_order):
    parameters = branch_parameters(1499, input_order)
    field = CyclotomicField(120)
    eta, imag = field.zeta() ** 15, field.zeta() ** 30
    alpha = 21 + 23 * (eta + eta ** 3)
    if input_order == 585:
        theta = field.zeta() ** 24
        sqrt_five = 1 + 2 * theta + 2 * theta ** 4
        selected = 18 + imag * sqrt_five
        if sqrt_five ** 2 != field.element([5]):
            raise ValueError("the exact fifth-root expression for sqrt5 failed")
    else:
        selected = (5 + 4 * imag) if input_order == 729 else (19 + 2 * imag)
    if alpha * alpha.conjugate() != field.element([1499]):
        raise ValueError("the large source norm witness failed")
    if selected * selected.conjugate() != field.element([parameters["selected_metric_norm"]]):
        raise ValueError("the selected-space norm witness failed")
    if input_order == 585:
        small = CyclotomicField(20)
        zeta = small.zeta()
        witness = 18 + zeta ** 5 * (1 + 2 * zeta ** 4 + 2 * zeta ** 16)
        if witness * witness.conjugate() != small.element([329]):
            raise ValueError("the norm329 subfield witness failed")
        fixed_base_image = sum((int(c) * zeta ** (73 * j) for j, c in enumerate(witness.coefficients)), small.zero())
        if fixed_base_image != witness.conjugate() or fixed_base_image == witness or (73 - 1) % 24:
            raise ValueError("the metric witness did not have its required essential fifth dependence")
        return {"selected_norm": 329,
                "selected_scalar_formula": "18+i*sqrt5",
                "selected_scalar_subfield_conductor": 20,
                "selected_scalar_subfield_coefficients": [int(c) for c in witness.coefficients],
                "fifth_root_formula_for_sqrt5": "1+2*zeta5+2*zeta5^(-1)",
                "source_norm": 1499, "source_scalar_formula": "21+23*(zeta8+zeta8^3)",
                "both_scalar_norms_verified": True,
                "automorphism_fixing_phase24_exponent": 73,
                "that_automorphism_changes_witness_to_distinct_conjugate": True,
                "scalar_metric_is_not_a_root_completion": True,
                "phase24_metric_obstruction": {
                    "metric_order": parameters["selected_rows"], "prime": 47,
                    "prime_residue_mod24": 47 % 24,
                    "determinant_norm_prime_valuation": parameters["selected_rows"],
                    "self_conjugate_prime_has_even_norm_valuation": True,
                    "all_integral_metric_matrices_over_Q_zeta24_excluded": True,
                    "unit_root_row_and_column_gauges_do_not_remove_this_requirement": True}}
    return {"selected_norm": parameters["selected_metric_norm"],
            "selected_scalar_gaussian_coefficients": [5, 4] if input_order == 729 else [19, 2],
            "source_norm": 1499, "source_scalar_formula": "21+23*(zeta8+zeta8^3)",
            "both_scalar_norms_verified": True,
            "scalar_metric_is_not_a_root_completion": True}


@cache
def small_metric_inputs():
    """Verify the actual complete sources used by the finite mixing branches."""
    from experiments.scalar_lift.phase120_common_base import source_thirty_five

    data = json.loads((ROOT / "certificates/scalar-lift/large-joint-targets.json").read_text())
    seven = np.asarray(next(row["exponents"] for row in data["verified_sources"] if row["order"] == 7))
    sources = {2: np.array([[0, 0], [0, 60]]),
               5: np.array([[24 * ((3 * (i + j) ** 2 + 1) % 5) for j in range(5)] for i in range(5)]),
               7: (20 * seven + 104) % 120,
               35: np.asarray(source_thirty_five(), dtype=np.int64)}
    field = CyclotomicField(120)
    imag, omega = field.zeta() ** 30, field.zeta() ** 40
    scalars = {2: 1 + imag, 5: 1 + 2 * imag, 7: 3 + 2 * omega,
               35: (1 + 2 * imag) * (3 + 2 * omega)}
    reports = []
    for order, source in sources.items():
        gram = exact_column_gram(source, 120)
        expected = np.zeros_like(gram)
        expected[np.arange(order), np.arange(order), 0] = order
        if (not np.array_equal(gram, expected)
                or scalars[order] * scalars[order].conjugate() != field.element([order])):
            raise ValueError("a normalized smaller-Butson metric input failed")
        reports.append({"order": order, "phase": 120, "complete_Gram_verified": True,
                        "integral_normalizer_norm_equals_order": True,
                        "normalized_source_is_unitary_over_the_required_field": True})
    return reports


def root_pair_control(input_order, descriptor=None):
    """Two actual orthogonal admitted-root rows of length1499, not a source."""
    problem = boundary_problem(input_order, descriptor)
    s, delta = problem["selected_rows"], problem["selected_metric_norm"]
    if delta > s:
        raise ValueError("the simple unit/zero root-sum row does not fit this branch")
    x = np.full(s, 24, dtype=np.int64)
    y = np.r_[np.full(delta, 64, dtype=np.int64), np.full(s - delta, 84, dtype=np.int64)]
    boundary = problem["T_L"][0]
    rows = np.array([np.r_[x, y, boundary], np.r_[y, x, boundary]])
    roots = root_vectors(120)
    gram = np.array([[roots[(first - second) % 120].sum(axis=0) for second in rows] for first in rows])
    expected = np.zeros_like(gram)
    expected[0, 0, 0] = expected[1, 1, 0] = 1499
    if not np.array_equal(gram, expected) or np.any(rows % 3 == 2) or np.any(rows % 5 == 0):
        raise ValueError("the actual-length admitted root pair failed its full exact Gram")
    return {"input_order": input_order, "source_order": 1499, "selected_rows": s,
            "unit_root_sum_positions": delta, "zero_root_sum_positions": s - delta,
            "E_row_squared_norm": delta, "F_row_squared_norm": 3 * delta + 4 * (s - delta),
            "pair_phases": {"unit_sum": [24, 64], "zero_sum": [24, 84]},
            "full_two_by_two_Gram_of_length1499_rows_verified": True,
            "all_entries_in_raw_phase120_lift_alphabet": True,
            "boundary_input_row_index": problem["selected_row_indices"][0],
            "two_row_exponent_sha256": hashlib.sha256(rows.astype("<i8").tobytes()).hexdigest(),
            "cross_row_completion_and_corner_still_missing": True}


def branch_prefilter(source_order, input_order, selection, metric_descriptor):
    parameters = branch_parameters(source_order, input_order)
    rows = row_selection(input_order, selection, source_order)
    metric = dict(metric_descriptor)
    kind, delta = metric.get("kind"), parameters["selected_metric_norm"]
    if kind in ("scalar", "monomial", "nonuniform_spectral_monomial"):
        status = "EXCLUDED_ROOT_SUM_ENTRY_BOUND" if delta > 4 else "OPEN_ROOT_COMPLETION"
        reason = "E has a nonzero entry with squared modulus delta, but every sum of two unit roots has squared modulus at most4."
    elif kind == "normalized_small_Butson_blocks":
        block = int(metric["block_order"])
        if parameters["selected_rows"] % block:
            raise ValueError("the proposed complete blocks do not tile the selected space")
        if delta > 4 * block:
            status, reason = "EXCLUDED_ROOT_SUM_ENTRY_BOUND", "Each E entry has squared modulus delta/block_order>4."
        elif delta % block:
            status, reason = "EXCLUDED_ROOT_SUM_INTEGRAL_NORM", "The rational squared modulus delta/block_order is not an integer, whereas the squared norm of an integral root sum is an algebraic integer."
        else:
            status, reason = "OPEN_ROOT_COMPLETION", "This entrywise gate does not settle the simultaneous corner and root split."
    elif kind == "dense_integral_root_split":
        status, reason = "OPEN_ROOT_COMPLETION", "The actual two-row root gate passes; a complete coupled E,F,D remains missing."
    else:
        raise ValueError("unsupported finite metric descriptor")
    return {**parameters, "row_selection_descriptor": dict(selection),
            "selected_row_indices_sha256": hashlib.sha256(json.dumps(rows).encode()).hexdigest(),
            "metric_descriptor": metric, "status": status, "reason": reason,
            "actual_row_subset_selected": True, "source_constructed": False}


def initial_catalogue(input_orders=(567, 729)):
    result = []
    selections = [{"kind": "stride", "step": step, "shift": shift}
                  for step, shift in ((1, 0), (1, 1), (2, 7))]
    for d in input_orders:
        block_orders = {567: (2,), 585: (), 729: (5, 7, 35)}[d]
        metrics = [{"kind": "monomial", "root_phase_shift": 0},
                   {"kind": "nonuniform_spectral_monomial", "phase_slope": 1, "phase_shift": 15},
                   {"kind": "nonuniform_spectral_monomial", "phase_slope": 7, "phase_shift": 24},
                   *[{"kind": "normalized_small_Butson_blocks", "block_order": k}
                     for k in block_orders],
                   {"kind": "dense_integral_root_split"}]
        for selection in selections:
            for metric in metrics:
                result.append(branch_prefilter(1499, d, selection, metric))
    return result


def certificate(input_orders=(567, 729)):
    catalogue = initial_catalogue(input_orders)
    return {"schema": "rectangular-smaller-Butson-sources-v1", "source_order": 1499, "phase": 120,
            "verified_small_metric_inputs": small_metric_inputs(),
            "branches": [{"parameters": branch_parameters(1499, d), "known_input": source_certificate(d),
                          "metric_norms": metric_certificate(d), "actual_root_pair": root_pair_control(d)}
                         for d in input_orders],
            "finite_branch_catalogue": catalogue,
            "finite_branch_count": len(catalogue),
            "excluded_finite_branches": sum(row["status"].startswith("EXCLUDED") for row in catalogue),
            "open_dense_completion_branches": sum(row["status"] == "OPEN_ROOT_COMPLETION" for row in catalogue),
            "new_characteristic_zero_source_or_real_order_constructed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--output", type=Path)
    action.add_argument("--check", type=Path)
    parser.add_argument("--input-order", type=int, choices=(567, 585, 729))
    args = parser.parse_args()
    report = certificate((args.input_order,) if args.input_order is not None else (567, 729))
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    elif json.loads(args.check.read_text()) != report:
        raise ValueError("the rectangular smaller-source certificate failed replay")
    print(json.dumps({key: report[key] for key in ("finite_branch_count", "excluded_finite_branches",
                                                   "open_dense_completion_branches")}))


if __name__ == "__main__":
    main()
