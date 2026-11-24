"""Independent exact trace/CRT and small orthogonal-design substitution audit.

This imports neither the frame-search verifier nor the endpoint executor.
The cyclic group algebra is a product of cyclotomic fields, not a field.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from sympy import Poly, Symbol, cyclotomic_poly, divisors, mobius, totient

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "experiments/scalar_lift/williamson_inputs.json"
if not INPUT.exists():
    INPUT = Path(__file__).resolve().with_name("williamson_inputs.json")
Z = Symbol("z")


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def circulant(row):
    row = np.asarray(row, dtype=np.int64)
    i = np.arange(len(row))
    return row[(i[None, :] - i[:, None]) % len(row)]


def gs(rows):
    """Literal Goethals--Seidel array, independently implemented."""
    a, b, c, d = map(circulant, rows)
    rev = np.eye(len(a), dtype=np.int64)[:, ::-1]
    return np.block([[a, b @ rev, c @ rev, d @ rev],
                     [-b @ rev, a, -d.T @ rev, c.T @ rev],
                     [-c @ rev, d.T @ rev, a, -b.T @ rev],
                     [-d @ rev, -c.T @ rev, b.T @ rev, a]])


def od_control():
    # Orbit representatives0,1,3; mu=2; states0,1,2.
    t = np.zeros((4, 7), dtype=np.int64)
    t[0, 0] = 1
    t[0, [1, 2, 4]] = -1
    t[1, 3], t[2, 6], t[3, 5] = 1, 1, 1
    for shift in range(7):
        require(int(np.sum(t * np.roll(t, shift, axis=1))) == (7 if shift == 0 else 0),
                "independent T7 periodic correlation failed")
    a, b, c, d = t
    rotations = (t, np.array([-b, a, d, -c]),
                 np.array([-c, -d, a, b]), np.array([-d, c, -b, a]))
    linear = [np.eye(4, dtype=np.int64),
              np.array([[0,-1,0,0],[1,0,0,0],[0,0,0,1],[0,0,-1,0]], dtype=np.int64),
              np.array([[0,0,-1,0],[0,0,0,-1],[1,0,0,0],[0,1,0,0]], dtype=np.int64),
              np.array([[0,0,0,-1],[0,0,1,0],[0,-1,0,0],[1,0,0,0]], dtype=np.int64)]
    for i in range(4):
        require(np.array_equal(linear[i] @ t, rotations[i]), "rotation formula mismatch")
        require(np.array_equal(linear[i].T @ linear[i], np.eye(4, dtype=np.int64)),
                "formal rotation diagonal coefficient failed")
        for j in range(i):
            require(not np.any(linear[i].T @ linear[j] + linear[j].T @ linear[i]),
                    "formal rotation mixed coefficient failed")
    x = [gs(v) for v in rotations]
    eye = np.eye(28, dtype=np.int64)
    require(np.all(sum(v * v for v in x) == 1), "OD supports are not disjoint and complete")
    for i in range(4):
        require(np.array_equal(x[i] @ x[i].T, 7 * eye), "OD diagonal coefficient failed")
        for j in range(i):
            require(not np.any(x[i] @ x[j].T + x[j] @ x[i].T),
                    "OD mixed coefficient failed")
    n = [v @ x[0].T for v in x]
    require(not np.array_equal(n[1] @ n[2], 7 * n[3])
            and not np.array_equal(n[1] @ n[2], -7 * n[3]),
            "Clifford/quaternion distinction control changed")
    volume_trace = int(np.trace(n[1] @ n[2] @ n[3])) // (7 ** 3)
    return x, {
        "p": 7, "multiplier": 2, "multiplier_order": 3, "states": [0, 1, 2],
        "literal_T_rows": t.tolist(), "all_7_periodic_correlations_verified": True,
        "rational_dimension": 28, "disjoint_complete_coefficient_supports": True,
        "coefficient_weight": 7, "all_10_OD_Gram_coefficients_verified": True,
        "real_4x4_rotation_coefficient_matrices": [a.tolist() for a in linear],
        "all_10_formal_rotation_norm_coefficients_verified": True,
        "brief_quaternion_scope_control": {
            "U_j": "X_j X_0^T/7",
            "U1_U2_equals_either_sign_of_U3": False,
            "trace_U1_U2_U3": volume_trace,
            "conclusion": "the OD identities alone do not identify this GS action with a Hamilton representation"}}


def companion(d):
    polynomial = Poly(cyclotomic_poly(d, Z), Z)
    degree = polynomial.degree()
    coefficients = [int(polynomial.nth(i)) for i in range(degree)]
    action = np.zeros((degree, degree), dtype=np.int64)
    if degree > 1:
        action[1:, :-1] = np.eye(degree - 1, dtype=np.int64)
    action[:, -1] = -np.array(coefficients, dtype=np.int64)
    return polynomial, action


def ramanujan(d, exponent):
    return sum(int(a * mobius(d // a)) for a in divisors(d) if exponent % a == 0)


def matrix_polynomial(coefficients, action):
    out = np.zeros_like(action)
    identity = np.eye(len(action), dtype=np.int64)
    for value in reversed(coefficients):
        out = out @ action + int(value) * identity
    return out


def audit_order(entry, x):
    w = entry["order"]
    rows = np.array([[1 if c == "+" else -1 for c in row]
                     for row in entry["first_rows_pm"]], dtype=np.int64)
    require(rows.shape == (4, w) and np.all(np.abs(rows) == 1), "invalid Williamson signs")
    williamson = [circulant(row) for row in rows]
    require(all(np.array_equal(a, a.T) for a in williamson), "Williamson symmetry failed")
    require(np.array_equal(sum(a @ a.T for a in williamson), 4 * w * np.eye(w, dtype=np.int64)),
            "full Williamson summed Gram failed")
    source_trace = np.zeros((w, w), dtype=np.int64)
    crt = np.zeros((w, w), dtype=np.int64)
    # beta_j=c^{-j}: rho(c)[r,s]=1 iff s-r=1 mod w.
    group_action = circulant([0, 1] + [0] * (w - 2))
    factors = []
    offset = 0
    for d in divisors(w):
        polynomial, action = companion(int(d))
        r = int(totient(d))
        require(np.array_equal(np.linalg.matrix_power(action, d), np.eye(r, dtype=np.int64)),
                "cyclotomic action order failed")
        inverse = np.linalg.matrix_power(action, d - 1)
        local_crt = np.column_stack([np.linalg.matrix_power(inverse, j)[:, 0] for j in range(w)])
        trace = np.array([[ramanujan(int(d), j - i) for j in range(r)]
                          for i in range(r)], dtype=np.int64)
        require(np.array_equal(action.T @ trace, trace @ inverse),
                "cyclotomic conjugation/trace compatibility failed")
        require(np.array_equal(local_crt @ group_action, action @ local_crt),
                "CRT multiplication action failed")
        source_trace[offset:offset+r, offset:offset+r] = trace
        crt[offset:offset+r] = local_crt
        f = [matrix_polynomial(row, action) for row in rows]
        require(all(np.array_equal(a, matrix_polynomial(row, inverse))
                    for a, row in zip(f, rows)), "factor polynomial is not self-conjugate")
        require(np.array_equal(sum(a @ a for a in f), 4 * w * np.eye(r, dtype=np.int64)),
                "factor-wise complete Williamson norm failed")
        for a, group_matrix in zip(f, williamson):
            require(np.array_equal(local_crt @ group_matrix, a @ local_crt),
                    "actual source polynomial action differs from its sign circulant")
        values = sorted({tuple(map(int, sign * a[:, 0])) for a in f for sign in (-1, 1)})
        factors.append({
            "cyclotomic_index": int(d), "field_degree": r,
            "field": f"Q(zeta_{d})", "cyclotomic_polynomial": str(polynomial.as_expr()),
            "basis": [f"zeta_{d}^{j}" for j in range(r)],
            "trace_metric_numerator": trace.tolist(), "trace_metric_denominator": w,
            "actual_source_entry_coordinates": [list(v) for v in values],
            "zero_is_in_actual_source_alphabet": [0] * r in [list(v) for v in values],
            "zero_preserved_by_the_restriction": True,
            "actual_restricted_scalar_entry_alphabet": sorted({int(v) for a in f for v in a.ravel()}
                                                              | {int(-v) for a in f for v in a.ravel()}),
            "all_used_source_actions_verified": True,
            "conjugation_is_trace_adjoint_verified": True,
            "complete_factor_norm_verified": f"sum_j f_j*conjugate(f_j)={4*w}",
            "module_rank": 28, "rational_dimension": 28 * r,
            "weighted_source_Gram_scale": 28 * w})
        offset += r
    require(offset == w, "CRT dimensions do not sum to w")
    require(np.array_equal(crt.T @ source_trace @ crt, w * np.eye(w, dtype=np.int64)),
            "normalized CRT metric does not recover the Euclidean group metric")
    endpoint = sum(np.kron(a, b) for a, b in zip(x, williamson))
    n = 28 * w
    require(np.all(np.abs(endpoint) == 1), "literal endpoint is not flat")
    require(np.array_equal(endpoint @ endpoint.T, n * np.eye(n, dtype=np.int64)),
            "complete independent small endpoint Gram failed")
    return {"williamson_order": w, "literal_first_rows_pm": entry["first_rows_pm"],
            "all_Williamson_symmetry_and_full_Gram_entries_verified": True,
            "CRT_factors": factors,
            "CRT_map_basis": "group basis beta_j=c^(-j) to factor power bases",
            "all_actual_polynomial_actions_verified": True,
            "S_transpose_G_CRT_S_equals_group_identity": True,
            "normalizing_denominator": w, "total_algebra_dimension": w,
            "weighted_source_order_over_group_algebra": 28,
            "source_and_endpoint_Gram_scale": n, "endpoint_order": n,
            "endpoint_alphabet": [-1, 1], "all_endpoint_Gram_entries_checked": n * n,
            "endpoint_int8_row_major_sha256": sha(endpoint.astype(np.int8).tobytes())}


def certificate():
    data = json.loads(INPUT.read_text())
    entries = {row["order"]: row for row in data["entries"]}
    x, od = od_control()
    return {
        "schema": "independent-williamson-cyclotomic-trace-audit-v1",
        "implementation": str(Path(__file__).resolve().relative_to(ROOT)),
        "implementation_sha256": sha(Path(__file__).read_bytes()),
        "literal_Williamson_input": str(INPUT.relative_to(ROOT)),
        "literal_Williamson_input_sha256": sha(INPUT.read_bytes()),
        "source_algebra": "Q[C_w]=product_{d|w} Q(zeta_d), not a field",
        "involution": "c -> c^(-1); factor-wise complex conjugation",
        "normalized_trace": "tau(a)=coefficient_1(a)=(1/w)*sum_{d|w}Tr_Q(zeta_d)/Q(a_d)",
        "group_basis": "1,c^(-1),...,c^(-(w-1))",
        "group_trace_metric": "I_w",
        "factor_trace_metrics": "G_d[i,j]=(1/w)*RamanujanSum_d(j-i)",
        "general_source": "D_T(f_0(c),...,f_3(c))=sum_j X_j*f_j(c)",
        "source_class": "complete weighted design over a commutative semisimple algebra; not a Butson source",
        "OD_formula": "D_T(x)=GS((t0*x0-t1*x1-t2*x2-t3*x3, t1*x0+t0*x1-t3*x2+t2*x3, t2*x0+t3*x1+t0*x2-t1*x3, t3*x0-t2*x1+t1*x2+t0*x3))",
        "general_OD_derivation": "The displayed real4x4 rotation L(x) has L(x)^T*L(x)=sum_j x_j^2*I. Hence sum_i V_i*V_i^T=sum_j x_j^2*sum_i T_i*T_i^T=p*sum_j x_j^2*I. The Goethals-Seidel identity for commuting real circulants gives GS(V)GS(V)^T=I4 tensor sum_i V_i*V_i^T. Coefficient comparison gives the ten OD identities for every periodic T tuple.",
        "general_Gram_proof": "XiXi^T=pI and XiXj^T+XjXi^T=0; self-conjugate commuting Williamson polynomials give D_T(f)D_T(f)^*=p*sum_j f_j*f_j^*=4pwI",
        "general_flatness_proof": "The Xi supports partition every matrix position, so every algebra-valued entry is +/-one f_j; rho(f_j) is its literal sign circulant.",
        "restriction_order_and_scale": "source order4p over Q[C_w], scalar dimensionw, endpoint order4pw and Gram scale4pw",
        "OD_small_control": od,
        "finite_independent_controls": [audit_order(entries[w], x) for w in (5, 9, 15, 17, 21)],
        "no_large_endpoint_allocated": True,
        "novelty_or_nonexistence_claimed_by_this_audit": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = certificate()
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"cyclotomic_trace_orders": [5, 9, 15, 17, 21],
                      "OD_Gram_coefficients_checked": 10,
                      "full_small_endpoint_orders": [140, 252, 420, 476, 588],
                      "all_exact_checks_passed": True}))


if __name__ == "__main__":
    main()
