"""Exact finite-field rank comparison at the known Hadamard order 885924.

The actual scalar endpoint is inequivalent to the standard Paley I matrix.
The endpoint itself is not allocated: its rank follows from the certified
frame identities and a 702-by-702 dual-number source matrix. This is an
executable rank audit and mathematical comparison, not a Lean rank proof.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

try:
    from flint import nmod_mat
except ImportError:  # The exact rank-factor audit needs only NumPy.
    nmod_mat = None

from experiments.scalar_lift import hadamard_frame_search as frame
from experiments.scalar_lift.tripled_unreal_scalar_endpoints import ROOTS, components, verified_source

ROOT = Path(__file__).resolve().parents[2]
MAX_DIMENSION = 702


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def mod3_matrix(matrix):
    a = np.asarray(matrix)
    require(a.ndim == 2 and 0 < min(a.shape) and max(a.shape) <= MAX_DIMENSION
            and np.issubdtype(a.dtype, np.integer),
            "a nonempty integer matrix with dimensions at most 702 is required")
    return np.asarray(a % 3, dtype=np.int16)


def matrix_hash(matrix):
    return digest(mod3_matrix(matrix).astype(np.uint8).tobytes())


def dual_number_matrix(exponents):
    """Multiplication by a+b*zeta6, with zeta6=-1+epsilon, epsilon^2=0.

    With A=a-b and B=b modulo3, the basis(epsilon,1) gives [[A,B],[0,A]].
    """
    e = np.asarray(exponents)
    require(e.ndim == 2 and e.shape[0] == e.shape[1] and len(e) <= MAX_DIMENSION//2
            and np.issubdtype(e.dtype, np.integer) and np.all((0 <= e) & (e < 6)),
            "a square sixth-root exponent matrix of order at most351 is required")
    a, b = ROOTS[e, 0], ROOTS[e, 1]
    diagonal, upper = (a-b) % 3, b % 3
    return np.block([[diagonal, upper], [np.zeros_like(diagonal), diagonal]]).astype(np.int16)


def verify_rank_factorization(matrix, u, v, left_inverse, pivot_columns):
    """Verify both rank bounds without trusting the elimination algorithm."""
    d = mod3_matrix(matrix)
    factors = [np.asarray(x) for x in (u, v, left_inverse)]
    require(all(x.ndim == 2 and np.issubdtype(x.dtype, np.integer)
                and np.all((0 <= x) & (x < 3)) for x in factors),
            "rank-factor entries must be integer residues0,1,2")
    u, v, left_inverse = (np.asarray(x, dtype=np.int16) for x in factors)
    pivots = list(pivot_columns)
    rank = len(pivots)
    require(len(set(pivots)) == rank and all(type(i) is int and 0 <= i < d.shape[1] for i in pivots),
            "invalid pivot columns")
    require(u.shape == (d.shape[0], rank) and v.shape == (rank, d.shape[1])
            and left_inverse.shape == (rank, d.shape[0]), "rank-factor shape mismatch")
    require(np.array_equal(u, d[:, pivots]), "selected-column identity failed")
    # At most702 terms, each at most4: int16 products cannot overflow.
    require(np.array_equal(u @ v % 3, d), "rank factorization failed")
    require(np.array_equal(left_inverse @ u % 3, np.eye(rank, dtype=np.int16)),
            "left-inverse identity failed")
    return rank


def rank_factorization(matrix):
    """Construct an exact mod3 rank certificate using independent NumPy elimination."""
    d = mod3_matrix(matrix)
    a, elementary = d.copy(), np.eye(d.shape[0], dtype=np.int16)
    pivots = []
    rank = 0
    for column in range(d.shape[1]):
        nonzero = np.flatnonzero(a[rank:, column])
        if not len(nonzero):
            continue
        row = rank + int(nonzero[0])
        a[[rank, row]], elementary[[rank, row]] = a[[row, rank]], elementary[[row, rank]]
        if a[rank, column] == 2:
            a[rank] = 2*a[rank] % 3
            elementary[rank] = 2*elementary[rank] % 3
        factors = a[:, column].copy()
        factors[rank] = 0
        a = (a-factors[:, None]*a[rank]) % 3
        elementary = (elementary-factors[:, None]*elementary[rank]) % 3
        pivots.append(column)
        rank += 1
        if rank == d.shape[0]:
            break
    require(not np.any(a[rank:]), "elimination left nonzero residual rows")
    u, v, left_inverse = d[:, pivots], a[:rank], elementary[:rank]
    verify_rank_factorization(d, u, v, left_inverse, pivots)
    return u, v, left_inverse, pivots


def checked_rank(matrix):
    """Use optional FLINT, or replay the exact NumPy factorization."""
    a = mod3_matrix(matrix)
    if nmod_mat is None:
        return len(rank_factorization(a)[3])
    return int(nmod_mat(a.shape[0], a.shape[1], a.reshape(-1).tolist(), 3).rank())


def printed_pair_matrix(exponents):
    """Return [[A,B],[0,A]] for the unreal source before the global gauge.

    The admitted exponents are those of zeta6^(-1)*(A+sqrt(-3)B)/2.
    Its actual dual action is [[A,A+2B],[0,A]]=F*C*T*C modulo3,
    where F=[[I,I],[0,I]], C=diag(2I,I), and T is the returned matrix.
    """
    actual = dual_number_matrix(exponents)
    e = np.asarray(exponents, dtype=np.int64)
    require(set(map(int, np.unique(e))) <= {0, 1, 3, 4}, "the admitted frame palette is required")
    a, b = components((e+1) % 6)
    n = len(e)
    require(np.array_equal(actual[:n, :n], a % 3)
            and np.array_equal(actual[n:, n:], a % 3)
            and np.array_equal(actual[:n, n:], (a+2*b) % 3)
            and not np.any(actual[n:, :n]), "actual-source gauge binding failed")
    return np.block([[a % 3, b % 3], [np.zeros_like(a), a % 3]]).astype(np.int16), a, b


def small_frame_controls():
    """Direct checks on literal frames and endpoints, of order at most468."""
    t7 = np.zeros((4, 7), dtype=np.int64)
    t7[0, 0], t7[0, [1, 2, 4]] = 1, -1
    t7[1, 3] = t7[2, 6] = t7[3, 5] = 1
    e3 = np.array([[3, 1, 1], [1, 3, 1], [1, 1, 3]], dtype=np.int64)
    controls = []
    for name, t, orders in (("T1", np.array([[1], [0], [0], [0]]), (13, 39, 117)),
                            ("T7", t7, (3, 13))):
        a = frame.D4 @ t
        p, q = frame.gs(a-t).astype(np.int64), frame.gs(a+t).astype(np.int64)
        width = len(p)
        eye = np.eye(width, dtype=np.int64)
        require(np.array_equal(p @ p.T, width*eye) and np.array_equal(q @ q.T, width*eye)
                and np.array_equal(p @ q.T+q @ p.T, width*eye), "small exact frame Gram failed")
        nilpotent = (q @ p.T*pow(width, -1, 3)+eye) % 3
        require(not np.any(nilpotent @ nilpotent % 3)
                and checked_rank(nilpotent) == width//2
                and checked_rank(p+q) == width//2, "small frame module check failed")
        for order in orders:
            e = e3 if order == 3 else verified_source(order)[0]
            blocks = {0: p, 1: q, 3: -p, 4: -q}
            require(order*width <= MAX_DIMENSION, "dense comparison control exceeds dimension bound")
            endpoint = np.block([[blocks[int(x)] for x in row] for row in e])
            source_rank = checked_rank(dual_number_matrix(e))
            printed, _, _ = printed_pair_matrix(e)
            require(checked_rank(printed) == source_rank, "printed-pair gauge comparison failed")
            actual = checked_rank(endpoint)
            require(actual == width//2*source_rank, "direct endpoint rank disagrees with module formula")
            controls.append({"frame": name, "frame_width": width, "source_order": order,
                             "control_endpoint_order": order*width, "dual_source_rank": source_rank,
                             "direct_endpoint_rank": actual, "predicted_endpoint_rank": width//2*source_rank,
                             "all_integer_frame_Grams_and_nilpotent_checks_passed": True})
    return controls


def audit():
    frame_path = ROOT / "certificates/scalar-lift/hadamard-7572-frame-reuse.json"
    source_path = ROOT / "certificates/scalar-lift/tripled-unreal-scalar-endpoints.json"
    frame_proof = json.loads(frame_path.read_text())["frame_proof"]
    for flag in ("full_P_integer_Gram_verified", "full_Q_integer_Gram_verified",
                 "full_symmetric_cross_identity_verified"):
        require(frame_proof[flag] is True, "completed frame proof is missing")
    width = frame_proof["rational_dimension"]
    require(width == 2524 and width % 2 == 0 and math.gcd(width, 3) == 1,
            "frame dimension does not have the stated modular properties")
    source_record = json.loads(source_path.read_text())
    retained = source_record["complete_endpoint_recipes"][2]
    for key in ("rational_dimension", "P_packed_sha256", "Q_packed_sha256"):
        require(source_record["fixed_frame"][key] == frame_proof[key],
                "completed endpoint uses a different frame")
    require(source_record["actual_used_scalar_images"] == {"0": "P", "1": "Q", "3": "-P", "4": "-Q"},
            "completed endpoint uses a different scalar-image convention")
    e, source_proof = verified_source(351)
    require(source_proof["exponent_sha256"] == retained["source"]["exponent_sha256"]
            and retained["endpoint_order"] == 885924 and retained["endpoint_Gram_proved_exactly"],
            "source differs from the completed known-order endpoint")
    actual_dual = dual_number_matrix(e)
    d, printed_a, printed_b = printed_pair_matrix(e)
    u, v, left_inverse, pivots = rank_factorization(d)
    rank = verify_rank_factorization(d, u, v, left_inverse, pivots)
    require(checked_rank(d) == rank == checked_rank(actual_dual), "supplemental rank check disagrees")
    order = width*len(e)
    paley_parameter = order-1
    require(paley_parameter % 4 == 3 and all(paley_parameter % k for k in range(2, math.isqrt(paley_parameter)+1))
            and order % 3 == 0, "Paley comparison hypotheses failed")
    ours, paley = width//2*rank, order//2
    require(rank == 195 and ours == 246090 and paley == 442962 and ours != paley,
            "the stated rank separation was not established")
    return {
        "schema": "known-order-Paley-comparison-v1",
        "implementation_sha256": digest(Path(__file__).read_bytes()),
        "order": order, "order_status": "known", "field": "F3",
        "invariant": "Matrix rank modulo3; preserved by signed row and column permutations and transposition.",
        "source": {"order": len(e), "actual_admitted_exponents": sorted(map(int, np.unique(e))),
                   "exponent_int64_little_endian_sha256": source_proof["exponent_sha256"],
                   "completed_endpoint_certificate": str(source_path.relative_to(ROOT)),
                   "completed_endpoint_certificate_sha256": digest(source_path.read_bytes()),
                   "complete_exact_source_Gram_replayed": True},
        "frame": {"width": width, "certificate": str(frame_path.relative_to(ROOT)),
                  "certificate_sha256": digest(frame_path.read_bytes()),
                  "P_packed_sha256": frame_proof["P_packed_sha256"],
                  "Q_packed_sha256": frame_proof["Q_packed_sha256"],
                  "complete_prior_P_Q_cross_Gram_proof_reused": True,
                  "matrix_M": "M=P+Q; M*M^T=3dI", "determinant_3_valuation": width//2,
                  "rank_M_mod3": width//2,
                  "rank_proof": "M*M^T=0 mod3 gives rank<=d/2; 3^(d-rank) divides detM and v3(detM)=d/2 give rank>=d/2.",
                  "nilpotent": "N=QP^(-1)+I; N^2=0 and rank(N)=d/2",
                  "free_dual_number_module_multiplicity": width//2},
        "printed_pair_source": {"matrix_order": len(d), "basis": ["epsilon", "1"],
                        "ring": "F3[epsilon]/(epsilon^2), zeta6=-1+epsilon",
                        "matrix_recipe": "T=[[A3,B3],[0,A3]] modulo3, for the printed unreal pair U3=(A3+sqrt(-3)*B3)/2",
                        "actual_admitted_source": "zeta6^(-1)*U3",
                        "actual_dual_action_binding": "D_actual=[[A3,A3+2B3],[0,A3]]=F*C*T*C; F=[[I,I],[0,I]], C=diag(2I,I), both invertible",
                        "all_entries_of_actual_source_binding_checked": True,
                        "A3_sign_entries_int8_sha256": digest(printed_a.astype(np.int8).tobytes()),
                        "B3_sign_entries_int8_sha256": digest(printed_b.astype(np.int8).tobytes()),
                        "matrix_mod3_uint8_sha256": matrix_hash(d), "rank": rank,
                        "actual_dual_matrix_mod3_uint8_sha256": matrix_hash(actual_dual),
                        "pivot_columns": pivots, "selected_column_identity_verified": True,
                        "rank_factorization": {"recipe": "U=T[:,pivot_columns], T=UV, LU=I_rank",
                                               "U_dimensions": list(u.shape), "V_dimensions": list(v.shape),
                                               "L_dimensions": list(left_inverse.shape),
                                               "U_mod3_uint8_sha256": matrix_hash(u),
                                               "V_mod3_uint8_sha256": matrix_hash(v),
                                               "L_mod3_uint8_sha256": matrix_hash(left_inverse),
                                               "all_entries_of_both_product_identities_verified": True,
                                               "maximum_integer_product_coefficient": 4*MAX_DIMENSION},
                        "mandatory_exact_checks": ["NumPy mod3 Gaussian elimination", "selected-column identity T[:,pivot_columns]=U", "rank-factor identity UV=T", "left-inverse identity LU=I_rank"],
                        "optional_cross_check": "When python-flint is installed, its finite-field rank is also checked; the certificate replay does not require it."},
        "endpoint_rank": ours,
        "endpoint_rank_formula": "H=rho(source)*(I tensor P); P invertible mod3, and rho is1262 copies of the regular dual-number action.",
        "Paley_I": {"prime_parameter": paley_parameter, "rank_mod3": paley,
                    "primary_reference": "https://doi.org/10.1002/sapm1933121311",
                    "rank_proof": "K*K^T=0 and K+K^T=2I modulo3 imply rank(K)<=order/2 and rank(K)>=order/2."},
        "small_literal_frame_controls": small_frame_controls(),
        "largest_matrix_dimension_used": MAX_DIMENSION,
        "full_885924_endpoint_allocated": False,
        "conclusion": "The actual completed endpoint is inequivalent to the standard Paley I matrix of the same known order.",
        "new_Hadamard_order_claimed": False, "global_new_equivalence_class_claimed": False,
        "rank_comparison_kernel_checked_in_Lean": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path, help="replay and compare a retained JSON certificate")
    args = parser.parse_args()
    result = audit()
    if args.check:
        require(result == json.loads(args.check.read_text()), "retained rank certificate does not match the exact replay")
    text = json.dumps(result, indent=2)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    elif not args.check:
        print(text, end="")


if __name__ == "__main__":
    main()
