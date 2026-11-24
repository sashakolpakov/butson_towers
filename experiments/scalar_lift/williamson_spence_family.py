"""Literal infinite Williamson families from Spence's exact recurrence.

The sources are group-circulant sign quadruples on C_((q+1)/2) x
(F_q,+)^r. No dense large source or endpoint is needed for entry access.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from functools import cache
import hashlib
import json
from math import isqrt
import operator
from pathlib import Path

import numpy as np
ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "experiments/scalar_lift/williamson_inputs.json"
if not INPUT.exists():
    INPUT = Path(__file__).resolve().with_name("williamson_inputs.json")
PRIMARY = "https://doi.org/10.1017/S1446788700020255"
PARAMETERS = {13: (13, 1), 17: (17, 1), 61: (61, 1), 101: (101, 1), 121: (11, 2)}
FRAME_WIDTH = 2524
MODULUS = 998244353


def divisors(n):
    small = [d for d in range(1, isqrt(n)+1) if n % d == 0]
    return sorted(set(small + [n//d for d in small]))


def totient(n):
    result, remaining, prime = n, n, 2
    while prime*prime <= remaining:
        if remaining % prime == 0:
            result -= result//prime
            while remaining % prime == 0:
                remaining //= prime
        prime += 1
    if remaining > 1:
        result -= result//remaining
    return result


def require(value, message):
    if not value:
        raise ValueError(message)


def integer(value, label):
    require(not isinstance(value, (bool, np.bool_)), f"Boolean {label}")
    try:
        return operator.index(value)
    except TypeError as exc:
        raise ValueError(f"integer {label} required") from exc


def sha(value):
    return hashlib.sha256(value).hexdigest()


def circulant(row):
    row = np.asarray(row, dtype=np.int64)
    i = np.arange(len(row))
    return row[(i[None, :] - i[:, None]) % len(row)]


def field_subtract(q, a, b):
    p, degree = PARAMETERS[q]
    if degree == 1:
        return (a - b) % p
    return ((a % p - b % p) % p) + p * ((a // p - b // p) % p)


def character(q, a):
    p, degree = PARAMETERS[q]
    if a == 0:
        return 0
    if degree == 1:
        norm = a
    else:
        # GF(121)=GF(11)[u]/(u^2-2); 2 is a nonsquare modulo11.
        norm = ((a % p) ** 2 - 2 * (a // p) ** 2) % p
    value = pow(int(norm), (p - 1) // 2, p)
    require(value in (1, p - 1), "quadratic character failed")
    return 1 if value == 1 else -1


@cache
def paley(q):
    require(q in PARAMETERS, "unsupported prime-power parameter")
    row = np.array([character(q, a) for a in range(q)], dtype=np.int64)
    Q = np.array([[row[field_subtract(q, j, i)] for j in range(q)]
                  for i in range(q)], dtype=np.int64)
    require(np.array_equal(Q, Q.T) and not np.any(Q.sum(axis=1))
            and np.array_equal(Q @ Q, q * np.eye(q, dtype=np.int64) - np.ones((q, q), dtype=np.int64)),
            "complete Paley matrix identities failed")
    row.setflags(write=False)
    return row


@cache
def starter(q):
    require(q in PARAMETERS, "unsupported prime-power parameter")
    w = (q + 1) // 2
    if q in (13, 17):
        # Independently checked symmetric-circulant conference halves.
        rr, ss = {
            13: ([0,-1,-1,1,1,-1,-1], [1,-1,1,1,1,1,-1]),
            17: ([0,-1,1,-1,-1,-1,-1,1,-1], [1,-1,-1,1,1,1,1,-1,-1]),
        }[q]
        R, S = circulant(rr), circulant(ss)
        require(np.array_equal(R @ R + S @ S, q * np.eye(w, dtype=np.int64)),
                "conference-half starter Gram failed")
        row_i = np.eye(w, dtype=np.int64)[0]
        rows = np.array([ss, ss, np.array(rr) + row_i, np.array(rr) - row_i])
        origin = {"recipe": "D=E=S,F=I+R,G=-I+R; symmetric circulants R^2+S^2=qI",
                  "R_first_row": rr, "S_first_row": ss,
                  "complete_conference_half_Gram_verified": True}
    else:
        record = next(x for x in json.loads(INPUT.read_text())["entries"] if x["order"] == w)
        raw = np.array([[1 if c == "+" else -1 for c in s]
                        for s in record["first_rows_pm"]], dtype=np.int64)
        indices, signs = {61: ((0,3,1,2), (1,-1,1,1)),
                          101: ((0,1,2,3), (1,1,1,1)),
                          121: ((0,1,2,3), (1,1,1,-1))}[q]
        rows = raw[list(indices)] * np.array(signs)[:, None]
        origin = {"recipe": "signed re-pairing of exact retained Williamson input",
                  "input": str(INPUT.relative_to(ROOT)), "input_sha256": sha(INPUT.read_bytes()),
                  "raw_order": w, "channel_indices": list(indices), "channel_signs": list(signs)}
    matrices = [circulant(a) for a in rows]
    eye = np.eye(w, dtype=np.int64)
    require(np.all(np.abs(rows) == 1) and all(np.array_equal(a, a.T) for a in matrices),
            "starter signs/symmetry failed")
    for i, a in enumerate(matrices):
        for b in matrices[:i]:
            require(np.array_equal(a @ b, b @ a), "starter commutativity failed")
    d, e, f, g = matrices
    require(np.array_equal((d-e) @ (d-e) + (f-g) @ (f-g), 4 * eye),
            "strong Spence minus starter identity failed")
    require(np.array_equal((d+e) @ (d+e) + (f+g) @ (f+g), 4*q * eye),
            "strong Spence plus starter identity failed")
    require(np.array_equal(sum(a @ a for a in matrices), 4*w * eye),
            "starter Williamson Gram failed")
    proof = {"q": q, "order": w, "origin": origin,
             "literal_first_rows_pm": ["".join("+" if v == 1 else "-" for v in a) for a in rows],
             "all_signs_symmetry_and_commutativity_verified": True,
             "complete_strong_minus_Gram": "4I", "complete_strong_plus_Gram": f"{4*q}I",
             "complete_summed_Gram": f"{4*w}I"}
    rows.setflags(write=False)
    return rows, proof


@dataclass(frozen=True)
class SpenceSource:
    q: int
    r: int

    def __post_init__(self):
        q, r = integer(self.q, "q"), integer(self.r, "r")
        require(q in PARAMETERS and r >= 0, "require supported q and r>=0")
        object.__setattr__(self, "q", q)
        object.__setattr__(self, "r", r)
        starter(q)
        paley(q)

    @property
    def order(self):
        return (self.q + 1) // 2 * self.q ** self.r

    @property
    def group_dimensions(self):
        p, degree = PARAMETERS[self.q]
        return ((self.q + 1) // 2,) + (p,) * (degree * self.r)

    def entries(self, i, j):
        i, j = integer(i, "row"), integer(j, "column")
        require(0 <= i < self.order and 0 <= j < self.order, "source index out of range")
        labels = []
        for _ in range(self.r):
            i, a = divmod(i, self.q)
            j, b = divmod(j, self.q)
            labels.append((field_subtract(self.q, b, a), a == b))
        rows, _ = starter(self.q)
        values = tuple(map(int, rows[:, (j-i) % rows.shape[1]]))
        for difference, diagonal in reversed(labels):
            chi = int(paley(self.q)[difference])
            out = []
            for a, b in (values[:2], values[2:]):
                u, v = (a-b)//2, (a+b)//2
                out.extend((u + v*(chi + diagonal), u + v*(chi - diagonal)))
            values = tuple(out)
        require(all(v in (-1, 1) for v in values), "recursive source entry left sign alphabet")
        return values

    def entry(self, channel, i, j):
        channel = integer(channel, "channel")
        require(0 <= channel < 4, "source channel out of range")
        return self.entries(i, j)[channel]

    def row(self, channel, i):
        channel, i = integer(channel, "channel"), integer(i, "row")
        require(0 <= channel < 4 and 0 <= i < self.order, "source channel/row out of range")
        return np.fromiter((self.entry(channel, i, j) for j in range(self.order)),
                           dtype=np.int8, count=self.order)

    def coefficient_rows(self):
        return np.array([self.entries(0, j) for j in range(self.order)], dtype=np.int64).T


@cache
def ntt_permutation(n):
    rev, value = np.zeros(n, dtype=np.int64), np.arange(n, dtype=np.int64)
    for _ in range(n.bit_length() - 1):
        rev = (rev << 1) | (value & 1)
        value >>= 1
    return rev


@cache
def ntt_roots(length, inverse):
    root = pow(3, (MODULUS-1)//length, MODULUS)
    if inverse:
        root = pow(root, MODULUS-2, MODULUS)
    values, value = [], 1
    for _ in range(length//2):
        values.append(value)
        value = value * root % MODULUS
    return np.array(values, dtype=np.int64)


def ntt(values, inverse=False):
    a = np.asarray(values, dtype=np.int64)[:, ntt_permutation(values.shape[1])].copy() % MODULUS
    length = 2
    while length <= a.shape[1]:
        v = a.reshape(a.shape[0], -1, length)
        u = v[:, :, :length//2].copy()
        z = v[:, :, length//2:] * ntt_roots(length, inverse) % MODULUS
        v[:, :, :length//2] = (u + z) % MODULUS
        v[:, :, length//2:] = (u - z) % MODULUS
        length *= 2
    if inverse:
        a = a * pow(a.shape[1], MODULUS-2, MODULUS) % MODULUS
    return a


def exact_group_norms(rows, dimensions, q, r):
    """Check every group-ring coefficient by exact NTT plus an integer bound."""
    size = int(np.prod(dimensions))
    targets = (4*size, 4*q**r, 4*q**(r+1))
    require(rows.shape == (4, size) and np.all(np.abs(rows) == 1),
            "exact convolution requires sign coefficient rows")
    require(8*size < MODULUS and all(0 <= a <= 8*size for a in targets),
            "unsupported exact convolution size/target bound")
    coordinates = np.array(np.unravel_index(np.arange(size), dimensions))
    reverse = np.zeros(size, dtype=np.int64)
    encoded = np.zeros(size, dtype=np.int64)
    coefficient_volume, ordinary_stride = 1, 1
    for axis in range(len(dimensions)-1, -1, -1):
        n = dimensions[axis]
        reverse += ((-coordinates[axis]) % n) * ordinary_stride
        encoded += coordinates[axis] * coefficient_volume
        ordinary_stride *= n
        coefficient_volume *= 2*n-1
    require(np.array_equal(rows, rows[:, reverse]), "group polynomials are not self-conjugate")
    length = 1 << (coefficient_volume-1).bit_length()
    require(length <= 1 << 23, "NTT size exceeds modulus two-power factor")
    padded = np.zeros((4, length), dtype=np.int64)
    padded[:, encoded] = rows
    transformed = ntt(padded)
    norm = sum(a*a % MODULUS for a in transformed) % MODULUS
    minus = sum(((transformed[i]-transformed[i+1]) % MODULUS) ** 2 % MODULUS
                for i in (0, 2)) % MODULUS
    plus = sum(((transformed[i]+transformed[i+1]) % MODULUS) ** 2 % MODULUS
               for i in (0, 2)) % MODULUS
    products = ntt(np.array([norm, minus, plus]), inverse=True)
    require(not np.any(products[:, coefficient_volume:]), "linear-convolution padding alias")
    folded = np.zeros((3, size), dtype=np.int64)
    indices = np.arange(coefficient_volume, dtype=np.int64)
    dest = np.zeros(coefficient_volume, dtype=np.int64)
    stride = 1
    for n in reversed(dimensions):
        dest += (indices % (2*n-1) % n) * stride
        indices //= 2*n-1
        stride *= n
    for i in range(3):
        np.add.at(folded[i], dest, products[i, :coefficient_volume])
    folded %= MODULUS
    expected = np.zeros_like(folded)
    expected[:, 0] = targets
    require(np.array_equal(folded, expected), "complete exact group-ring norm failed")
    return {"all_group_coefficients_per_identity_checked": size,
            "complete_identities_checked": ["sum_j A_j^2=4wI", "(D-E)^2+(F-G)^2=4q^rI", "(D+E)^2+(F+G)^2=4q^(r+1)I"],
            "method": "integer NTT of Kronecker-encoded multivariate polynomials, followed by exact group folding",
            "modulus": MODULUS, "NTT_length": length,
            "absolute_integer_coefficient_bound": 8*size,
            "bound_strictly_below_modulus": True,
            "all_polynomials_self_conjugate": True,
            "integer_exactness": "target residual coefficients have magnitude at most8w<modulus; zero residues therefore prove zero integers"}


def cyclotomic_factors(source):
    p, _ = PARAMETERS[source.q]
    multiplicity = (source.q**source.r - 1)//(p-1)
    out = [{"index": int(d), "field_degree": int(totient(d)), "multiplicity": 1}
           for d in divisors((source.q+1)//2)]
    if multiplicity:
        out += [{"index": int(p*d), "field_degree": int(totient(p*d)), "multiplicity": multiplicity}
                for d in divisors((source.q+1)//2)]
    require(sum(f["field_degree"]*f["multiplicity"] for f in out) == source.order,
            "cyclotomic dimensions do not sum to group order")
    return out


def finite_member(q):
    source = SpenceSource(q, 1)
    rows = source.coefficient_rows()
    norms = exact_group_norms(rows, source.group_dimensions, q, 1)
    selected = sorted({0, 1, source.order//2, source.order-1})
    for i in selected:
        for j in selected:
            require(source.entries(i, j) == source.entries(j, i), "literal symmetry control failed")
            require(all(source.entry(k, i, j) == int(source.row(k, i)[j]) for k in range(4)),
                    "literal entry/row control failed")
    dense = q in (13, 17)
    if dense:
        coords = np.array(np.unravel_index(np.arange(source.order), source.group_dimensions))
        difference = np.zeros((source.order, source.order), dtype=np.int64)
        for c, d in zip(coords, source.group_dimensions):
            difference = difference*d + (c[None, :]-c[:, None]) % d
        matrices = rows[:, difference]
        require(np.array_equal(sum(a @ a.T for a in matrices),
                               4*source.order*np.eye(source.order, dtype=np.int64)),
                "independent complete dense expanded-source Gram failed")
    return {"q": q, "r": 1, "source_order": source.order,
            "group_dimensions": list(source.group_dimensions),
            "literal_first_rows_int8_row_major_sha256": sha(rows.astype(np.int8).tobytes()),
            "actual_source_component_entry_alphabet": [-1, 1],
            "complete_exact_source_norm_checks": norms,
            "full_dense_summed_source_Gram_executed": dense,
            "dense_summed_Gram_entries_checked": source.order**2 if dense else 0,
            "literal_entry_and_row_indices_checked": selected,
            "cyclotomic_factors": cyclotomic_factors(source),
            "source_algebra_dimension": source.order,
            "weighted_OD_source_order_over_algebra": 2524,
            "weighted_OD_source_Gram_scale": FRAME_WIDTH*source.order,
            "literal_T631_endpoint_order_and_Gram_scale": FRAME_WIDTH*source.order,
            "endpoint_allocation_performed": False,
            "order_novelty_scope": "classified separately in the completed endpoint ledger; no family-wide novelty inference"}


def certificate():
    return {"schema": "constructive-spence-williamson-scalar-families-v1",
            "primary_reference": {"author": "Edward Spence", "title": "An infinite family of Williamson matrices",
                                  "year": 1977, "doi_url": PRIMARY, "location": "Theorem4, printedp254; all four recurrence lines checked against the PDF image"},
            "implementation": str(Path(__file__).resolve().relative_to(ROOT)),
            "implementation_sha256": sha(Path(__file__).read_bytes()),
            "starters": [starter(q)[1] for q in PARAMETERS],
            "Paley_parameters": [{"q": q, "prime": p, "degree": e,
                                  "field_model": f"GF({p})" if e == 1 else "GF(11)[u]/(u^2-2)",
                                  "literal_quadratic_character_first_row": paley(q).tolist(),
                                  "all_q_squared_Q_Gram_entries_verified": True,
                                  "identities": "Q=Q^T,Q^2=qI-J,QJ=JQ=0"}
                                 for q, (p, e) in PARAMETERS.items()],
            "recurrence": {"D_next": "(D-E)/2 tensor J + (D+E)/2 tensor(Q+I)",
                           "E_next": "(D-E)/2 tensor J + (D+E)/2 tensor(Q-I)",
                           "F_next": "(F-G)/2 tensor J + (F+G)/2 tensor(Q+I)",
                           "G_next": "(F-G)/2 tensor J + (F+G)/2 tensor(Q-I)"},
            "independent_induction_proof": {
                "M": "(D-E)^2+(F-G)^2", "P": "(D+E)^2+(F+G)^2",
                "M_next": "P tensor I", "P_next": "qP tensor I+(qM-P) tensor J",
                "invariants": "M_r=4q^rI,P_r=4q^(r+1)I; sumfour A_j^2=(M_r+P_r)/2=4w_rI",
                "symmetry_and_commutativity": "inherited since previous components commute and I,J,Q commute and are symmetric",
                "entry_test": "for paired signs a,b, exactly one of(a-b)/2,(a+b)/2 is nonzero; on the diagonal the Q term vanishes, off the diagonal the I term vanishes, giving only signs",
                "literal_entry_cost": "O(r) arithmetic and O(r) labels, no dense source matrix"},
            "scalar_restriction": {
                "algebra": "Q[C_w0 x (C_p)^(e*r)],q=p^e,w0=(q+1)/2",
                "cyclotomic_product": "product_{d|w0}Q(zeta_d) times product_{d|w0}Q(zeta_(p*d))^((q^r-1)/(p-1))",
                "dimension": "w_r=w0*q^r", "conjugation": "invert every group element",
                "trace": "normalized regular trace: coefficient of identity=(1/w_r)*sum of all cyclotomic field traces, with recorded multiplicities",
                "metric": "I_w in the inverse-group basis; factor power-basis blocks are RamanujanSum_d(j-i)/w_r",
                "complete_weighted_OD_source": "D_T(A_0,...,A_3), order4*631 over the group algebra and Gram4*631*w_rI",
                "actual_source_entry_alphabet": "+/-the four recorded group polynomials A_j; weighted, not asserted Butson",
                "literal_restriction_entry_images": "+/-the sign group-circulant matrices A_j",
                "endpoint_entry_alphabet": [-1, 1],
                "endpoint_order_and_Gram_scale": "2524*w_r",
                "zero_policy": "zero entries in individual cyclotomic factors are retained; the complete inverse-group basis reconstruction has only sign entries"},
            "families": [{"q": q, "source_order_formula": f"{(q+1)//2}*{q}^r",
                          "real_endpoint_order_formula": f"{FRAME_WIDTH*((q+1)//2)}*{q}^r",
                          "parameter_range": "all integers r>=0", "construction_complete": True,
                          "every_order_new_claimed": False} for q in PARAMETERS],
            "first_expanded_members": [finite_member(q) for q in PARAMETERS],
            "source_method_priority_claimed": False,
            "new_order_status_scope": "The constructive infinite families are established here; novelty of individual real orders is recorded separately."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--q", type=int, choices=PARAMETERS)
    parser.add_argument("--r", type=int, default=1)
    parser.add_argument("--entry", type=int, nargs=3, metavar=("CHANNEL", "ROW", "COLUMN"))
    args = parser.parse_args()
    if args.entry is not None:
        require(args.q is not None, "--entry requires --q")
        source = SpenceSource(args.q, args.r)
        print(source.entry(*args.entry))
        return
    report = certificate()
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"complete_family_parameters": list(PARAMETERS),
                      "first_expanded_source_orders": [x["source_order"] for x in report["first_expanded_members"]],
                      "first_expanded_endpoint_orders": [x["literal_T631_endpoint_order_and_Gram_scale"] for x in report["first_expanded_members"]]}))


if __name__ == "__main__":
    main()
