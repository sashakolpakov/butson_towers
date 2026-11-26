"""Additional exact Spence sources, without changing the original five families.

Six prime Paley parameters have literal strong starters among the retained
Williamson inputs. Their recurrence is classical. This module certifies
the sources and scalar endpoints; it does not infer finite-order novelty.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from functools import cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from experiments.scalar_lift import williamson_spence_family as original
from experiments.scalar_lift.williamson_scalar_endpoints import (
    ImplicitScalarEndpoint, group_algebra_description, verified_t_design,
)

ROOT = Path(__file__).resolve().parents[2]
INPUTS = ROOT / "experiments/scalar_lift/williamson_inputs.json"
NEAR_INPUTS = ROOT / "experiments/scalar_lift/near_williamson_inputs.json"
PARAMETERS = {
    5: (3, (1, -1, 1, 1)),
    53: (27, (1, -1, 1, 1)),
    97: (49, (1, 1, 1, -1)),
    109: (55, (1, 1, 1, 1)),
    113: (57, (1, 1, 1, 1)),
    137: (69, (1, -1, 1, 1)),
}
HIGHER_ORIGINAL_MEMBERS = ((13, 2), (13, 3), (61, 2), (17, 2))
require = original.require
integer = original.integer


def sha(data):
    return hashlib.sha256(data).hexdigest()


@cache
def starter(q):
    q = integer(q, "q")
    require(q in PARAMETERS, "unsupported additional Spence parameter")
    width, signs = PARAMETERS[q]
    path = NEAR_INPUTS if width == 69 else INPUTS
    data = json.loads(path.read_text())
    record = next(x for x in data["entries"] if x.get("width", x.get("order")) == width)
    require(not record.get("first_matrix_reversed", False),
            "the strong starter requires ordinary symmetric circulants")
    words = record["first_rows_pm"]
    require(len(words) == 4 and all(len(w) == width and set(w) <= {"+", "-"} for w in words),
            "invalid literal starter sign words")
    raw = np.array([[1 if x == "+" else -1 for x in word] for word in words], dtype=np.int64)
    raw_hash = sha(raw.astype(np.int8).tobytes())
    expected = (record["verification"]["literal_rows_int8_row_major_sha256"]
                if width == 69 else record["first_rows_int8_row_major_sha256"])
    require(raw_hash == expected, "literal starter source hash failed")
    rows = raw * np.array(signs, dtype=np.int64)[:, None]
    matrices = [original.circulant(row) for row in rows]
    eye = np.eye(width, dtype=np.int64)
    require(all(np.array_equal(a, a.T) for a in matrices), "starter symmetry failed")
    require(all(np.array_equal(a @ b, b @ a) for a in matrices for b in matrices),
            "starter commutation failed")
    d, e, f, g = matrices
    require(np.array_equal((d-e) @ (d-e) + (f-g) @ (f-g), 4*eye),
            "strong minus starter Gram failed")
    require(np.array_equal((d+e) @ (d+e) + (f+g) @ (f+g), 4*q*eye),
            "strong plus starter Gram failed")
    require(np.array_equal(sum(a @ a for a in matrices), 4*width*eye),
            "summed starter Gram failed")
    proof = {
        "q": q, "order": width, "channel_permutation": [0, 1, 2, 3],
        "channel_signs": list(signs), "input_file": str(path.relative_to(ROOT)),
        "input_file_sha256": sha(path.read_bytes()), "input_rows_int8_sha256": raw_hash,
        "literal_first_rows_pm": ["".join("+" if x == 1 else "-" for x in row) for row in rows],
        "all_signs_symmetry_and_commutation_verified": True,
        "complete_strong_minus_Gram_scale": 4,
        "complete_strong_plus_Gram_scale": 4*q,
        "complete_summed_Gram_scale": 4*width,
        "primary_source": record["source"] if width == 69 else data["source_reference"],
        "primary_literal_selection": record.get("literal_selection"),
    }
    rows.setflags(write=False)
    return rows, proof


@cache
def paley(q):
    q = integer(q, "q")
    require(q in PARAMETERS and q % 4 == 1
            and all(q % d for d in range(2, math.isqrt(q)+1)),
            "a supported prime Paley parameter is required")
    row = np.array([0] + [1 if pow(x, (q-1)//2, q) == 1 else -1
                          for x in range(1, q)], dtype=np.int64)
    matrix = original.circulant(row)
    require(np.array_equal(matrix, matrix.T) and not np.any(matrix.sum(axis=1))
            and np.array_equal(matrix @ matrix, q*np.eye(q, dtype=np.int64)
                               - np.ones((q, q), dtype=np.int64)),
            "complete prime Paley identities failed")
    row.setflags(write=False)
    return row


@dataclass(frozen=True)
class AdditionalSpenceSource(original.SpenceSource):
    def __post_init__(self):
        q, r = integer(self.q, "q"), integer(self.r, "r")
        require(q in PARAMETERS and r >= 0, "require an additional supported q and r>=0")
        object.__setattr__(self, "q", q)
        object.__setattr__(self, "r", r)
        starter(q)
        paley(q)

    @property
    def group_dimensions(self):
        return ((self.q+1)//2,) + (self.q,)*self.r

    def entries(self, i, j):
        i, j = integer(i, "row"), integer(j, "column")
        require(0 <= i < self.order and 0 <= j < self.order, "source index out of range")
        labels = []
        for _ in range(self.r):
            i, a = divmod(i, self.q)
            j, b = divmod(j, self.q)
            labels.append(((b-a) % self.q, a == b))
        rows, _ = starter(self.q)
        values = tuple(map(int, rows[:, (j-i) % rows.shape[1]]))
        for difference, diagonal in reversed(labels):
            chi = int(paley(self.q)[difference])
            out = []
            for a, b in (values[:2], values[2:]):
                u, v = (a-b)//2, (a+b)//2
                out.extend((u+v*(chi+diagonal), u+v*(chi-diagonal)))
            values = tuple(out)
        require(all(x in (-1, 1) for x in values), "recursive source entry left sign alphabet")
        return values


def source(q, r=1):
    q, r = integer(q, "q"), integer(r, "r")
    cls = AdditionalSpenceSource if q in PARAMETERS else original.SpenceSource
    return cls(q, r)


def endpoint(q, r=1):
    """Literal entry access at arbitrary r, with the proved source recurrence."""
    matrices, _ = verified_t_design(631)
    return ImplicitScalarEndpoint(matrices, source(q, r))


def finite_member(q, r=1):
    """Complete finite group-norm replay and allocation-free endpoint proof.

    Full source rows are materialized only within the exact NTT size bound.
    Endpoint controls use individual entries, never a full expanded row.
    """
    s = source(q, r)
    volume = math.prod(2*n-1 for n in s.group_dimensions)
    require(8*s.order < original.MODULUS and volume <= 1 << 23,
            "finite member exceeds the exact NTT replay bound; use implicit entry access")
    rows = s.coefficient_rows()
    norms = original.exact_group_norms(rows, s.group_dimensions, s.q, s.r)
    selected = sorted({0, 1, s.order//2, s.order-1})
    for i in selected:
        for j in selected:
            left, right = np.unravel_index(i, s.group_dimensions), np.unravel_index(j, s.group_dimensions)
            difference = 0
            for a, b, n in zip(left, right, s.group_dimensions, strict=True):
                difference = difference*n + int((b-a) % n)
            require(s.entries(i, j) == tuple(map(int, rows[:, difference]))
                    and s.entries(i, j) == s.entries(j, i), "literal group-development control failed")
    matrices, t_proof = verified_t_design(631)
    h = ImplicitScalarEndpoint(matrices, s)
    endpoint_indices = sorted({0, 1, s.order, h.order//2, h.order-1})
    entries = [[h.entry(i, j) for j in endpoint_indices] for i in endpoint_indices]
    require(all(x in (-1, 1) for row in entries for x in row), "literal endpoint sign control failed")
    m, p = 4, 4*s.q
    trace = []
    for level in range(s.r):
        require(s.q*m-p == 0, "Spence recurrence mixed coefficient failed")
        trace.append({"r": level, "minus_Gram_scale": m, "plus_Gram_scale": p,
                      "mixed_J_coefficient": 0, "next_minus_Gram_scale": p,
                      "next_plus_Gram_scale": s.q*p})
        m, p = p, s.q*p
    require((m+p)//2 == 4*s.order, "source recurrence Gram scale failed")
    return {
        "q": s.q, "r": s.r, "source_order": s.order,
        "group_dimensions": list(s.group_dimensions),
        "strong_starter": (starter(s.q)[1] if s.q in PARAMETERS else original.starter(s.q)[1]),
        "source_recurrence_reference": original.PRIMARY,
        "literal_first_rows_int8_row_major_sha256": sha(rows.astype(np.int8).tobytes()),
        "actual_source_entry_alphabet": [-1, 1], "complete_exact_source_norm_checks": norms,
        "recurrence_trace": trace, "source_algebra": group_algebra_description(s.group_dimensions),
        "source_entry_control_indices": selected,
        "weighted_OD_source_order_over_algebra": 2524,
        "weighted_OD_source_Gram_scale": h.order,
        "T631_orthogonal_design_proof": json.loads(json.dumps(t_proof)),
        "endpoint_order": h.order, "endpoint_Gram_scale": h.order,
        "endpoint_entry_recipe": "sum_k X[k,i,j]*A[k,a,b], row=(i,a), column=(j,b)",
        "endpoint_entry_control_indices": endpoint_indices, "endpoint_entry_control_values": entries,
        "full_endpoint_Gram_verified_algebraically": True,
        "endpoint_Gram_proof": "Complete T631 polynomial OD identity, exact commuting self-adjoint group-source norm, and disjoint sign supports.",
        "full_dense_source_or_endpoint_allocation_performed": s.r == 0,
        "allocation_scope": "Only step zero uses dense starter components; all positive-step source norms use group rows and no dense endpoint is allocated.",
        "finite_order_novelty_inferred": False,
    }


def certificate():
    return {
        "schema": "additional-spence-scalar-families-v1",
        "implementation_sha256": sha(Path(__file__).read_bytes()),
        "original_recurrence_implementation_sha256": sha(Path(original.__file__).read_bytes()),
        "scope": "Six newly implemented classical source recurrences and exact finite scalar endpoints; finite-order novelty is classified separately.",
        "construction_reference": original.PRIMARY,
        "construction_reference_location": "Spence (1977), Theorem 4, printed p.254",
        "original_five_parameter_registry_unchanged": sorted(original.PARAMETERS),
        "families": [{"q": q, "base_source_order": PARAMETERS[q][0],
                      "endpoint_order_formula": f"{2524*PARAMETERS[q][0]}*{q}^r",
                      "parameter_range": "all integers r>=0", "literal_entry_recipe_complete": True}
                     for q in PARAMETERS],
        "strong_starters": [starter(q)[1] for q in PARAMETERS],
        "Paley_checks": [{"q": q, "prime": q, "first_row_int8_sha256": sha(paley(q).astype(np.int8).tobytes()),
                          "complete_integer_symmetry_zero_sum_and_square_identities_verified": True}
                         for q in PARAMETERS],
        "first_expanded_members": [finite_member(q) for q in PARAMETERS],
        "higher_original_family_members": [finite_member(q, r) for q, r in HIGHER_ORIGINAL_MEMBERS],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--q", type=int)
    parser.add_argument("--r", type=int, default=1)
    parser.add_argument("--entry", nargs=3, type=int, metavar=("CHANNEL", "ROW", "COLUMN"))
    args = parser.parse_args()
    if args.entry is not None:
        require(args.q is not None, "--entry requires --q")
        print(source(args.q, args.r).entry(*args.entry))
        return
    result = certificate() if args.q is None else finite_member(args.q, args.r)
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
