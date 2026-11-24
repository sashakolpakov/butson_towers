"""Small exact reference executor for norm-preserving scalar descent.

This transports the complete operator and its Hermitian form.  It never
optimizes a Gram residual and never labels a weighted child a Butson matrix.
Large sources should retain the same state as a lazy block DAG instead of
materializing the reference implementation's dense coefficient arrays.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass

import numpy as np
from sympy import Matrix, eye, kronecker_product, totient

from experiments.scalar_lift.relative_scalar_restriction import (
    RationalArray,
    conjugate_transpose,
    edge,
    regular_matrix,
)
from experiments.scalar_lift.vanishing_collision_search import integral_roots


def _from_regular(matrix, phase, rows, columns):
    h = int(totient(phase))
    coefficients = [[[matrix[i*h+a, j*h] for a in range(h)]
                     for j in range(columns)] for i in range(rows)]
    result = RationalArray.exact(coefficients)
    if regular_matrix(result, phase) != matrix:
        raise ValueError("matrix is not a regular lower-field operator")
    return result


def _restrict_array(value, restriction):
    coefficients = value.to_sympy()
    actions = restriction.actions.to_sympy()[:restriction.upper_degree]
    expanded = np.einsum("rce,eabu->racbu", coefficients, actions)
    rows, columns = coefficients.shape[:2]
    return RationalArray.exact(expanded.reshape(rows*restriction.width,
                                                columns*restriction.width,
                                                restriction.lower_degree))


@dataclass(frozen=True)
class ExactNormState:
    phase: int
    operator: RationalArray
    metric: RationalArray
    source_order: int
    relative_degree: int
    parent: str | None
    provenance: dict

    @property
    def order(self):
        return self.operator.numerators.shape[0]

    def definition(self):
        return {"schema": "exact-norm-descent-v1", "phase": self.phase,
                "operator": self.operator.to_json(), "metric": self.metric.to_json(),
                "source_order": self.source_order, "relative_degree": self.relative_degree,
                "parent": self.parent, "provenance": self.provenance}

    @property
    def digest(self):
        text = json.dumps(self.definition(), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(text.encode()).hexdigest()

    def verify(self):
        m = regular_matrix(self.operator, self.phase)
        adjoint = regular_matrix(conjugate_transpose(self.operator, self.phase), self.phase)
        g = regular_matrix(self.metric, self.phase)
        if regular_matrix(conjugate_transpose(self.metric, self.phase), self.phase) != g:
            raise ValueError("transported metric is not Hermitian")
        if adjoint*g*m != self.source_order*g:
            raise ValueError("complete transported norm identity failed")
        return True

    def restrict(self, lower_phase, basis_exponents=None):
        restriction = edge(self.phase, lower_phase, basis_exponents)
        m = _restrict_array(self.operator, restriction)
        restricted_g = _restrict_array(self.metric, restriction)
        local_g = regular_matrix(restriction.metric, lower_phase)
        # Transfer of a NONIDENTITY incoming form; trace transitivity is
        # essential.  A naive Kronecker product would discard its entries.
        transferred = kronecker_product(eye(self.order), local_g) * regular_matrix(restricted_g, lower_phase)
        g = _from_regular(transferred, lower_phase, m.numerators.shape[0], m.numerators.shape[0])
        child = ExactNormState(int(lower_phase), m, g, self.source_order,
                               self.relative_degree*restriction.width, self.digest,
                               {"operation": "restriction_of_scalars",
                                "edge_sha256": restriction.metadata()["exact_definition_sha256"],
                                "upper_phase": self.phase, "lower_phase": int(lower_phase),
                                "relative_basis_exponents": list(restriction.basis_exponents),
                                "metric_operation": "Tr_upper/lower of incoming Hermitian form"})
        child.verify()
        return child

    def endpoint(self, input_basis, output_basis, amplification=None):
        """Audit two supplied rational lattices; solve no norm equations.

        Ordinary target Gram follows iff the transported forms satisfy
        G_out=c I and G_in=amplification*c I.  Flatness is checked separately.
        """
        if self.phase not in (1, 2):
            raise ValueError("a rational endpoint is required")
        s, t = Matrix(input_basis), Matrix(output_basis)
        if s.shape != (self.order, self.order) or t.shape != s.shape:
            raise ValueError("two square bases of the complete endpoint are required")
        if not s.det() or not t.det():
            raise ValueError("endpoint bases must be nonsingular")
        m = regular_matrix(self.operator, self.phase)
        g = regular_matrix(self.metric, self.phase)
        b = t.inv()*m*s
        g_in, g_out = s.T*g*s, t.T*g*t
        if b.T*g_out*b != self.source_order*g_in:
            raise ValueError("two-basis norm transport failed")
        amplification = self.relative_degree if amplification is None else int(amplification)
        c = g_out[0, 0]
        normalized = c > 0 and g_out == c*eye(self.order) and g_in == amplification*c*eye(self.order)
        if normalized and b.T*b != self.source_order*amplification*eye(self.order):
            raise ValueError("ordinary normalized norm failed")
        defects = sum(value not in (-1, 1) for value in b)
        return {"state_sha256": self.digest, "dimension": self.order,
                "source_norm_scale": self.source_order, "amplification": amplification,
                "weighted_norm_verified": True, "Euclidean_norm_verified": bool(normalized),
                "Euclidean_norm_scale": self.source_order*amplification if normalized else None,
                "flat_entry_defects": defects, "sign_Hadamard_verified": bool(normalized and defects == 0),
                "matrix": [[str(value) for value in row] for row in b.tolist()],
                "input_metric": [[str(value) for value in row] for row in g_in.tolist()],
                "output_metric": [[str(value) for value in row] for row in g_out.tolist()]}


def from_butson(exponents, phase, source_id="explicit"):
    """Verify the actual complete source before admitting the root node."""
    raw = np.asarray(exponents)
    if raw.ndim != 2 or raw.shape[0] != raw.shape[1] or np.any(raw != raw.astype(np.int64)):
        raise ValueError("a square exact exponent matrix is required")
    exponents = raw.astype(np.int64) % int(phase)
    roots = np.asarray(integral_roots(phase), dtype=object)
    n, h = len(exponents), int(totient(phase))
    metric = np.zeros((n, n, h), dtype=object)
    metric[np.arange(n), np.arange(n), 0] = 1
    state = ExactNormState(int(phase), RationalArray.exact(roots[exponents]),
                           RationalArray.exact(metric), n, 1, None,
                           {"operation": "verified_actual_Butson_source", "source_id": source_id,
                            "exponents": exponents.tolist()})
    state.verify()
    return state


def positive_control():
    """Actual Fourier BH(4,8), exact 8->4->2 descent, then known H16."""
    exponents = 2*np.outer(np.arange(4), np.arange(4)) % 8
    source = from_butson(exponents, 8, "Fourier4-embedded-in-eighth-roots")
    middle, direct = source.restrict(4), source.restrict(2)
    terminal = middle.restrict(2)
    local_h = Matrix([[1, 1], [1, -1]])
    frame = kronecker_product(local_h, local_h)
    report = terminal.endpoint(kronecker_product(eye(4), frame), eye(16))
    # For the tower basis (0,2,1,3), compare literal direct and tower forms.
    change = eye(4)[:, [0, 2, 1, 3]]
    total = kronecker_product(eye(4), change)
    md, mt = regular_matrix(direct.operator, 2), regular_matrix(terminal.operator, 2)
    gd, gt = regular_matrix(direct.metric, 2), regular_matrix(terminal.metric, 2)
    if mt != total.inv()*md*total or gt != total.T*gd*total:
        raise ValueError("tower and direct restriction do not match by exact basis change")
    return {"source": source.definition(), "middle": middle.definition(),
            "terminal": terminal.definition(), "endpoint": report,
            "direct_and_tower_conjugacy_verified": True,
            "new_unknown_order_constructed": False}


if __name__ == "__main__":
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check", type=Path)
    args = parser.parse_args()
    result = positive_control()
    if args.check and json.loads(args.check.read_text()) != result:
        raise ValueError("exact norm descent certificate differs")
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"known_order": result["endpoint"]["dimension"],
                      "sign_Hadamard_verified": result["endpoint"]["sign_Hadamard_verified"],
                      "all_stages_exact_weighted_norm": True,
                      "new_unknown_order_constructed": False}))
