"""Exact finite-phase certificates for complex generalized weighing matrices.

The first target is the omitted lifting calculation in Egan's published
nonexistence proof for ``CGW(11, 5; 4)``.  The support reduction and uniqueness
of the ``2-(11,5,2)`` design are literature inputs.  Everything after that
reduction is checked here by a small affine system over ``Z/4Z``.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from itertools import combinations
from typing import Any

from .coefficients import CyclotomicField

ExponentEntry = int | None
ExponentMatrix = tuple[tuple[ExponentEntry, ...], ...]
SupportRows = tuple[tuple[int, ...], ...]

CGW_11_5_4_SCHEMA = "hadamard.cgw-11-5-4-classification.v1"
CGW_EQUIVALENCE_CONTRACT = {
    "name": "classical CGW monomial equivalence",
    "definition": "W2 = P W1 Q* for independent monomial P,Q over mu_4",
    "row_permutations": True,
    "column_permutations": True,
    "row_phase_scalings": True,
    "column_phase_scalings": True,
    "transpose_quotiented": False,
    "coefficient_galois_quotiented": False,
    "entrywise_conjugation_quotiented": False,
}
CGW_11_5_4_SOURCE = {
    "author": "Ronan Egan",
    "title": (
        "A survey of complex generalized weighing matrices and a construction "
        "of quantum error-correcting codes"
    ),
    "result": "Proposition 2.15: there is no CGW(11,5;4)",
    "doi": "10.1016/j.disc.2024.114201",
    "arxiv": "https://arxiv.org/abs/2309.07522",
    "published": "Discrete Mathematics 348 (2025), 114201",
    "support_input": (
        "every hypothetical support is permutation-equivalent to the unique "
        "2-(11,5,2) design"
    ),
    "omitted_step_reproduced_here": (
        "the paper states that the unique support has no mu_4 lift but omits "
        "the lifting calculation"
    ),
    "support_reduction_machine_reproved": False,
}
CGW_11_5_4_PROOF_OBLIGATIONS = {
    "published_support_reduction_cited": True,
    "canonical_unique_support_verified_as_biplane": True,
    "every_supported_phase_variable_accounted": True,
    "every_row_pair_orthogonality_equation_encoded": True,
    "two_term_mu4_vanishing_reduction_exact": True,
    "affine_mod4_contradiction_verified": True,
    "classical_monomial_equivalence_accounted": True,
}

_ORDER = 11
_WEIGHT = 5
_PHASE = 4
_BIPLANE_DIFFERENCE_SET = (0, 1, 2, 4, 7)


def _canonical_json_digest(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def cyclic_biplane_support() -> SupportRows:
    """Return the cyclic representative of the unique ``2-(11,5,2)`` design."""

    return tuple(
        tuple(sorted((row + difference) % _ORDER for difference in _BIPLANE_DIFFERENCE_SET))
        for row in range(_ORDER)
    )


def support_matrix(support: SupportRows) -> tuple[tuple[int, ...], ...]:
    order = len(support)
    return tuple(
        tuple(int(column in support[row]) for column in range(order))
        for row in range(order)
    )


def _support_graph_connected(support: SupportRows) -> bool:
    order = len(support)
    columns_to_rows = {
        column: tuple(row for row in range(order) if column in support[row])
        for column in range(order)
    }
    pending = [("row", 0)]
    seen = set(pending)
    while pending:
        kind, index = pending.pop()
        neighbors = (
            (("column", column) for column in support[index])
            if kind == "row"
            else (("row", row) for row in columns_to_rows[index])
        )
        for neighbor in neighbors:
            if neighbor not in seen:
                seen.add(neighbor)
                pending.append(neighbor)
    return len(seen) == 2 * order


def audit_biplane_support(support: SupportRows | None = None) -> dict[str, Any]:
    """Check the complete incidence identities of the canonical support."""

    rows = cyclic_biplane_support() if support is None else support
    if len(rows) != _ORDER:
        raise ValueError("the biplane support must have eleven rows")
    if any(
        len(row) != _WEIGHT
        or len(set(row)) != _WEIGHT
        or any(column < 0 or column >= _ORDER for column in row)
        for row in rows
    ):
        raise ValueError("the biplane support has a malformed row")
    column_weights = tuple(
        sum(column in row for row in rows) for column in range(_ORDER)
    )
    row_intersections = tuple(
        len(set(rows[first]) & set(rows[second]))
        for first, second in combinations(range(_ORDER), 2)
    )
    columns = tuple(
        tuple(row for row in range(_ORDER) if column in rows[row])
        for column in range(_ORDER)
    )
    column_intersections = tuple(
        len(set(columns[first]) & set(columns[second]))
        for first, second in combinations(range(_ORDER), 2)
    )
    difference_multiplicities = {
        residue: sum(
            (left - right) % _ORDER == residue
            for left in _BIPLANE_DIFFERENCE_SET
            for right in _BIPLANE_DIFFERENCE_SET
            if left != right
        )
        for residue in range(1, _ORDER)
    }
    matrix = support_matrix(rows)
    gram = tuple(
        tuple(
            sum(matrix[row][column] * matrix[other][column] for column in range(_ORDER))
            for other in range(_ORDER)
        )
        for row in range(_ORDER)
    )
    expected_gram = tuple(
        tuple(_WEIGHT if row == other else 2 for other in range(_ORDER))
        for row in range(_ORDER)
    )
    valid = (
        column_weights == (_WEIGHT,) * _ORDER
        and row_intersections == (2,) * 55
        and column_intersections == (2,) * 55
        and set(difference_multiplicities.values()) == {2}
        and gram == expected_gram
        and _support_graph_connected(rows)
    )
    if not valid:
        raise ValueError("the supplied support is not the cyclic 2-(11,5,2) biplane")
    return {
        "order": _ORDER,
        "row_weight": _WEIGHT,
        "column_weight": _WEIGHT,
        "row_pair_count": len(row_intersections),
        "row_pair_intersection_values": sorted(set(row_intersections)),
        "column_pair_count": len(column_intersections),
        "column_pair_intersection_values": sorted(set(column_intersections)),
        "difference_multiplicities": {
            str(residue): multiplicity
            for residue, multiplicity in difference_multiplicities.items()
        },
        "incidence_gram_identity": "S S^T = 3 I_11 + 2 J_11",
        "connected_bipartite_support": True,
        "edge_count": _ORDER * _WEIGHT,
        "row_column_gauge_rank": 2 * _ORDER - 1,
        "gauge_fixed_phase_variable_count": (
            _ORDER * _WEIGHT - (2 * _ORDER - 1)
        ),
        "support_digest": _canonical_json_digest(matrix),
        "valid": True,
    }


def _validate_exponent_matrix(matrix: ExponentMatrix, phase: int) -> None:
    if phase <= 0 or not matrix or any(len(row) != len(matrix) for row in matrix):
        raise ValueError("a CGW exponent matrix must be nonempty and square")
    for row in matrix:
        for entry in row:
            if entry is not None and not isinstance(entry, int):
                raise ValueError("a CGW entry must be None or an integer exponent")


def is_cgw_exponent_matrix(
    matrix: ExponentMatrix,
    phase: int,
    weight: int,
) -> bool:
    """Check ``W W* = weight I`` exactly in ``Q(zeta_phase)``."""

    _validate_exponent_matrix(matrix, phase)
    order = len(matrix)
    if any(sum(entry is not None for entry in row) != weight for row in matrix):
        return False
    if any(
        sum(matrix[row][column] is not None for row in range(order)) != weight
        for column in range(order)
    ):
        return False
    field = CyclotomicField(phase)
    root = field.zeta()
    for first in range(order):
        for second in range(order):
            inner = field.zero()
            for column in range(order):
                left = matrix[first][column]
                right = matrix[second][column]
                if left is not None and right is not None:
                    inner += root ** ((left - right) % phase)
            expected = field.coerce(weight if first == second else 0)
            if inner != expected:
                return False
    return True


def apply_cgw_monomial_transform(
    matrix: ExponentMatrix,
    phase: int,
    *,
    row_permutation: tuple[int, ...],
    column_permutation: tuple[int, ...],
    row_phases: tuple[int, ...],
    column_phases: tuple[int, ...],
) -> ExponentMatrix:
    """Apply ``D_r P_r W P_c^T D_c*`` in exponent coordinates."""

    _validate_exponent_matrix(matrix, phase)
    order = len(matrix)
    if (
        sorted(row_permutation) != list(range(order))
        or sorted(column_permutation) != list(range(order))
        or len(row_phases) != order
        or len(column_phases) != order
    ):
        raise ValueError("invalid CGW monomial transformation data")
    out = []
    for target_row, source_row in enumerate(row_permutation):
        row = []
        for target_column, source_column in enumerate(column_permutation):
            entry = matrix[source_row][source_column]
            row.append(
                None
                if entry is None
                else (
                    entry
                    + row_phases[target_row]
                    - column_phases[target_column]
                )
                % phase
            )
        out.append(tuple(row))
    return tuple(out)


@dataclass(frozen=True, slots=True)
class BiplanePhaseConstraint:
    row_pair: tuple[int, int]
    overlap_columns: tuple[int, int]
    terms: tuple[tuple[int, int], ...]
    rhs: int = 2

    def vector(self, variable_count: int, modulus: int = _PHASE) -> tuple[int, ...]:
        row = [0] * variable_count
        for variable, coefficient in self.terms:
            row[variable] = (row[variable] + coefficient) % modulus
        return tuple(row)

    def to_json(self, edges: tuple[tuple[int, int], ...]) -> dict[str, Any]:
        return {
            "row_pair": list(self.row_pair),
            "overlap_columns": list(self.overlap_columns),
            "terms": [
                {
                    "variable_index": variable,
                    "edge": list(edges[variable]),
                    "coefficient": coefficient,
                }
                for variable, coefficient in self.terms
            ],
            "rhs_mod_4": self.rhs,
        }


def biplane_phase_system(
    support: SupportRows | None = None,
) -> tuple[tuple[tuple[int, int], ...], tuple[BiplanePhaseConstraint, ...]]:
    """Build every row-pair orthogonality congruence for a ``mu_4`` lift."""

    rows = cyclic_biplane_support() if support is None else support
    audit_biplane_support(rows)
    edges = tuple(
        (row, column)
        for row in range(_ORDER)
        for column in rows[row]
    )
    edge_index = {edge: index for index, edge in enumerate(edges)}
    constraints = []
    for first, second in combinations(range(_ORDER), 2):
        left_column, right_column = sorted(
            set(rows[first]) & set(rows[second])
        )
        constraints.append(
            BiplanePhaseConstraint(
                (first, second),
                (left_column, right_column),
                (
                    (edge_index[(first, left_column)], 1),
                    (edge_index[(second, left_column)], -1),
                    (edge_index[(first, right_column)], -1),
                    (edge_index[(second, right_column)], 1),
                ),
            )
        )
    return edges, tuple(constraints)


def find_modular_contradiction(
    constraints: tuple[BiplanePhaseConstraint, ...],
    variable_count: int,
    modulus: int = _PHASE,
) -> dict[str, Any]:
    """Find a row-combination witness ``0 = nonzero (mod modulus)``.

    The deterministic elimination uses only unit pivots.  For the biplane
    system it exposes a contradiction after rank 34.
    """

    if modulus != 4:
        raise ValueError("the current exact witness finder is specialized to modulus four")
    equation_count = len(constraints)
    rows = [
        [*constraint.vector(variable_count, modulus), constraint.rhs % modulus]
        for constraint in constraints
    ]
    combinations_matrix = [
        [int(row == column) for column in range(equation_count)]
        for row in range(equation_count)
    ]
    pivot_row = 0
    pivot_columns = []
    for column in range(variable_count):
        source = next(
            (
                row
                for row in range(pivot_row, equation_count)
                if rows[row][column] % 2
            ),
            None,
        )
        if source is None:
            continue
        rows[pivot_row], rows[source] = rows[source], rows[pivot_row]
        combinations_matrix[pivot_row], combinations_matrix[source] = (
            combinations_matrix[source],
            combinations_matrix[pivot_row],
        )
        inverse = pow(rows[pivot_row][column], -1, modulus)
        rows[pivot_row] = [
            inverse * value % modulus for value in rows[pivot_row]
        ]
        combinations_matrix[pivot_row] = [
            inverse * value % modulus
            for value in combinations_matrix[pivot_row]
        ]
        for row in range(equation_count):
            if row == pivot_row:
                continue
            multiplier = rows[row][column]
            if not multiplier:
                continue
            rows[row] = [
                (value - multiplier * pivot) % modulus
                for value, pivot in zip(rows[row], rows[pivot_row], strict=True)
            ]
            combinations_matrix[row] = [
                (value - multiplier * pivot) % modulus
                for value, pivot in zip(
                    combinations_matrix[row],
                    combinations_matrix[pivot_row],
                    strict=True,
                )
            ]
        pivot_columns.append(column)
        pivot_row += 1
    contradiction_row = next(
        (
            row
            for row in range(pivot_row, equation_count)
            if not any(rows[row][:-1]) and rows[row][-1] % modulus
        ),
        None,
    )
    if contradiction_row is None:
        raise ValueError("unit-pivot elimination found no modular contradiction")
    multipliers = tuple(combinations_matrix[contradiction_row])
    return {
        "modulus": modulus,
        "unit_pivot_rank": pivot_row,
        "pivot_columns": pivot_columns,
        "sparse_multipliers": [
            {"constraint_index": index, "multiplier": multiplier}
            for index, multiplier in enumerate(multipliers)
            if multiplier
        ],
        "accumulated_lhs": rows[contradiction_row][:-1],
        "accumulated_rhs": rows[contradiction_row][-1],
    }


def verify_modular_contradiction(
    constraints: tuple[BiplanePhaseConstraint, ...],
    variable_count: int,
    witness: dict[str, Any],
) -> bool:
    modulus = int(witness.get("modulus", 0))
    if modulus != _PHASE:
        return False
    multipliers = [0] * len(constraints)
    records = witness.get("sparse_multipliers")
    if not isinstance(records, list):
        return False
    seen = set()
    for record in records:
        if not isinstance(record, dict):
            return False
        index = record.get("constraint_index")
        multiplier = record.get("multiplier")
        if (
            not isinstance(index, int)
            or isinstance(index, bool)
            or not 0 <= index < len(constraints)
            or index in seen
            or not isinstance(multiplier, int)
            or isinstance(multiplier, bool)
            or not 0 < multiplier < modulus
        ):
            return False
        seen.add(index)
        multipliers[index] = multiplier
    accumulated_lhs = [0] * variable_count
    accumulated_rhs = 0
    for multiplier, constraint in zip(multipliers, constraints, strict=True):
        if not multiplier:
            continue
        for variable, coefficient in constraint.terms:
            accumulated_lhs[variable] = (
                accumulated_lhs[variable] + multiplier * coefficient
            ) % modulus
        accumulated_rhs = (
            accumulated_rhs + multiplier * constraint.rhs
        ) % modulus
    return (
        not any(accumulated_lhs)
        and accumulated_rhs != 0
        and witness.get("accumulated_lhs") == accumulated_lhs
        and witness.get("accumulated_rhs") == accumulated_rhs
    )


def classify_cgw_11_5_4() -> dict[str, Any]:
    """Return the complete published-support-plus-exact-lift certificate."""

    support = cyclic_biplane_support()
    support_audit = audit_biplane_support(support)
    matrix = support_matrix(support)
    edges, constraints = biplane_phase_system(support)
    constraint_payloads = [constraint.to_json(edges) for constraint in constraints]
    witness = find_modular_contradiction(constraints, len(edges))
    if not verify_modular_contradiction(constraints, len(edges), witness):
        raise ArithmeticError("internally generated modular witness is invalid")
    witness["explanation"] = (
        "the selected row-pair congruences have zero total left side but "
        "right side 2 modulo 4"
    )
    witness["selected_equations"] = [
        {
            **record,
            "row_pair": constraint_payloads[record["constraint_index"]]["row_pair"],
            "overlap_columns": constraint_payloads[record["constraint_index"]][
                "overlap_columns"
            ],
        }
        for record in witness["sparse_multipliers"]
    ]
    return {
        "schema": CGW_11_5_4_SCHEMA,
        "object": "complex generalized weighing matrix",
        "parameters": {"order": _ORDER, "weight": _WEIGHT, "phase": _PHASE},
        "genuinely_non_hadamard": True,
        "zero_entries_per_row": _ORDER - _WEIGHT,
        "literature_status": "published nonexistence; omitted lift calculation reproduced",
        "novelty_claimed": False,
        "equivalence_contract": CGW_EQUIVALENCE_CONTRACT,
        "external_support_reduction": CGW_11_5_4_SOURCE,
        "canonical_support": {
            "construction": "cyclic development of the (11,5,2) difference set",
            "difference_set": list(_BIPLANE_DIFFERENCE_SET),
            "support_rows": [list(row) for row in support],
            "support_matrix": [list(row) for row in matrix],
            "audit": support_audit,
        },
        "phase_model": {
            "coefficient_group": "mu_4 = {1,i,-1,-i}",
            "edge_variable_count": len(edges),
            "edge_variables": [
                {"index": index, "row": row, "column": column}
                for index, (row, column) in enumerate(edges)
            ],
            "edge_variable_digest": _canonical_json_digest(edges),
            "gauge_fixed_variable_count": support_audit[
                "gauge_fixed_phase_variable_count"
            ],
            "raw_phase_assignment_count": str(_PHASE ** len(edges)),
            "gauge_fixed_phase_assignment_count": str(
                _PHASE ** support_audit["gauge_fixed_phase_variable_count"]
            ),
            "two_term_vanishing_rule": (
                "zeta_4^a + zeta_4^b = 0 iff a-b = 2 mod 4"
            ),
            "constraint_count": len(constraints),
            "constraints": constraint_payloads,
            "constraint_system_digest": _canonical_json_digest(constraint_payloads),
        },
        "inconsistency_certificate": witness,
        "completeness_accounting": {
            "support_classes_from_published_reduction": 1,
            "support_classes_machine_enumerated": False,
            "support_classes_checked_for_lifts": 1,
            "support_lift_systems": 1,
            "orthogonality_constraints_checked": len(constraints),
            "inconsistent_lift_systems": 1,
            "consistent_lift_systems": 0,
            "labeled_lift_count": 0,
            "classical_monomial_class_count": 0,
            "extended_class_count": 0,
        },
        "proof_obligations": CGW_11_5_4_PROOF_OBLIGATIONS,
        "machine_checked_scope": (
            "complete mu_4 lifting problem for the unique 2-(11,5,2) support"
        ),
        "full_classification_basis": (
            "Egan Proposition 2.15 support reduction plus the exact no-lift "
            "certificate recorded here"
        ),
        "lift_classification_complete": True,
        "classification_complete_conditional_on_published_support_reduction": True,
        "support_classification_machine_complete": False,
        "end_to_end_machine_classification_complete": False,
        "classification_complete": True,
        "nonexistence_certified": True,
    }


def verify_cgw_11_5_4_certificate(payload: dict[str, Any]) -> dict[str, Any]:
    """Verify the stored support, full congruence system, and contradiction."""

    if payload.get("schema") != CGW_11_5_4_SCHEMA:
        raise ValueError("invalid CGW(11,5;4) certificate schema")
    if payload.get("parameters") != {
        "order": _ORDER,
        "weight": _WEIGHT,
        "phase": _PHASE,
    }:
        raise ValueError("wrong CGW parameters")
    if payload.get("equivalence_contract") != CGW_EQUIVALENCE_CONTRACT:
        raise ValueError("wrong CGW equivalence contract")
    if payload.get("external_support_reduction") != CGW_11_5_4_SOURCE:
        raise ValueError("wrong published support-reduction dependency")
    support = cyclic_biplane_support()
    matrix = support_matrix(support)
    support_payload = payload.get("canonical_support", {})
    if (
        support_payload.get("difference_set") != list(_BIPLANE_DIFFERENCE_SET)
        or support_payload.get("support_rows") != [list(row) for row in support]
        or support_payload.get("support_matrix") != [list(row) for row in matrix]
        or support_payload.get("audit") != audit_biplane_support(support)
    ):
        raise ValueError("canonical biplane support or its audit is inconsistent")
    edges, constraints = biplane_phase_system(support)
    constraint_payloads = [constraint.to_json(edges) for constraint in constraints]
    phase_model = payload.get("phase_model", {})
    if (
        phase_model.get("coefficient_group") != "mu_4 = {1,i,-1,-i}"
        or phase_model.get("edge_variable_count") != len(edges)
        or phase_model.get("edge_variables")
        != [
            {"index": index, "row": row, "column": column}
            for index, (row, column) in enumerate(edges)
        ]
        or phase_model.get("edge_variable_digest") != _canonical_json_digest(edges)
        or phase_model.get("gauge_fixed_variable_count") != 34
        or phase_model.get("raw_phase_assignment_count") != str(4**55)
        or phase_model.get("gauge_fixed_phase_assignment_count") != str(4**34)
        or phase_model.get("two_term_vanishing_rule")
        != "zeta_4^a + zeta_4^b = 0 iff a-b = 2 mod 4"
        or phase_model.get("constraint_count") != len(constraints)
        or phase_model.get("constraints") != constraint_payloads
        or phase_model.get("constraint_system_digest")
        != _canonical_json_digest(constraint_payloads)
    ):
        raise ValueError("CGW phase model is inconsistent")
    witness = payload.get("inconsistency_certificate")
    if not isinstance(witness, dict) or not verify_modular_contradiction(
        constraints,
        len(edges),
        witness,
    ):
        raise ValueError("invalid CGW modular contradiction witness")
    derived_witness = find_modular_contradiction(constraints, len(edges))
    for field in (
        "modulus",
        "unit_pivot_rank",
        "pivot_columns",
        "sparse_multipliers",
        "accumulated_lhs",
        "accumulated_rhs",
    ):
        if witness.get(field) != derived_witness[field]:
            raise ValueError("stored CGW contradiction is not the deterministic witness")
    expected_selected_equations = [
        {
            **record,
            "row_pair": constraint_payloads[record["constraint_index"]]["row_pair"],
            "overlap_columns": constraint_payloads[record["constraint_index"]][
                "overlap_columns"
            ],
        }
        for record in witness["sparse_multipliers"]
    ]
    if (
        witness.get("selected_equations") != expected_selected_equations
        or witness.get("explanation")
        != (
            "the selected row-pair congruences have zero total left side but "
            "right side 2 modulo 4"
        )
    ):
        raise ValueError("stored CGW contradiction explanation is inconsistent")
    accounting = payload.get("completeness_accounting", {})
    if accounting != {
        "support_classes_from_published_reduction": 1,
        "support_classes_machine_enumerated": False,
        "support_classes_checked_for_lifts": 1,
        "support_lift_systems": 1,
        "orthogonality_constraints_checked": 55,
        "inconsistent_lift_systems": 1,
        "consistent_lift_systems": 0,
        "labeled_lift_count": 0,
        "classical_monomial_class_count": 0,
        "extended_class_count": 0,
    }:
        raise ValueError("CGW completeness accounting is inconsistent")
    if payload.get("proof_obligations") != CGW_11_5_4_PROOF_OBLIGATIONS:
        raise ValueError("CGW proof obligations are unresolved")
    if any(
        payload.get(field) is not True
        for field in (
            "genuinely_non_hadamard",
            "lift_classification_complete",
            "classification_complete_conditional_on_published_support_reduction",
            "classification_complete",
            "nonexistence_certified",
        )
    ):
        raise ValueError("CGW classification completion flags are unresolved")
    if any(
        payload.get(field) is not False
        for field in (
            "support_classification_machine_complete",
            "end_to_end_machine_classification_complete",
        )
    ):
        raise ValueError("CGW machine-completeness scope is overstated")
    if payload.get("novelty_claimed") is not False:
        raise ValueError("the reproduction certificate must not claim novelty")
    return {
        "schema": "hadamard.cgw-certificate-verification.v1",
        "parameters": payload["parameters"],
        "support_digest": support_payload["audit"]["support_digest"],
        "constraint_system_digest": phase_model["constraint_system_digest"],
        "witness_term_count": len(witness["sparse_multipliers"]),
        "accumulated_rhs_mod_4": witness["accumulated_rhs"],
        "classical_monomial_class_count": 0,
        "lift_classification_complete": True,
        "classification_complete_conditional_on_published_support_reduction": True,
        "support_classification_machine_complete": False,
        "end_to_end_machine_classification_complete": False,
        "certificate_valid": True,
    }


__all__ = [
    "BiplanePhaseConstraint",
    "CGW_11_5_4_SCHEMA",
    "CGW_11_5_4_SOURCE",
    "CGW_11_5_4_PROOF_OBLIGATIONS",
    "CGW_EQUIVALENCE_CONTRACT",
    "ExponentMatrix",
    "apply_cgw_monomial_transform",
    "audit_biplane_support",
    "biplane_phase_system",
    "classify_cgw_11_5_4",
    "cyclic_biplane_support",
    "find_modular_contradiction",
    "is_cgw_exponent_matrix",
    "support_matrix",
    "verify_cgw_11_5_4_certificate",
    "verify_modular_contradiction",
]
