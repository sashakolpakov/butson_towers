"""Exact odd-integer-relation obstructions to uniform scalar-to-sign maps.

The map may use any field tower, uniform domain/codomain scalar bases, and
representation multiplicity. Nonuniform coordinate maps are not excluded.
No source construction, gauge search, or numerical matrix operation occurs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np

from experiments.scalar_lift.fifth_source_assembly import root_vectors
from experiments.scalar_lift.vanishing_known_catalogue import rotation_representatives

SCOPE = ("every uniform Q-linear map of the actual scalar alphabet to equally sized "
         "sign blocks, including all scalar-restriction towers, uniform input/output "
         "basis choices, and representation multiplicities")


def _small_prime_support(value):
    for prime in (2, 3, 5):
        while value % prime == 0:
            value //= prime
    return value == 1


def verify_relation(phase, used, terms):
    """An integer relation and odd coefficient sum, checked in characteristic zero."""
    phase = int(phase)
    used = {int(e) % phase for e in used}
    coefficients = Counter()
    for term in terms:
        exponent = int(term["exponent"]) % phase
        coefficient = int(term["coefficient"])
        if exponent not in used:
            raise ValueError("a relation root does not occur in the actual alphabet")
        coefficients[exponent] += coefficient
    coefficients = {e: c for e, c in sorted(coefficients.items()) if c}
    roots = root_vectors(phase)
    total = sum((c * roots[e] for e, c in coefficients.items()),
                np.zeros(roots.shape[1], dtype=np.int64))
    coefficient_sum = sum(coefficients.values())
    if np.any(total) or coefficient_sum % 2 != 1:
        raise ValueError("not an exact odd integer relation")
    return {"terms": [{"exponent": e, "coefficient": c} for e, c in coefficients.items()],
            "zero_integral_vector": total.tolist(), "coefficient_sum": coefficient_sum,
            "odd_coefficient_sum": True,
            "signed_weight": sum(abs(c) for c in coefficients.values())}


def alphabet_obstruction(phase, used):
    """Find an exact signed odd relation; SKIP never excludes another mechanism.

    For conductors supported on {2,3,5}, this is a complete test for this
    particular obstruction, using the proved complete minimal q30 catalogue.
    It does not mistake a mod-2 relation for an integral characteristic-zero
    relation. For other conductors the same positive witnesses remain valid,
    while a negative result makes no completeness claim.
    """
    if type(phase) is not int or phase < 2:
        raise ValueError("phase must be an integer at least two")
    used = sorted({int(e) % phase for e in used})
    complete = _small_prime_support(phase)
    common = math.lcm(phase, 30)
    base = {"phase": phase, "used_alphabet": used, "scope": SCOPE,
            "nonuniform_coordinate_maps_excluded": False,
            "source_row_column_gauges_searched": False,
            "field_towers_or_bases_enumerated": False}
    if common > 3840:
        return {**base, "status": "SKIP", "excluded": False,
                "reason": "finite catalogue embedding budget exceeded",
                "complete_for_odd_integer_relations": False}
    scale = common // phase
    signed = {}
    for exponent in used:
        signed.setdefault(exponent * scale, (exponent, 1))
    for exponent in used:
        signed.setdefault((exponent * scale + common // 2) % common, (exponent, -1))
    tested = 0
    for profile in rotation_representatives(30):
        if len(profile) % 2 == 0:
            continue
        offsets = [e * (common // 30) for e in profile]
        # Every representative contains zero, so its rotation must itself
        # belong to the signed alphabet. This exhausts all possible matches.
        for rotation in sorted(signed):
            tested += 1
            atom = [(rotation + e) % common for e in offsets]
            if not all(e in signed for e in atom):
                continue
            terms = [{"exponent": signed[e][0], "coefficient": signed[e][1]} for e in atom]
            relation = verify_relation(phase, used, terms)
            return {**base, "status": "EXCLUDED", "excluded": True,
                    "relation": relation, "embedding_modulus": common,
                    "minimal_positive_atom": atom, "profiles_checked": tested,
                    "complete_for_odd_integer_relations": complete,
                    "reason": "exact integer scalar relation has odd coefficient sum"}
    return {**base, "status": "SKIP", "excluded": False,
            "reason": "no odd integer relation" if complete else "no tested odd relation",
            "complete_for_odd_integer_relations": complete, "profiles_checked": tested,
            "embedding_modulus": common,
            "sign_flattenability_proved": False}


def source_audit(source):
    """Bind the alphabet gate to a fixed verified registry source and occurrences."""
    if source.proof.get("full_Gram_verified") is not True:
        raise ValueError("the fixed input needs its own completed Gram proof")
    gate = alphabet_obstruction(source.phase, source.used_alphabet)
    occurrences = []
    if gate["excluded"]:
        for term in gate["relation"]["terms"]:
            exponent = term["exponent"]
            if source.exponents is not None:
                position = np.argwhere(source.exponents == exponent)[0]
                row, column = map(int, position)
            elif source.id in ("character-73", "character-241"):
                # Its zero row is the literal group-ring coefficient table.
                row = 0
                column = next(j for j in range(source.order) if source.entry(0, j) == exponent)
            else:
                raise ValueError("no literal occurrence accessor registered for this implicit source")
            if source.entry(row, column) != exponent:
                raise ArithmeticError("actual source occurrence failed")
            occurrences.append({**term, "row": row, "column": column,
                                "actual_entry_exponent": exponent})
    return {"source_id": source.id, "order": source.order, "phase": source.phase,
            "source_hash": source.source_hash, "full_source_Gram_previously_verified": True,
            "source_proof_sha256": hashlib.sha256(json.dumps(
                source.proof, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
            "alphabet_gate": gate, "actual_root_occurrences": occurrences}


def certificate():
    from experiments.scalar_lift.forward_source_registry import load_source

    identifiers = ("character-73", "character-241", "conference-1000-q8", "common120-35")
    reports = [source_audit(load_source(identifier)) for identifier in identifiers]
    if [row["alphabet_gate"]["excluded"] for row in reports] != [True, True, False, False]:
        raise ArithmeticError("fixed-source obstruction controls changed")
    return {"schema": "uniform-scalar-restriction-obstructions-v1", "sources": reports,
            "proof": "docs/SCALAR_RESTRICTION_OBSTRUCTIONS.md",
            "classification_reused": "docs/VANISHING_KNOWN_CATALOGUE.md",
            "complete_scope": "odd integer relations for conductors supported on primes2,3,5",
            "new_source_or_Hadamard_constructed": False,
            "nonuniform_coordinate_maps_excluded": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument("--report", type=Path)
    actions.add_argument("--check", type=Path)
    parser.add_argument("--all-registry", action="store_true")
    args = parser.parse_args()
    result = registry_certificate() if args.all_registry else certificate()
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    elif json.loads(args.check.read_text()) != result:
        raise ValueError("the obstruction certificate failed exact replay")
    reports = result["source_audits"] if args.all_registry else result["sources"]
    print(json.dumps({"fixed_sources_checked": len(reports),
                      "uniform_scalar_branches_excluded": [r["source_id"] for r in reports
                                                           if r["alphabet_gate"]["excluded"]]}))


def registry_certificate():
    """Read-only audit of every currently registered, completed fixed input."""
    from experiments.scalar_lift.forward_source_registry import load_source, sources

    descriptors = sources()
    reports = [source_audit(load_source(descriptor)) for descriptor in descriptors]
    excluded = [row["source_id"] for row in reports if row["alphabet_gate"]["excluded"]]
    retained = [row["source_id"] for row in reports if not row["alphabet_gate"]["excluded"]]
    return {"schema": "registered-scalar-restriction-obstructions-v1",
            "source_descriptors": descriptors, "source_count": len(reports),
            "source_audits": reports, "excluded_uniform_scalar_branches": excluded,
            "retained_after_this_gate": retained,
            "counts": {"excluded": len(excluded), "retained": len(retained)},
            "retained_does_not_assert_sign_flattenability": True,
            "no_source_construction_or_modification": True,
            "proof": "docs/SCALAR_RESTRICTION_OBSTRUCTIONS.md"}


if __name__ == "__main__":
    main()
