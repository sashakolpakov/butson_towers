"""Complete minimal positive vanishing sums in mu30 and mu120.

Square-free relative order and the complete 2pq classification imply that
these catalogues have no missing larger-weight atoms. Original matrix
columns remain distinct in the lazy indexed matcher, even at equal roots.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from functools import cache
from itertools import combinations, product
from pathlib import Path

SUPPORTED_MODULI = (30, 120)
PHI30 = (1, 1, 0, -1, -1, -1, 0, 1, 1)
SQUARE_FREE_REFERENCE = "https://arxiv.org/pdf/math/9511209"
CLASSIFICATION_REFERENCE = "https://arxiv.org/html/2008.11268v2#S3"


def _modulus(value):
    if type(value) is not int or value not in SUPPORTED_MODULI:
        raise ValueError("this complete catalogue supports exactly phases 30 and 120")
    return value


def cyclotomic_polynomial(modulus):
    modulus = _modulus(modulus)
    scale = modulus // 30
    return tuple(PHI30[index // scale] if index % scale == 0 else 0
                 for index in range(8 * scale + 1))


@cache
def root_vectors(modulus):
    """Integral power-basis vectors, computed without numerical arithmetic."""
    polynomial = cyclotomic_polynomial(modulus)
    degree = len(polynomial) - 1
    current = (1,) + (0,) * (degree - 1)
    table = []
    for _ in range(modulus):
        table.append(current)
        shifted = [0, *current[:-1]]
        current = tuple(value - current[-1] * polynomial[j]
                        for j, value in enumerate(shifted))
    if current != table[0]:
        raise ValueError("the exact root table did not return to one")
    return tuple(table)


def exact_sum(exponents, modulus):
    roots = root_vectors(modulus)
    total = [0] * len(roots[0])
    for exponent in exponents:
        root = roots[exponent % modulus]
        total = [a + b for a, b in zip(total, root, strict=True)]
    return tuple(total)


def canonical_rotation(profile, modulus=30):
    profile = tuple(profile)
    if not profile:
        raise ValueError("a nonempty root multiset is required")
    return min(tuple(sorted((value - shift) % modulus for value in profile))
               for shift in profile)


@cache
def minimal_profiles(modulus=120):
    """All literal root sets, including rotations; no weight truncation."""
    modulus = _modulus(modulus)
    if modulus == 120:
        # Four cosets of mu30 in mu120; rotations within each coset are
        # already present in the complete phase30 catalogue.
        return tuple(sorted((tuple(4 * exponent + offset for exponent in profile)
                             for offset in range(4) for profile in minimal_profiles(30)),
                            key=lambda item: (len(item), item)))
    profiles = set()
    for prime in (2, 3, 5):
        for rotation in range(30):
            profiles.add(tuple(sorted((rotation + 30 * j // prime) % 30
                                      for j in range(prime))))
    # theta=zeta30^6, omega=zeta30^10. Replace theta^j by
    # -omega*theta^j and -omega^2*theta^j at a proper nonempty subset J.
    for size in range(1, 5):
        for selected in combinations(range(5), size):
            chosen = set(selected)
            exponents = [6 * j for j in range(5) if j not in chosen]
            exponents += [(6 * j + offset) % 30 for j in chosen for offset in (5, 25)]
            for rotation in range(30):
                profiles.add(tuple(sorted((value + rotation) % 30 for value in exponents)))
    return tuple(sorted(profiles, key=lambda item: (len(item), item)))


@cache
def rotation_representatives(modulus=30):
    modulus = _modulus(modulus)
    return tuple(sorted({canonical_rotation(profile, modulus)
                         for profile in minimal_profiles(modulus)},
                        key=lambda item: (len(item), item)))


def iter_indexed_blocks(residues, modulus=120):
    """Lazily match all minimal atoms to original, distinct column positions.

    No permutation of the matrix columns is quotient out. Equal root values
    in different positions give different choices. Output can be large;
    callers should consume this iterator without materializing it blindly.
    """
    modulus = _modulus(modulus)
    values = tuple(residues)
    if any(type(value) is not int for value in values):
        raise ValueError("root exponents must be exact integers")
    buckets = [[] for _ in range(modulus)]
    for index, value in enumerate(values):
        buckets[value % modulus].append(index)
    for profile in minimal_profiles(modulus):
        choices = [buckets[value] for value in profile]
        if all(choices):
            for columns in product(*choices):
                yield tuple(sorted(columns))


def proper_subsum_audit(profile, modulus=30):
    """Independent exhaustive minimality check for one short representative."""
    roots = root_vectors(modulus)
    sums = [(0,) * len(roots[0])]
    for exponent in profile:
        root = roots[exponent % modulus]
        sums += [tuple(a + b for a, b in zip(value, root, strict=True))
                 for value in sums]
    if any(sums[-1]) or any(not any(value) for value in sums[1:-1]):
        raise ValueError("the proposed atom has a nonzero total or a proper zero subsum")
    return len(sums) - 2


def certificate():
    representatives = rotation_representatives(30)
    checked = sum(proper_subsum_audit(profile) for profile in representatives)
    reports = []
    for modulus in SUPPORTED_MODULI:
        profiles = minimal_profiles(modulus)
        if any(any(exact_sum(profile, modulus)) for profile in profiles):
            raise ValueError("a catalogue atom failed exact cyclotomic arithmetic")
        reports.append({
            "modulus": modulus,
            "literal_minimal_profiles": [list(profile) for profile in profiles],
            "literal_profile_count": len(profiles),
            "weight_histogram": {str(k): v for k, v in sorted(Counter(map(len, profiles)).items())},
            "rotation_class_count": len(rotation_representatives(modulus)),
            "maximum_minimal_weight": max(map(len, profiles)),
            "all_exact_sums_verified": True,
            "larger_minimal_weights_excluded_by_theorem": True,
        })
    if ([row["literal_profile_count"] for row in reports] != [211, 844]
            or len(representatives) != 9 or checked != 1370):
        raise ValueError("the complete phase30/120 catalogue counts changed")
    return {
        "schema": "vanishing-known-catalogue-v1",
        "proof": "docs/VANISHING_KNOWN_CATALOGUE.md",
        "references": {
            "square_free_relative_order": {"url": SQUARE_FREE_REFERENCE, "result": "Lam-Leung Corollary 3.2"},
            "complete_2pq_classification": {"url": CLASSIFICATION_REFERENCE, "result": "Christie-Dykema-Klep Theorem 3.2"},
        },
        "phase30_rotation_representatives": [list(profile) for profile in representatives],
        "proper_subsums_independently_checked": checked,
        "catalogues": reports,
        "multiplicities_in_every_minimal_profile": 1,
        "original_column_labels_preserved_by_indexed_matcher": True,
        "general_modulus_or_weight12_catalogue_modified": False,
        "Butson_source_constructed": False,
        "orthogonality_of_multiple_rows_not_implied_by_one_vanishing_sum": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument("--output", type=Path)
    actions.add_argument("--check", type=Path)
    args = parser.parse_args()
    report = certificate()
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    elif json.loads(args.check.read_text()) != report:
        raise ValueError("the complete vanishing catalogue failed exact replay")
    print(json.dumps({"phase30_profiles": 211, "phase120_profiles": 844,
                      "rotation_classes": 9, "maximum_minimal_weight": 9,
                      "Butson_source_constructed": False}))


if __name__ == "__main__":
    main()
