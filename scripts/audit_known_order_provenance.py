#!/usr/bin/env python3
"""Audit the positive prior-construction provenance of the 31 known orders.

Run without arguments to regenerate the certificate, or use ``--check`` to
require that its complete contents agree with the current public inputs.
Only the Python standard library is needed. The audit checks ledger links,
source hashes and integer construction hypotheses. It does not reconstruct
prior matrices, download the cited papers, or test matrix equivalence.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from math import isqrt
from pathlib import Path
import re
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = "certificates/scalar-lift/known-order-provenance.json"
INVENTORY = "certificates/scalar-lift/completed-order-inventory.json"
EXPECTED_KNOWN_COUNT = 31
REFERENCES = {
    "paley": {
        "authors": "R. E. A. C. Paley",
        "title": "On Orthogonal Matrices",
        "journal": "Journal of Mathematics and Physics",
        "volume": "12", "year": 1933, "pages": "311–320",
        "doi": "10.1002/sapm1933121311",
        "url": "https://doi.org/10.1002/sapm1933121311",
    },
    "cooperwallis": {
        "authors": "Joan Cooper and Jennifer Wallis",
        "title": "A construction for Hadamard arrays",
        "journal": "Bulletin of the Australian Mathematical Society",
        "volume": "7", "year": 1972, "pages": "269–277",
        "doi": "10.1017/S0004972700045019",
        "url": "https://doi.org/10.1017/S0004972700045019",
        "location": "Theorem 3 and Corollary 4, pp. 274–275; odd T inputs through 19",
    },
    "turyn1972": {
        "authors": "R. J. Turyn",
        "title": "An infinite class of Williamson matrices",
        "journal": "Journal of Combinatorial Theory, Series A",
        "volume": "12", "year": 1972, "pages": "319–321",
        "doi": "10.1016/0097-3165(72)90095-7",
        "url": "https://doi.org/10.1016/0097-3165(72)90095-7",
    },
    "whiteman1973": {
        "authors": "A. L. Whiteman",
        "title": "An infinite family of Hadamard matrices of Williamson type",
        "journal": "Journal of Combinatorial Theory, Series A",
        "volume": "14", "year": 1973, "pages": "334–340",
        "doi": "10.1016/0097-3165(73)90010-1",
        "url": "https://doi.org/10.1016/0097-3165(73)90010-1",
        "location": "New proof of Turyn's WT((q+1)/2) theorem for prime powers q ≡ 1 mod 4",
    },
    "seberry-yamada-1990": {
        "authors": "J. Seberry and M. Yamada",
        "title": ("On the products of Hadamard matrices, Williamson matrices and "
                  "other orthogonal matrices using M-structures"),
        "journal": "Journal of Combinatorial Mathematics and Combinatorial Computing",
        "volume": "7", "year": 1990, "pages": "97–137",
        "url": "https://documents.uow.edu.au/~jennie/WEBPDF/139_1990.pdf",
        "location": "Theorem 4, printed p. 105",
    },
    "survey": {
        "authors": "R. Egan, I. S. Kotsireas, S. London and P. Ó Catháin",
        "title": "Recent progress in the construction of Hadamard matrices",
        "arxiv": "2610.27175v2",
        "url": "https://arxiv.org/html/2610.27175v2#S4.T4",
        "location": ("Theorem 2, Tables 3–4 and the paragraph following Table 4: "
                     "Paley's sixteen-block arrangement and its product with T-matrices"),
    },
    "djk2014": {
        "authors": "D. Ž. Đoković, O. Golubitsky and I. S. Kotsireas",
        "title": "Some new orders of Hadamard and skew-Hadamard matrices",
        "journal": "Journal of Combinatorial Designs",
        "volume": "22", "year": 2014, "pages": "270–277",
        "doi": "10.1002/jcd.21358", "arxiv": "1301.3671v2",
        "url": "https://arxiv.org/html/1301.3671v2#S3.SS3",
        "location": "Section 3.3: explicit cyclic supplementary difference sets on Z/631Z",
    },
    "djokovic2009": {
        "authors": "D. Ž. Đoković",
        "title": "Hadamard matrices of small order and Yang conjecture",
        "arxiv": "0912.5091",
        "url": "https://arxiv.org/html/0912.5091#S3",
        "location": ("Section 2: BS(r,s) → T(r+s); Section 3: explicit "
                     "2159 = 17·127, with earlier WT127 from good matrices (1993)"),
    },
    "best2013": {
        "authors": "D. Best, D. Ž. Đoković, H. Kharaghani and H. Ramp",
        "title": "Turyn-type sequences: Classification, Enumeration and Construction",
        "journal": "Journal of Combinatorial Designs",
        "volume": "21", "year": 2013, "pages": "24–35",
        "doi": "10.1002/jcd.21318", "arxiv": "1206.4107",
        "url": "https://arxiv.org/html/1206.4107",
        "location": "TT(38) exists; the introduction and Section 5 give TT(n) → BS(2n−1,n)",
    },
}


class AuditError(ValueError):
    """A public input or a derived construction hypothesis failed its audit."""


def require(condition: bool, message: str) -> None:
    """Check an invariant even when Python runs with optimization enabled."""
    if not condition:
        raise AuditError(message)


def is_prime(n: int) -> bool:
    """Use complete trial division; the largest present input is below 24 million."""
    return type(n) is int and n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def select(value: Any, selector: str) -> Any:
    """Resolve the ledger's dotted-key and integer-list-index selectors strictly."""
    position = 0
    for match in re.finditer(r"(?:^|\.)([^.\[\]]+)|\[(\d+)\]", selector):
        require(match.start() == position, f"Invalid selector syntax: {selector!r}")
        key, index = match.groups()
        try:
            value = value[int(index)] if index is not None else value[key]
        except (KeyError, IndexError, TypeError) as exc:
            raise AuditError(f"Unresolved selector: {selector!r}") from exc
        position = match.end()
    require(position == len(selector) and position > 0,
            f"Invalid selector syntax: {selector!r}")
    return value


class Sources:
    """Read each referenced public JSON input once and retain its exact hash."""

    def __init__(self, root: Path):
        self.root = root.resolve()
        self.values: dict[str, Any] = {}
        self.hashes: dict[str, str] = {}

    def load(self, relative: str) -> Any:
        if relative not in self.values:
            path = (self.root / relative).resolve()
            require(path.is_relative_to(self.root), f"Input escapes repository: {relative}")
            raw = path.read_bytes()
            self.values[relative] = json.loads(raw)
            self.hashes[relative] = hashlib.sha256(raw).hexdigest()
        return self.values[relative]


def collect_priors(sources: Sources) -> dict[int, list[dict[str, Any]]]:
    """Collect the current detailed statuses and the three historical controls."""
    selected: dict[int, list[dict[str, Any]]] = {}

    def add(order: int, prior: dict[str, Any], path: str, selector: str) -> None:
        selected.setdefault(order, []).append({
            "raw_prior": prior, "path": path, "selector": selector,
        })

    for path in (
        "experiments/scalar_lift/williamson_order_status.json",
        "experiments/scalar_lift/williamson_family_order_status.json",
        "experiments/scalar_lift/near_williamson_order_status.json",
        "experiments/scalar_lift/additional_spence_order_status.json",
    ):
        data = sources.load(path)
        entries = data.get("entries", data)
        for key, row in entries.items():
            if isinstance(row, dict) and row.get("status") == "known":
                selector = f"entries.{key}" if "entries" in data else key
                add(row["order"], row["prior_construction"], path, selector)

    path = "certificates/scalar-lift/higher-order-records.json"
    old = sources.load(path)
    selector = "periodic_T127_record.ordinary_H1524_status"
    require(select(old, selector) == "known_PaleyI_q1523", "H1524 prior status changed")
    add(1524, {"family": "Paley I", "parameters": {"q": 1523}}, path, selector)
    for index, row in enumerate(old["repeat_endpoint_statuses"]):
        if row.get("output_order") not in (3252, 6564):
            continue
        proofs = [(j, p) for j, p in enumerate(row["positive_knownness_proofs"])
                  if p["construction"] == "PaleyI"]
        require(len(proofs) == 1, f"Expected one Paley I proof for H{row['output_order']}")
        j, proof = proofs[0]
        add(row["output_order"], {"family": "Paley I", "parameters": {"q": proof["q"]}},
            path, f"repeat_endpoint_statuses[{index}].positive_knownness_proofs[{j}]")

    path = "certificates/scalar-lift/tripled-unreal-scalar-endpoints.json"
    rows = sources.load(path)["complete_endpoint_recipes"]
    matches = [(i, row) for i, row in enumerate(rows) if row["endpoint_order"] == 885924]
    require(len(matches) == 1, "Expected one H885924 endpoint recipe")
    index, row = matches[0]
    proof = row["known_order_construction"]
    add(885924, {"family": proof["method"], "parameters": {"q": proof["prime_power"]}},
        path, f"complete_endpoint_recipes[{index}].known_order_construction")
    return selected


def normalize_prior(order: int, raw: dict[str, Any]) -> dict[str, Any]:
    """Check the positive integer hypotheses for one cited earlier construction."""
    family = raw["family"]
    parameters = raw.get("parameters", raw)
    q = parameters.get("q", parameters.get("prime_power_parameter"))
    checks: dict[str, bool] = {}
    note = ""
    if family in ("Paley I", "PaleyI"):
        kind = "Paley I"
        checks = {"q_is_prime": is_prime(q), "q_mod4_is3": q % 4 == 3,
                  "order_identity": order == q + 1}
        result_parameters = {"q": q}
        refs = ["paley"]
    elif family in ("Paley II", "PaleyII"):
        kind = "Paley II"
        checks = {"q_is_prime": is_prime(q), "q_mod4_is1": q % 4 == 1,
                  "order_identity": order == 2 * (q + 1)}
        result_parameters = {"q": q}
        refs = ["paley"]
    elif "sixteen-block" in family:
        kind = "Paley-block product"
        t = parameters.get("t", parameters.get("T_order"))
        s = parameters.get("s", (q + 1) // 4)
        checks = {"q_is_prime": is_prime(q), "q_mod8_is3": q % 8 == 3,
                  "base_order_identity": q + 1 == 4 * s,
                  "order_identity": order == 4 * s * t,
                  "T_input_classical_odd_through19": t % 2 == 1 and 1 <= t <= 19}
        result_parameters = {"q": q, "s": s, "T_order": t}
        refs = ["seberry-yamada-1990", "survey", "cooperwallis"]
    elif "Goethals" in family:
        kind = "Goethals–Seidel / cyclic SDS"
        v = parameters["cyclic_group_order"]
        sizes = [315, 330, 330, 330]
        lam = sum(sizes) - v
        checks = {"published_cyclic_group_order": v == 631,
                  "order_identity": order == 4 * v,
                  "SDS_parameter_identity": sum(k * (k - 1) for k in sizes) == lam * (v - 1)}
        result_parameters = {"v": v, "block_sizes": sizes, "lambda": lam}
        refs = ["djk2014"]
        note = ("Section 3.3 gives λ=674; the introduction has a typographical 647. "
                "Existence of the literal SDSs is taken from the cited paper, not inferred "
                "from the parameter identity. Its SDS inequivalence results do not "
                "compare the resulting Hadamard matrices with our representatives.")
    elif "Cooper" in family:
        t = parameters.get("T_order", parameters.get("known_T_order"))
        w = parameters.get("Williamson_order", parameters.get("known_Williamson_type_order"))
        result_parameters = {"T_order": t, "Williamson_order": w}
        checks = {"order_identity": order == 4 * t * w}
        if q is not None:
            kind = "Cooper–Wallis / Turyn–Whiteman"
            result_parameters["q"] = q
            checks.update(q_is_prime=is_prime(q), q_mod4_is1=q % 4 == 1,
                          Williamson_order_identity=2 * w == q + 1)
            refs = ["cooperwallis", "turyn1972", "whiteman1973"]
        else:
            kind = "Cooper–Wallis / prior WT127"
            checks["published_prior_WT127_with_T17"] = (w, t) == (127, 17)
            refs = ["cooperwallis", "djokovic2009"]
        if t == 113:
            refs += ["best2013", "djokovic2009"]
            checks["TT38_to_BS75_38_to_T113"] = 2 * 38 - 1 + 38 == t
            result_parameters["T_source"] = "TT(38) → BS(75,38) → T113"
        else:
            checks["T_input_classical_odd_through19"] = t % 2 == 1 and 1 <= t <= 19
    else:
        raise AuditError(f"Unrecognized prior-construction family for H{order}: {family!r}")
    failed = [name for name, passed in checks.items() if passed is not True]
    require(not failed, f"H{order} failed integer hypotheses: {', '.join(failed)}")
    return {"prior_construction": kind, "parameters": result_parameters,
            "reference_keys": refs, "integer_hypotheses_checked": checks, "note": note}


def build_report(root: Path = ROOT) -> dict[str, Any]:
    """Rebuild the complete certificate from the current repository inputs."""
    sources = Sources(root)
    ledger = sources.load(INVENTORY)
    known_rows = [row for row in ledger["orders"] if row["status"] == "known"]
    orders = [row["order"] for row in known_rows]
    require(len(orders) == EXPECTED_KNOWN_COUNT, "Known-order count differs from the audited 31")
    require(orders == sorted(set(orders)), "Known-order ledger is duplicated or unsorted")
    require(orders == ledger["completed_known_real_orders"], "Known-order list differs from ledger rows")
    require(ledger["summary"]["known_recovered_orders"] == len(orders), "Known-order summary differs")
    priors = collect_priors(sources)
    require(set(priors) == set(orders), "Detailed known statuses differ from the completed known-order set")
    records = []
    endpoint_selectors = status_selectors = prior_proofs = 0
    endpoint_files: set[str] = set()
    for row in known_rows:
        order = row["order"]
        require(row["construction_complete"] is True, f"H{order} is not marked complete")
        require(bool(row["positive_certificates"]), f"H{order} lacks a positive endpoint certificate")
        require(bool(row["status_sources"]), f"H{order} lacks a status source")
        for certificate in row["positive_certificates"]:
            endpoint = select(sources.load(certificate["path"]), certificate["selector"])
            possible_orders = [endpoint[key] for key in (
                "order", "output_order", "endpoint_order", "known_order", "hadamard_order", "real_order")
                if type(endpoint.get(key)) is int]
            require(bool(possible_orders) and all(value == order for value in possible_orders),
                    f"H{order} endpoint selector has a different order")
            endpoint_selectors += 1
            endpoint_files.add(certificate["path"])
        for status in row["status_sources"]:
            value = select(sources.load(status["path"]), status["selector"])
            require(value.get("order") == order and value.get("status") == "known",
                    f"H{order} ledger status selector disagrees")
            status_selectors += 1
        proofs = [normalize_prior(order, item["raw_prior"]) for item in priors[order]]
        prior_proofs += len(proofs)
        record = {"order": order, **proofs[0]}
        record.update({
            "prior_status_provenance": [
                {"path": item["path"], "selector": item["selector"], "normalized_proof": proof}
                for item, proof in zip(priors[order], proofs)
            ],
            "our_positive_endpoint_certificates": row["positive_certificates"],
            "ledger_status_sources": row["status_sources"],
            "matrix_equivalence_tested_in_this_audit": False,
        })
        records.append(record)
    return {
        "schema": "known-completed-order-positive-provenance-v1",
        "producer": "scripts/audit_known_order_provenance.py",
        "source_inventory": INVENTORY,
        "inventory_checked_date": ledger["checked_date"],
        "scope": ("Positive prior existence-order provenance only. This audit checks ledger "
                  "bindings, input hashes and integer construction hypotheses. Prior matrix "
                  "existence is attributed to the cited constructions; those matrices are not "
                  "reconstructed here. No matrix-equivalence test is performed by this audit."),
        "known_order_count": len(records),
        "orders": orders,
        "construction_counts": dict(Counter(row["prior_construction"] for row in records)),
        "all_integer_hypotheses_passed": True,
        "endpoint_selectors_checked": endpoint_selectors,
        "endpoint_certificate_files_checked": len(endpoint_files),
        "ledger_status_selectors_checked": status_selectors,
        "detailed_prior_proofs_checked": prior_proofs,
        "rows": records,
        "references": REFERENCES,
        "source_sha256": dict(sorted(sources.hashes.items())),
    }


def canonical_bytes(report: dict[str, Any]) -> bytes:
    """Serialize the report reproducibly, without local paths or timestamps."""
    return (json.dumps(report, indent=2, ensure_ascii=False) + "\n").encode("utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=ROOT, help="repository root (default: this checkout)")
    parser.add_argument("--output", type=Path, help=f"certificate path (default: REPO/{CERTIFICATE})")
    parser.add_argument("--check", action="store_true", help="require the existing certificate to match exactly")
    args = parser.parse_args(argv)
    output = args.output if args.output is not None else args.repo / CERTIFICATE
    try:
        report = build_report(args.repo)
        expected = canonical_bytes(report)
        if args.check:
            require(output.read_bytes() == expected,
                    f"Certificate differs from current inputs: {output}; regenerate without --check")
        else:
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_bytes(expected)
    except (AuditError, OSError, ValueError, KeyError, TypeError) as exc:
        print(f"Known-order provenance audit failed: {exc}", file=sys.stderr)
        return 1
    print(json.dumps({"ok": True, "mode": "check" if args.check else "write",
                      "known_orders": report["known_order_count"],
                      "endpoint_selectors_checked": report["endpoint_selectors_checked"],
                      "ledger_status_selectors_checked": report["ledger_status_selectors_checked"],
                      "detailed_prior_proofs_checked": report["detailed_prior_proofs_checked"],
                      "construction_counts": report["construction_counts"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
