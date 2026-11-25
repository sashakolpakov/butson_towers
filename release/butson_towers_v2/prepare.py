#!/usr/bin/env python3
"""Build the explicit Hexagon source archive and validate its declarations.

This command makes no network requests and never submits a paper. Credentials
are read only to check that their bytes are absent from the explicit payload.
When present, the user's private administrative defaults are applied locally;
article-specific content and classifications remain in the metadata template.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BUNDLE = ROOT / "artifacts/hadamard-7572"
ANCILLARY = (
    "recipe.py", "t631.json", "certificate.json", "novelty.json",
    "independent-sequences.json", "independent-full-npz-audit.json",
    "independent-scalar-frame-audit.json",
)
WEIGHTED_ANCILLARY = {
    "williamson_scalar_endpoints.py": "experiments/scalar_lift/williamson_scalar_endpoints.py",
    "williamson_spence_family.py": "experiments/scalar_lift/williamson_spence_family.py",
    "williamson_cyclotomic_trace_audit.py": "experiments/scalar_lift/williamson_cyclotomic_trace_audit.py",
    "williamson_inputs.json": "experiments/scalar_lift/williamson_inputs.json",
    "williamson_order_status.json": "experiments/scalar_lift/williamson_order_status.json",
    "williamson_family_order_status.json": "experiments/scalar_lift/williamson_family_order_status.json",
    "williamson-scalar-endpoints.json": "certificates/scalar-lift/williamson-scalar-endpoints.json",
    "williamson-scalar-known-control.json": "certificates/scalar-lift/williamson-scalar-known-control.json",
    "williamson-family-endpoints.json": "certificates/scalar-lift/williamson-family-endpoints.json",
    "williamson-spence-family.json": "certificates/scalar-lift/williamson-spence-family.json",
    "williamson-cyclotomic-trace-audit.json": "certificates/scalar-lift/williamson-cyclotomic-trace-audit.json",
}
ORIGIN = "https://hexagonmath.org"
REQUIRED = (
    "title", "abstract", "subject", "mscCodes", "authorshipMode", "authors",
    "submittingContributorIndex", "sponsorName", "humanUnderstanding",
    "correctnessConfirmed", "communityInterestConfirmed", "citationsConfirmed",
    "license", "rightsConfirmed",
)
ADMINISTRATIVE_FIELDS = (
    "authors", "submittingContributorIndex", "sponsorName",
    "submittingOrganizationId", "license", "authorshipMode",
    "humanUnderstanding", "provenance", "correctnessConfirmed",
    "communityInterestConfirmed", "citationsConfirmed", "rightsConfirmed",
)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def apply_administrative_defaults(metadata, path=None):
    """Apply only fixed public request fields from the private local file."""
    path = ROOT / "hexagon_parameters.local" if path is None else Path(path)
    result = dict(metadata)
    if not path.exists():
        return result
    if path.stat().st_mode & 0o777 != 0o600:
        raise ValueError("The private Hexagon parameter file must have mode 0600")
    record = json.loads(path.read_text())
    if record.get("schema") != "hexagon-fixed-submission-parameters-v1":
        raise ValueError("Unsupported private Hexagon parameter schema")
    defaults = record["fixed_administrative_parameters"]
    if set(defaults) != set(ADMINISTRATIVE_FIELDS):
        raise ValueError("Unexpected fixed administrative parameter fields")
    paper_provenance = (result.get("provenance") or "").strip()
    fixed_provenance = defaults["provenance"].strip()
    result.update(defaults)
    if paper_provenance.startswith(fixed_provenance):
        result["provenance"] = paper_provenance
    elif paper_provenance:
        result["provenance"] = fixed_provenance + "\n\n" + paper_provenance
    return result


def main():
    source = (HERE / "main.tex").read_bytes()
    assert source == (ROOT / "manuscript/hadamard_7572_announcement.tex").read_bytes()
    build = json.loads((HERE / "build.json").read_text())
    assert digest(source) == build["source_sha256"], "Rebuild the changed source first"
    assert digest((HERE / "announcement.pdf").read_bytes()) == build["pdf_sha256"]
    certificate = json.loads((BUNDLE / "certificate.json").read_text())
    files = {"main.tex": source}
    for name in ANCILLARY:
        data = (BUNDLE / name).read_bytes()
        expected = certificate["file_sha256"].get(name)
        if expected is not None:
            assert digest(data) == expected, name
        files["anc/" + name] = data
    files["anc/frame-reuse.json"] = (ROOT / "certificates/scalar-lift/hadamard-7572-frame-reuse.json").read_bytes()
    files["anc/reuse.py"] = (HERE / "reuse.py").read_bytes()
    t127 = ROOT / "artifacts/t-matrices-127"
    t127_certificate = json.loads((t127 / "certificate.json").read_text())
    for name in ("t127.json", "recipe.py", "certificate.json", "independent-audit.json", "novelty.json"):
        data = (t127 / name).read_bytes()
        expected = t127_certificate["file_sha256"].get(name)
        if expected is not None:
            assert digest(data) == expected, name
        files["anc/t127/" + name] = data
    for name, path in WEIGHTED_ANCILLARY.items():
        files["anc/" + name] = (ROOT / path).read_bytes()
    weighted = json.loads(files["anc/williamson-scalar-endpoints.json"])
    weighted_control = json.loads(files["anc/williamson-scalar-known-control.json"])
    families = json.loads(files["anc/williamson-family-endpoints.json"])
    for report in (weighted, weighted_control, families):
        assert report["implementation_sha256"] == digest(files["anc/williamson_scalar_endpoints.py"])
        assert report["all_outputs_are_completed_literal_sign_matrix_recipes"]
    assert families["source_recurrence_implementation_sha256"] == digest(files["anc/williamson_spence_family.py"])
    assert weighted["Williamson_input_sha256"] == digest(files["anc/williamson_inputs.json"])
    assert weighted["status_record_sha256"] == digest(files["anc/williamson_order_status.json"])
    assert families["status_record_sha256"] == digest(files["anc/williamson_family_order_status.json"])
    proof_sources = json.loads((HERE / "proof-sources.json").read_text())
    proof_verification = (HERE / "lean-verification.json").read_bytes()
    proof_report = json.loads(proof_verification)
    assert proof_report["success"] and proof_report["proof_theorems_checked"] == 85
    assert not proof_report["project_axioms_or_placeholders_used"]
    assert "ButsonTowers.unreal_cubic_family" in proof_report["axioms_by_theorem"]
    assert proof_sources["source_package_sha256"] == digest(
        (ROOT / "release/lean_verification/butson-towers-lean-repro.zip").read_bytes())
    assert proof_verification == (ROOT / "release/lean_verification/verification.json").read_bytes()
    files["anc/proof-sources.json"] = (HERE / "proof-sources.json").read_bytes()
    files["anc/lean-verification.json"] = proof_verification
    for certificate_name, script_name in (
        ("williamson-spence-family.json", "williamson_spence_family.py"),
        ("williamson-cyclotomic-trace-audit.json", "williamson_cyclotomic_trace_audit.py"),
    ):
        assert json.loads(files["anc/" + certificate_name])["implementation_sha256"] == digest(files["anc/" + script_name])
    new_orders = {7572, 17668, 22716, 32812, 98436, 295308}
    known_orders = {1524, 3252, 6564, 885924}
    for report in (weighted, weighted_control, families):
        for item in report["outputs"]:
            status = item["reported_order_status"]["status"]
            assert status in {"known", "new_to_best_of_authors_knowledge"}
            (known_orders if status == "known" else new_orders).add(item["output_order"])
    assert len(new_orders) == 18 and len(known_orders) == 20 and not new_orders & known_orders
    ancillary_manifest = {
        "purpose": "Exact Hadamard scalar restrictions, six new output families and periodic T127 reconstruction",
        "proof_assistant_verification": {
            "verified_ranks": "all eighteen announced ranks",
            "verified_all_parameter_families": 6,
            "finite_witnesses_and_construction_proofs": True,
            "audit_record": "anc/lean-verification.json",
            "source_packet_pointer": "anc/proof-sources.json",
            "Lean_source_packet_included_in_this_archive": False,
        },
        "novelty_convention": "Known where an earlier construction was located; otherwise new to the best of the author's knowledge.",
        "reported_results": {
            "new_real_Hadamard_orders": sorted(new_orders),
            "new_output_family_count": 6,
            "family_novelty_scope": "The six output construction families are new to the best of the author's knowledge; finite ranks have individual known/new labels.",
            "new_periodic_T_matrix_orders": [127],
            "known_real_Hadamard_orders": sorted(known_orders),
            "infinite_real_order_families": [item["real_endpoint_order_formula"] for item in families["family_source_proof"]["families"]] + ["32812*3^r"],
            "classical_source_recurrences": {"Spence": 5, "unreal_tripling": 1},
            "literature_search_date": "2026-11-24",
        },
        "original_novelty_audits": "Retained as dated provenance; current author claims are recorded in reported_results.",
        "replay": "python3 anc/recipe.py --write h7572-reconstructed.npz",
        "further_endpoints_replay": "python3 anc/reuse.py",
        "periodic_T127_replay": "python3 anc/t127/recipe.py --write-hadamard h1524-reconstructed.npz",
        "weighted_scalar_replay": "python3 anc/williamson_scalar_endpoints.py --output weighted-replay.json",
        "infinite_families_replay": "python3 anc/williamson_scalar_endpoints.py --families --output families-replay.json",
        "optional_independent_trace_audit": "python3 anc/williamson_cyclotomic_trace_audit.py --output trace-replay.json",
        "T127_original_certificate_sha256": digest((t127 / "certificate.json").read_bytes()),
        "T127_omitted_reference_files": ["README.md"],
        "dependencies": ["Python>=3.10", "NumPy"],
        "optional_independent_trace_audit_additional_dependency": "SymPy",
        "reference_packed_sha256": certificate["matrix_storage"]["packed_sign_bits_sha256"],
        "original_certificate_sha256": digest((BUNDLE / "certificate.json").read_bytes()),
        "omitted_reference_bundle_files": ["README.md", "h7572.npz"],
        "omission_reason": "The archive-free replay rebuilds the matrix from the included sequence data; the original certificate retains reference hashes for its complete repository bundle.",
        "included_files": {name: digest(data) for name, data in files.items() if name.startswith("anc/")},
        "server_execution_required": False,
    }
    files["anc/manifest.json"] = (json.dumps(ancillary_manifest, indent=2) + "\n").encode()
    credential = ROOT / "hexagon_key.local"
    secret_values = []
    if credential.exists():
        raw = credential.read_text().strip()
        secret_values.append(raw.encode())
        match = re.fullmatch(r"(?:HEXAGON_API_KEY|HEXAGON_TOKEN)\s*=\s*(.+)", raw)
        if match:
            raw = match.group(1).strip().strip("\"'")
        secret_values.append(raw.removeprefix("Bearer ").encode())
    for name, data in files.items():
        assert not name.startswith("/") and ".." not in Path(name).parts
        assert ".local" not in Path(name).parts and not name.endswith(".local")
        assert len(data) <= 4 * 1024 * 1024
        assert not any(secret and secret in data for secret in secret_values), "Credential material in payload"
    archive = HERE / "source.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as out:
        for name, data in files.items():
            info = zipfile.ZipInfo(name, (2026, 11, 24, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            out.writestr(info, data)
    with zipfile.ZipFile(archive) as saved:
        assert saved.testzip() is None and set(saved.namelist()) == set(files)
        assert all(saved.read(name) == data for name, data in files.items())
    metadata_path = HERE / "metadata.template.json"
    metadata = apply_administrative_defaults(json.loads(metadata_path.read_text()))
    write_json(metadata_path, metadata)
    contract = json.loads((HERE / "api-contract.json").read_text())
    assert set(contract["metadata_required_fields"]) == set(REQUIRED)
    missing = [name for name in REQUIRED if metadata.get(name) is None]
    invalid = []
    enums = {
        "license": {"CC BY 4.0", "CC BY-SA 4.0", "CC0 1.0"},
        "authorshipMode": {"no-human-author-claimed", "mixed", "human"},
        "humanUnderstanding": {"all", "some", "none"},
    }
    for name, allowed in enums.items():
        if metadata.get(name) is not None and metadata[name] not in allowed:
            invalid.append(name)
    for name in ("correctnessConfirmed", "communityInterestConfirmed", "citationsConfirmed"):
        if metadata.get(name) is not None and metadata[name] is not True:
            invalid.append(name)
    if metadata.get("rightsConfirmed") is not None and type(metadata["rightsConfirmed"]) is not bool:
        invalid.append("rightsConfirmed")
    if metadata.get("humanUnderstanding") == "some" and not metadata.get("provenance"):
        invalid.append("provenance")
    ready = not missing and not invalid
    write_json(HERE / "submission-manifest.json", {
        "schema": "hexagon-butson-towers-source-preparation-v2", "origin": ORIGIN,
        "release": "butson_towers_v2", "supersedes_local_release": "hadamard_7572_announcement",
        "source_archive": "source.zip", "source_archive_sha256": digest(archive.read_bytes()),
        "source_archive_members": [{"path": name, "bytes": len(data), "sha256": digest(data)} for name, data in files.items()],
        "upload_request": {"kind": "source-archive", "filename": "source.zip", "relativePath": "source.zip", "contentType": "application/zip", "sizeBytes": archive.stat().st_size},
        "entry_point_after_expansion": "main.tex", "reference_pdf": "announcement.pdf",
        "reference_pdf_sha256": build["pdf_sha256"], "pages": build["pages"],
        "metadata_file": "metadata.template.json", "metadata_sha256": digest((HERE / "metadata.template.json").read_bytes()),
        "api_contract": "api-contract.json", "api_schema_sha256": contract["schema_sha256"],
        "api_schema_version": contract["schema_version"],
        "ancillary_extensions_used": sorted({Path(name).suffix for name in files if name.startswith("anc/")}),
        "formal_proof_packet": proof_sources,
        "credential_included": False, "copyrighted_primary_source_copies_included": False,
    })
    write_json(HERE / "readiness.json", {
        "ready_to_create_private_draft": ready,
        "ready_for_submission": ready and metadata.get("rightsConfirmed") is True,
        "missing_metadata_fields": missing, "invalid_metadata_fields": invalid,
        "paper_pages": build["pages"], "local_typesetting_passed": True,
        "credential_excluded_from_archive": True, "remote_draft_created": False,
        "api_commit_performed": False,
    })
    if ready:
        write_json(HERE / "draft-request.json", metadata)
        request_bytes = json.dumps(metadata, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
        write_json(HERE / "draft-headers.json", {
            "Idempotency-Key": "butson-towers-v2-" + digest(request_bytes + archive.read_bytes()),
            "Content-Type": "application/json",
        })
    elif (HERE / "draft-request.json").exists():
        (HERE / "draft-request.json").unlink()
    if not ready and (HERE / "draft-headers.json").exists():
        (HERE / "draft-headers.json").unlink()
    print(json.dumps({"archive": str(archive.relative_to(ROOT)), "members": len(files),
                      "ready_to_create_private_draft": ready, "missing_metadata_fields": missing}))


if __name__ == "__main__":
    main()
