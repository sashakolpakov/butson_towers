#!/usr/bin/env python3
"""Check this release's hashes, metadata and extracted standalone Python recipes."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import zipfile

from build import prose, title_argument

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    start = time.monotonic()
    archive = HERE / "source.zip"
    manifest = json.loads((HERE / "submission-manifest.json").read_text())
    metadata = json.loads((HERE / "metadata.template.json").read_text())
    build = json.loads((HERE / "build.json").read_text())
    require(digest(archive.read_bytes()) == manifest["source_archive_sha256"], "archive hash mismatch")
    require(digest((HERE / "announcement.pdf").read_bytes()) == build["pdf_sha256"]
            == manifest["reference_pdf_sha256"], "PDF hash mismatch")
    require(digest((HERE / "metadata.template.json").read_bytes()) == manifest["metadata_sha256"],
            "metadata hash mismatch")
    outputs = {}
    replay_env = os.environ.copy()
    replay_env.pop("PYTHONPATH", None)
    replay_env["PYTHONNOUSERSITE"] = "1"
    with zipfile.ZipFile(archive) as saved, tempfile.TemporaryDirectory(prefix="butson-v3-replay-") as temp:
        records = manifest["source_archive_members"]
        require(len(saved.namelist()) == len(records)
                and set(saved.namelist()) == {r["path"] for r in records}, "archive members differ")
        for record in records:
            name = record["path"]
            require(not Path(name).is_absolute() and ".." not in Path(name).parts, "unsafe archive path")
            require(not any(p.casefold().endswith(".local") for p in Path(name).parts), "private archive path")
            data = saved.read(name)
            require(len(data) == record["bytes"] and digest(data) == record["sha256"], "member hash mismatch")
        source = saved.read("main.tex")
        require(source == (HERE / "main.tex").read_bytes()
                == (ROOT / "manuscript/hadamard_7572_announcement.tex").read_bytes(), "source mismatch")
        text = source.decode()
        require(prose(title_argument(text)) == metadata["title"], "title mismatch")
        abstract = prose(text.split(r"\begin{abstract}", 1)[1].split(r"\end{abstract}", 1)[0])
        require(abstract == metadata["abstract"], "abstract mismatch")
        payload = json.loads(saved.read("anc/manifest.json"))
        for name, expected in payload["included_files"].items():
            require(digest(saved.read(name)) == expected, "ancillary manifest mismatch")
        expected_orders = payload["reported_results"]["new_real_Hadamard_orders"]
        require(len(expected_orders) == 18, "announced rank count mismatch")
        require(len(payload["reported_results"]["known_real_Hadamard_orders"]) == 20, "known rank count mismatch")
        require(len(payload["reported_results"]["infinite_real_order_families"]) == 6, "family count mismatch")
        proof_sources = json.loads(saved.read("anc/proof-sources.json"))
        proof_report = json.loads(saved.read("anc/lean-verification.json"))
        require(proof_report["success"] and proof_report["proof_theorems_checked"] == 85,
                "Lean verification record incomplete")
        require(proof_sources["source_package_sha256"] == digest(
            (ROOT / "release/lean_verification/butson-towers-lean-repro.zip").read_bytes()),
            "separate Lean packet hash mismatch")
        fresh = json.loads(saved.read("anc/lean-rerun.json"))
        pair = json.loads(saved.read("anc/sixth-root-verification.json"))
        pair_pointer = json.loads(saved.read("anc/sixth-root-proof-sources.json"))
        require(fresh["success"] and fresh["proof_theorems_checked"] == 85,
                "fresh baseline kernel audit missing")
        require(pair["success"] and pair["proof_theorems_checked"] == 10,
                "sixth-root supplement audit missing")
        require(pair_pointer["source_package_sha256"] == digest(
            (ROOT / "release/lean_sixth_root_pair/butson-towers-sixth-root-pair.zip").read_bytes()),
            "sixth-root supplement hash mismatch")
        require(saved.read("anc/manuscript_data.py") == (ROOT / "scripts/verify_manuscript_data.py").read_bytes(),
                "manuscript checker source differs")
        saved.extractall(temp)
        commands = (
            ("manuscript_data", ["anc/manuscript_data.py"]),
            ("h7572", ["anc/recipe.py", "--write", "h7572-reconstructed.npz"]),
            ("frame_reuse", ["anc/reuse.py"]),
            ("periodic_t127", ["anc/t127/recipe.py", "--write-hadamard", "h1524-reconstructed.npz"]),
            ("weighted", ["anc/williamson_scalar_endpoints.py", "--output", "weighted-replay.json"]),
            ("families", ["anc/williamson_scalar_endpoints.py", "--families", "--output", "families-replay.json"]),
            ("trace", ["anc/williamson_cyclotomic_trace_audit.py", "--output", "trace-replay.json"]),
        )
        for name, command in commands:
            result = subprocess.run([sys.executable, *command], cwd=temp, text=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=replay_env)
            require(result.returncode == 0, f"standalone {name} replay failed: {result.stdout}")
            outputs[name] = json.loads(result.stdout)
            print(f"Extracted {name} replay passed", flush=True)
        require(outputs["manuscript_data"]["all_checks_passed"],
                "manuscript finite calculations incomplete")
        require(len(outputs["weighted"]["orders"]) == 27 and outputs["weighted"]["complete_Gram_proofs"],
                "weighted source replay incomplete")
        require(len(outputs["families"]["orders"]) == 5 and outputs["families"]["complete_Gram_proofs"],
                "family replay incomplete")
        require(outputs["trace"]["all_exact_checks_passed"], "independent trace audit incomplete")
        require(outputs["frame_reuse"]["full_P_Q_and_symmetric_cross_integer_identities_verified"]
                and len(outputs["frame_reuse"]["outputs"]) == 6, "frame/cubic replay incomplete")
        require(outputs["h7572"]["full_endpoint_integer_Gram_verified"], "H7572 Gram incomplete")
        require(outputs["h7572"]["endpoint_row_pairs_checked"] == 28663806, "H7572 row-pair count mismatch")
        require(not outputs["h7572"]["reference_archive_checked"], "replay used a repository NPZ")
        require(outputs["periodic_t127"]["scalar_Hadamard_control"]["full_H1524_integer_Gram_verified"],
                "T127 control Gram incomplete")
    report = {
        "schema": "butson-towers-announcement-extracted-replay-v3",
        "source_zip_sha256": digest(archive.read_bytes()), "source_and_metadata_synchronized": True,
        "all_archive_members_verified": True, "all_replays_used_only_extracted_archive_files": True,
        "standalone_replay_environment": {
            "PYTHONPATH_removed": True,
            "PYTHONNOUSERSITE": "1",
            "working_directory": "fresh extracted archive temporary directory",
        },
        "full_H7572_integer_Gram_verified": True, "H7572_row_pairs_checked": 28663806,
        "full_H1524_integer_Gram_verified": True,
        "manuscript_finite_calculations_replayed": True,
        "frame_reuse_and_cubic_source_replays_passed": True,
        "all_27_weighted_sources_replayed": True, "five_all_parameter_family_recipes_replayed": True,
        "independent_cyclotomic_trace_audit_passed": True,
        "new_real_Hadamard_ranks": expected_orders,
        "known_real_Hadamard_ranks": payload["reported_results"]["known_real_Hadamard_orders"],
        "six_new_output_families": payload["reported_results"]["infinite_real_order_families"],
        "unchanged_separate_Lean_packet_hash_verified": True,
        "Lean_kernel_rerun_for_this_editorial_release": True,
        "included_completed_Lean_audit_theorems": proof_report["proof_theorems_checked"],
        "fresh_baseline_kernel_theorems": 85,
        "additional_sixth_root_matrix_kernel_theorems": 10,
        "remote_requests_performed": False, "elapsed_seconds": round(time.monotonic() - start, 2),
    }
    (HERE / "archive-verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"source_archive_independent_replay_passed": True,
                      "new_ranks": len(expected_orders), "families": 6}))


if __name__ == "__main__":
    main()
