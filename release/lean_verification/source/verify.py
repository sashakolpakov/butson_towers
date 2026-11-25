#!/usr/bin/env python3
"""Verify the pinned Butson Towers Lean existence proofs."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path
import re
import shutil
import subprocess
import sys

TOOLCHAIN = "leanprover/lean4:v4.32.0"
MATHLIB_REVISION = "81a5d257c8e410db227a6665ed08f64fea08e997"
ALLOWED_AXIOMS = frozenset({"propext", "Classical.choice", "Quot.sound"})
THEOREMS = (
    'ButsonTowers.inverse_metric_adjoint',
    'ButsonTowers.framed_zero',
    'ButsonTowers.frame_pair_transport',
    'ButsonTowers.frame_sum_transport',
    'ButsonTowers.weighted_pair_transport',
    'ButsonTowers.weighted_sum_transport',
    'ButsonTowers.weighted_gram_block',
    'ButsonTowers.frame_gram_block',
    'ButsonTowers.framed_matrix_flat',
    'ButsonTowers.framed_matrix_gram',
    'ButsonTowers.framed_matrix_hadamard',
    'ButsonTowers.Spence.doubled_minus_norm',
    'ButsonTowers.Spence.doubled_plus_norm',
    'ButsonTowers.Spence.sum_squares_from_strong_norms',
    'ButsonTowers.Spence.strong_norm_induction',
    'ButsonTowers.Spence.entry_step_sign',
    'ButsonTowers.TraceTowers.trace_metric_adjoint',
    'ButsonTowers.TraceTowers.inverse_trace_metric_adjoint',
    'ButsonTowers.TraceTowers.weighted_trace_gram_transport',
    'ButsonTowers.TraceTowers.basis_gram_congruence',
    'ButsonTowers.GroupAlgebra.regular_mul',
    'ButsonTowers.GroupAlgebra.regular_involution',
    'ButsonTowers.GroupAlgebra.regular_trace',
    'ButsonTowers.GroupAlgebra.inverse_basis_trace_metric',
    'ButsonTowers.GroupAlgebra.regularAlgHom_injective',
    'ButsonTowers.GroupAlgebra.regular_star',
    'ButsonTowers.GroupAlgebra.regularAlgHom_normalized_trace',
    'ButsonTowers.GroupAlgebra.group_algebra_identity_trace_metric',
    'ButsonTowers.GoethalsSeidel.gs_gram',
    'ButsonTowers.OrthogonalDesign.coefficient_gram',
    'ButsonTowers.OrthogonalDesign.coefficient_cross',
    'ButsonTowers.OrthogonalDesign.coefficients_signed_partition',
    'ButsonTowers.tensor_substitute_gram',
    'ButsonTowers.tensor_substitute_flat',
    'ButsonTowers.tensor_substitute_hadamard',
    'ButsonTowers.SpenceMatrix.strong_step',
    'ButsonTowers.SpenceMatrix.family_strong',
    'ButsonTowers.SpenceMatrix.family_norm',
    'ButsonTowers.SpenceMatrix.family_card',
    'ButsonTowers.SpenceMatrix.family_hadamard',
    'ButsonTowers.Unreal.next_amicable',
    'ButsonTowers.Unreal.next_norm',
    'ButsonTowers.Unreal.quad_amicable',
    'ButsonTowers.hadamard_reindex',
    'ButsonTowers.exists_from_t',
    'ButsonTowers.spence_exists_from_t',
    'ButsonTowers.Data.t127_codes_corr',
    'ButsonTowers.Data.t127_corr',
    'ButsonTowers.Data.t127_onehot',
    'ButsonTowers.Data.t631_codes_corr',
    'ButsonTowers.Data.t631_corr',
    'ButsonTowers.Data.t631_onehot',
    'ButsonTowers.Data.unreal13_norm',
    'ButsonTowers.hadamard_7572',
    'ButsonTowers.hadamard_17668',
    'ButsonTowers.hadamard_22716',
    'ButsonTowers.hadamard_32812',
    'ButsonTowers.hadamard_42908',
    'ButsonTowers.hadamard_58052',
    'ButsonTowers.hadamard_73196',
    'ButsonTowers.hadamard_78244',
    'ButsonTowers.hadamard_93388',
    'ButsonTowers.hadamard_98436',
    'ButsonTowers.hadamard_108532',
    'ButsonTowers.hadamard_128724',
    'ButsonTowers.hadamard_153964',
    'ButsonTowers.hadamard_229684',
    'ButsonTowers.hadamard_295308',
    'ButsonTowers.hadamard_386172',
    'ButsonTowers.hadamard_4772884',
    'ButsonTowers.hadamard_13001124',
    'ButsonTowers.spence_family_13',
    'ButsonTowers.spence_family_17',
    'ButsonTowers.spence_family_61',
    'ButsonTowers.spence_family_101',
    'ButsonTowers.spence_family_121',
    'ButsonTowers.unreal_cubic_family',
    'ButsonTowers.announced_orders_exist',
    'ButsonTowers.williamson_table_endpoints',
    'ButsonTowers.GoethalsSeidel.python_reflection_127',
    'ButsonTowers.GoethalsSeidel.python_reflection_631',
    'ButsonTowers.GSAlgebra.blocks_gram',
    'ButsonTowers.TraceTowers.finite_free_trace_transitivity',
    'ButsonTowers.TraceTowers.tower_action_coordinates',
    'ButsonTowers.TraceTowers.tower_trace_metric',
)
SOURCE_FILES = (
    'formal/.gitignore',
    'formal/lean-toolchain',
    'formal/lakefile.toml',
    'formal/lake-manifest.json',
    'formal/ButsonTowers.lean',
    'formal/AuditAxioms.lean',
    'formal/build.py',
    'formal/DATA_PROVENANCE.json',
    'formal/ButsonTowers/Construction.lean',
    'formal/ButsonTowers/Data/CorrRange.lean',
    'formal/ButsonTowers/Data/FiniteChecks.lean',
    'formal/ButsonTowers/Data/LiteralTree.lean',
    'formal/ButsonTowers/Data/Paley.lean',
    'formal/ButsonTowers/Data/PaleyChecks.lean',
    'formal/ButsonTowers/Data/PaleyPairChecks.lean',
    'formal/ButsonTowers/Data/RowChecks.lean',
    'formal/ButsonTowers/Data/StrongChecks.lean',
    'formal/ButsonTowers/Data/StrongStarters.lean',
    'formal/ButsonTowers/Data/T127.lean',
    'formal/ButsonTowers/Data/T631.lean',
    'formal/ButsonTowers/Data/T631Blocked.lean',
    'formal/ButsonTowers/Data/T631Blocks/Chunk0.lean',
    'formal/ButsonTowers/Data/T631Blocks/Chunk1.lean',
    'formal/ButsonTowers/Data/T631Blocks/Chunk2.lean',
    'formal/ButsonTowers/Data/T631Blocks/Chunk3.lean',
    'formal/ButsonTowers/Data/TreeLiterals.lean',
    'formal/ButsonTowers/Data/Unreal13.lean',
    'formal/ButsonTowers/Data/Williamson/Chunk0.lean',
    'formal/ButsonTowers/Data/Williamson/Chunk1.lean',
    'formal/ButsonTowers/Data/Williamson/Chunk2.lean',
    'formal/ButsonTowers/Data/Williamson.lean',
    'formal/ButsonTowers/Endpoints.lean',
    'formal/ButsonTowers/Existence.lean',
    'formal/ButsonTowers/FiniteGroupAlgebra.lean',
    'formal/ButsonTowers/GSAlgebra.lean',
    'formal/ButsonTowers/GoethalsSeidel.lean',
    'formal/ButsonTowers/GoethalsSeidelCore.lean',
    'formal/ButsonTowers/GoethalsSeidelEntry.lean',
    'formal/ButsonTowers/GroupAlgebra.lean',
    'formal/ButsonTowers/OrthogonalDesign.lean',
    'formal/ButsonTowers/OrthogonalDesignCore.lean',
    'formal/ButsonTowers/OrthogonalSubstitution.lean',
    'formal/ButsonTowers/ReflectionConvention.lean',
    'formal/ButsonTowers/ScalarFrame.lean',
    'formal/ButsonTowers/Spence.lean',
    'formal/ButsonTowers/SpenceHadamard.lean',
    'formal/ButsonTowers/SpenceMatrix.lean',
    'formal/ButsonTowers/TraceTowers.lean',
    'formal/ButsonTowers/UnrealTripling.lean',
    'docs/FORMALIZATION.md',
)
SUPPORT_FILES = ("README.md", "verify.py", "build_package.py")
PAYLOAD_FILES = SOURCE_FILES + SUPPORT_FILES


class VerificationError(RuntimeError):
    pass


def require(condition, message):
    if not condition:
        raise VerificationError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def regular_file(root, relative):
    require(relative in PAYLOAD_FILES or relative == "MANIFEST.json", "unexpected payload path")
    path = root / relative
    require(not any(part == ".local" or part.endswith(".local") for part in Path(relative).parts),
            "local configuration cannot be packaged")
    require(all(not p.is_symlink() for p in (path, *path.parents) if p != root.parent),
            f"symlink in payload path: {relative}")
    require(path.is_file(), f"missing regular payload file: {relative}")
    return path


def verify_payload(root):
    manifest = json.loads(regular_file(root, "MANIFEST.json").read_text())
    require(manifest["schema"] == "butson-towers-lean-reproducibility-v1", "unexpected manifest schema")
    require(manifest["lean_toolchain"] == TOOLCHAIN, "manifest toolchain mismatch")
    require(manifest["mathlib_revision"] == MATHLIB_REVISION, "manifest mathlib mismatch")
    require(manifest["expected_theorems"] == list(THEOREMS), "manifest theorem set mismatch")
    records = manifest["files"]
    require(len(records) == len(PAYLOAD_FILES)
            and {r["path"] for r in records} == set(PAYLOAD_FILES), "payload whitelist mismatch")
    for record in records:
        data = regular_file(root, record["path"]).read_bytes()
        require(len(data) == record["bytes"] and digest(data) == record["sha256"],
                f"payload changed: {record['path']}")
    formal = root / "formal"
    require((formal / "lean-toolchain").read_text().strip() == TOOLCHAIN, "toolchain file mismatch")
    locked = json.loads((formal / "lake-manifest.json").read_text())
    packages = locked["packages"]
    require(len({p["name"] for p in packages}) == len(packages), "duplicate dependency name")
    require(all(p["type"] == "git" and re.fullmatch(r"[0-9a-f]{40}", p["rev"])
                and re.fullmatch(r"[A-Za-z0-9_]+", p["name"]) for p in packages),
            "every dependency must have a locked Git revision and a simple package name")
    require(next(p for p in packages if p["name"] == "mathlib")["rev"] == MATHLIB_REVISION,
            "locked mathlib revision mismatch")
    expected_lean = {str(Path(p).relative_to("formal")) for p in SOURCE_FILES
                     if p.startswith("formal/") and p.endswith(".lean")}
    actual_lean = set()
    for directory, subdirs, files in os.walk(formal, followlinks=False):
        subdirs[:] = [name for name in subdirs if name != ".lake"]
        require(not any((Path(directory) / name).is_symlink() for name in subdirs),
                "unexpected symlinked formal source directory")
        actual_lean.update(str((Path(directory) / name).relative_to(formal))
                           for name in files if name.endswith(".lean"))
    require(actual_lean == expected_lean, "formal Lean source set differs from the whitelist")
    audit = (formal / "AuditAxioms.lean").read_text()
    require(tuple(re.findall(r"^#print axioms (\S+)\s*$", audit, re.MULTILINE)) == THEOREMS,
            "axiom audit does not name exactly the expected theorems")
    for relative in SOURCE_FILES:
        if relative.endswith(".lean"):
            text = (root / relative).read_text()
            require(not re.search(r"\b(sorry|admit|axiom|unsafe|native_decide)\b", text),
                    f"placeholder or unsupported proof feature in {relative}")
    return packages


def run(command, cwd, *, echo=True):
    if echo:
        print("+ " + " ".join(map(str, command)), flush=True)
    result = subprocess.run(list(map(str, command)), cwd=cwd, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if echo and result.stdout:
        print(result.stdout, end="" if result.stdout.endswith("\n") else "\n", flush=True)
    require(result.returncode == 0, f"command failed ({result.returncode}): {command[0]}")
    return result.stdout


def check_dependencies(directory, packages):
    records = []
    for package in packages:
        path = directory / package["name"]
        require(path.is_dir(), f"missing dependency: {package['name']}")
        head = run(["git", "-C", path, "rev-parse", "HEAD"], directory, echo=False).strip()
        require(head == package["rev"], f"Git HEAD does not match lock: {package['name']}")
        run(["git", "-C", path, "diff", "--quiet", "HEAD", "--"], directory, echo=False)
        records.append({"name": package["name"], "revision": head, "tracked_tree_clean": True})
    return records


def link_cache(formal, directory, packages):
    # Validate every source checkout before creating a single dependency link.
    records = check_dependencies(directory, packages)
    target = formal / ".lake" / "packages"
    target.mkdir(parents=True, exist_ok=True)
    for package in packages:
        path, source = target / package["name"], directory / package["name"]
        if path.exists() or path.is_symlink():
            require(path.is_symlink() and path.resolve() == source.resolve(),
                    f"dependency destination already exists: {package['name']}")
        else:
            path.symlink_to(source.resolve(), target_is_directory=True)
    return records


def parse_axioms(output):
    reports = {}
    pattern = re.compile(
        r"'([^']+)'\s+(?:depends\s+on\s+axioms:\s*\[([^\]]*)\]"
        r"|does\s+not\s+depend\s+on\s+any\s+axioms)")
    position = 0
    for match in pattern.finditer(output):
        require(not output[position:match.start()].strip(), "unexpected axiom-audit output")
        theorem = match.group(1)
        require(theorem in THEOREMS and theorem not in reports, "unexpected or duplicate audited theorem")
        axioms = set(filter(None, (s.strip() for s in (match.group(2) or "").split(","))))
        require(axioms <= ALLOWED_AXIOMS, f"unsupported axiom in {theorem}")
        reports[theorem] = sorted(axioms)
        position = match.end()
    require(not output[position:].strip(), "unexpected trailing axiom-audit output")
    require(set(reports) == set(THEOREMS), "missing theorem axiom reports")
    return reports


def self_test():
    valid = "\n".join(f"'{name}' depends on axioms: [propext, Quot.sound]" for name in THEOREMS)
    require(len(parse_axioms(valid)) == len(THEOREMS), "valid parser fixture failed")
    wrapped = valid.replace("[propext, Quot.sound]", "[\n  propext,\n  Quot.sound\n]")
    require(parse_axioms(wrapped) == parse_axioms(valid), "wrapped parser fixture failed")
    for invalid in ("unexpected\n" + valid, valid + "\nunexpected",
                    "\n".join(valid.splitlines()[:-1]), valid + "\n" + valid.splitlines()[0],
                    valid.replace("[propext, Quot.sound]", "[sorryAx]", 1),
                    valid.replace(THEOREMS[0], "Unrelated.theorem", 1)):
        try:
            parse_axioms(invalid)
        except VerificationError:
            continue
        raise VerificationError("invalid parser fixture accepted")
    print("Axiom-parser single-line, wrapped, missing, duplicate, foreign-axiom, foreign-theorem and extra-text checks passed.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--cache-dir", type=Path, help="existing matching .lake/packages checkout directory")
    mode.add_argument("--from-source", action="store_true", help="skip prebuilt downloads; compile dependencies")
    parser.add_argument("--lake", default="lake", help="Lake executable; its Lean version is checked")
    parser.add_argument("--report", type=Path, help="optional output JSON report")
    parser.add_argument("--self-test", action="store_true", help="test the exact axiom-output parser only")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    started = time.monotonic()
    root = Path(__file__).resolve().parent
    packages = verify_payload(root)
    formal = root / "formal"
    lake = shutil.which(args.lake)
    require(lake is not None, "Lake is not on PATH; install the pinned Lean toolchain with elan")
    version = run([lake, "--version"], formal).strip()
    require(re.search(r"\(Lean version 4\.32\.0\)", version) is not None, "Lake uses the wrong Lean version")
    if args.cache_dir:
        link_cache(formal, args.cache_dir.resolve(), packages)
    elif not args.from_source:
        run([lake, "exe", "cache", "get"], formal)
    else:
        # Materializes the locked sources without updating the manifest.
        run([lake, "--no-cache", "env", "lean", "--version"], formal)
    dependency_report = check_dependencies(formal / ".lake" / "packages", packages)
    build_command = [sys.executable, formal / "build.py", "--lake", lake]
    if args.from_source or args.cache_dir:
        build_command.append("--no-cache")
    run(build_command, formal)
    audited = parse_axioms(run([lake, "env", "lean", "AuditAxioms.lean"], formal))
    verify_payload(root)  # The build must not rewrite the portable lock or proof sources.
    check_dependencies(formal / ".lake" / "packages", packages)
    report = {"schema": "verified-butson-towers-lean-kernel-v1", "success": True,
              "payload_manifest_sha256": digest((root / "MANIFEST.json").read_bytes()),
              "scope": "concrete T127/T631 witnesses, all 27 Williamson sources, eighteen announced orders and five all-iteration Spence families",
              "lean_toolchain": TOOLCHAIN, "lake_version": version,
              "mathlib_revision": MATHLIB_REVISION, "dependency_pins": dependency_report,
              "proof_theorems_checked": len(audited), "axioms_by_theorem": audited,
              "elapsed_seconds": round(time.monotonic() - started, 2),
              "project_axioms_or_placeholders_used": False,
              "dependency_mode": "existing-pinned-cache" if args.cache_dir else
                                 "source-build" if args.from_source else "official-cache-download"}
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    try:
        main()
    except (VerificationError, OSError, ValueError, KeyError) as error:
        print(f"Verification failed: {error}", file=sys.stderr)
        sys.exit(1)
