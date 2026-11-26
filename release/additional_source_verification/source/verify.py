#!/usr/bin/env python3
"""Verify the current additional sources and their complete canonical Lean closure."""
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
THEOREMS = ('ButsonTowers.Data.nearWilliamson47_quad',
 'ButsonTowers.Data.nearWilliamson53_quad',
 'ButsonTowers.Data.nearWilliamson59_quad',
 'ButsonTowers.Data.nearWilliamson67_quad',
 'ButsonTowers.Data.nearWilliamson69_quad',
 'ButsonTowers.Data.nearWilliamson73_quad',
 'ButsonTowers.Data.nearWilliamson79_quad',
 'ButsonTowers.Data.nearWilliamson93_quad',
 'ButsonTowers.Data.nearWilliamson97_quad',
 'ButsonTowers.Data.nearWilliamson103_quad',
 'ButsonTowers.Data.nearWilliamson109_quad',
 'ButsonTowers.Data.nearWilliamson127_quad',
 'ButsonTowers.hadamard_118628',
 'ButsonTowers.hadamard_133772',
 'ButsonTowers.hadamard_148916',
 'ButsonTowers.hadamard_169108',
 'ButsonTowers.hadamard_174156',
 'ButsonTowers.hadamard_184252',
 'ButsonTowers.hadamard_199396',
 'ButsonTowers.hadamard_234732',
 'ButsonTowers.hadamard_244828',
 'ButsonTowers.hadamard_259972',
 'ButsonTowers.hadamard_275116',
 'ButsonTowers.hadamard_320548',
 'ButsonTowers.Data.additionalSpence5_strong',
 'ButsonTowers.Data.additionalSpence53_strong',
 'ButsonTowers.Data.additionalSpence97_strong',
 'ButsonTowers.Data.additionalSpence109_strong',
 'ButsonTowers.Data.additionalSpence113_strong',
 'ButsonTowers.Data.additionalSpence137_strong',
 'ButsonTowers.Data.additionalPaley5_isPaley',
 'ButsonTowers.Data.additionalPaley53_isPaley',
 'ButsonTowers.Data.additionalPaley97_isPaley',
 'ButsonTowers.Data.additionalPaley109_isPaley',
 'ButsonTowers.Data.additionalPaley113_isPaley',
 'ButsonTowers.Data.additionalPaley137_isPaley',
 'ButsonTowers.spence_family_5',
 'ButsonTowers.spence_family_53',
 'ButsonTowers.spence_family_97',
 'ButsonTowers.spence_family_109',
 'ButsonTowers.spence_family_113',
 'ButsonTowers.spence_family_137',
 'ButsonTowers.near_williamson_table_endpoints',
 'ButsonTowers.hadamard_2985892',
 'ButsonTowers.hadamard_38816596',
 'ButsonTowers.hadamard_291145924',
 'ButsonTowers.hadamard_15131380',
 'ButsonTowers.new_spence_family_109')
LEGACY_THEOREMS = ('ButsonTowers.inverse_metric_adjoint',
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
 'ButsonTowers.TraceTowers.tower_trace_metric')
PAIR_THEOREMS = ('ButsonTowers.SixthRootPair.transpose_mul_of_mul_transpose',
 'ButsonTowers.SixthRootPair.action_mul_transpose',
 'ButsonTowers.SixthRootPair.action_add_transpose',
 'ButsonTowers.SixthRootPair.action_transpose',
 'ButsonTowers.SixthRootPair.action_polynomial',
 'ButsonTowers.SixthRootPair.action_mul_left',
 'ButsonTowers.SixthRootPair.imaginaryAction_transpose',
 'ButsonTowers.SixthRootPair.imaginaryAction_square',
 'ButsonTowers.SixthRootPair.coordinateAction_mul',
 'ButsonTowers.SixthRootPair.coordinateAction_conjugate')
AUDITS = {"AuditAxioms.lean": LEGACY_THEOREMS, "AuditSixthRootPair.lean": PAIR_THEOREMS, "AuditAdditionalSources.lean": THEOREMS}
SOURCE_FILES = ('formal/.gitignore',
 'formal/lean-toolchain',
 'formal/lakefile.toml',
 'formal/lake-manifest.json',
 'formal/build.py',
 'formal/DATA_PROVENANCE.json',
 'formal/NEW_SOURCE_VERIFICATION.md',
 'formal/AuditAdditionalSources.lean',
 'formal/AuditAxioms.lean',
 'formal/AuditSixthRootPair.lean',
 'formal/ButsonTowers.lean',
 'formal/ButsonTowers/AdditionalResults.lean',
 'formal/ButsonTowers/AdditionalSpenceFamilies.lean',
 'formal/ButsonTowers/Construction.lean',
 'formal/ButsonTowers/Data/AdditionalSpence.lean',
 'formal/ButsonTowers/Data/AdditionalSpence/Q109.lean',
 'formal/ButsonTowers/Data/AdditionalSpence/Q113.lean',
 'formal/ButsonTowers/Data/AdditionalSpence/Q137.lean',
 'formal/ButsonTowers/Data/AdditionalSpence/Q5.lean',
 'formal/ButsonTowers/Data/AdditionalSpence/Q53.lean',
 'formal/ButsonTowers/Data/AdditionalSpence/Q97.lean',
 'formal/ButsonTowers/Data/CorrRange.lean',
 'formal/ButsonTowers/Data/FiniteChecks.lean',
 'formal/ButsonTowers/Data/LiteralTree.lean',
 'formal/ButsonTowers/Data/NearWilliamson.lean',
 'formal/ButsonTowers/Data/NearWilliamson/Core.lean',
 'formal/ButsonTowers/Data/NearWilliamson/W103.lean',
 'formal/ButsonTowers/Data/NearWilliamson/W109.lean',
 'formal/ButsonTowers/Data/NearWilliamson/W127.lean',
 'formal/ButsonTowers/Data/NearWilliamson/W47.lean',
 'formal/ButsonTowers/Data/NearWilliamson/W53.lean',
 'formal/ButsonTowers/Data/NearWilliamson/W59.lean',
 'formal/ButsonTowers/Data/NearWilliamson/W67.lean',
 'formal/ButsonTowers/Data/NearWilliamson/W69.lean',
 'formal/ButsonTowers/Data/NearWilliamson/W73.lean',
 'formal/ButsonTowers/Data/NearWilliamson/W79.lean',
 'formal/ButsonTowers/Data/NearWilliamson/W93.lean',
 'formal/ButsonTowers/Data/NearWilliamson/W97.lean',
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
 'formal/ButsonTowers/Data/Williamson.lean',
 'formal/ButsonTowers/Data/Williamson/Chunk0.lean',
 'formal/ButsonTowers/Data/Williamson/Chunk1.lean',
 'formal/ButsonTowers/Data/Williamson/Chunk2.lean',
 'formal/ButsonTowers/Endpoints.lean',
 'formal/ButsonTowers/Existence.lean',
 'formal/ButsonTowers/FiniteGroupAlgebra.lean',
 'formal/ButsonTowers/GSAlgebra.lean',
 'formal/ButsonTowers/GoethalsSeidel.lean',
 'formal/ButsonTowers/GoethalsSeidelCore.lean',
 'formal/ButsonTowers/GoethalsSeidelEntry.lean',
 'formal/ButsonTowers/GroupAlgebra.lean',
 'formal/ButsonTowers/NearWilliamsonEndpoints.lean',
 'formal/ButsonTowers/OrthogonalDesign.lean',
 'formal/ButsonTowers/OrthogonalDesignCore.lean',
 'formal/ButsonTowers/OrthogonalSubstitution.lean',
 'formal/ButsonTowers/ReflectionConvention.lean',
 'formal/ButsonTowers/ScalarFrame.lean',
 'formal/ButsonTowers/SixthRootPair.lean',
 'formal/ButsonTowers/Spence.lean',
 'formal/ButsonTowers/SpenceHadamard.lean',
 'formal/ButsonTowers/SpenceMatrix.lean',
 'formal/ButsonTowers/TraceTowers.lean',
 'formal/ButsonTowers/UnrealTripling.lean',
 'scripts/check_lean_literal_data.py',
 'scripts/check_lean_additional_sources.py',
 'artifacts/t-matrices-127/t127.json',
 'artifacts/hadamard-7572/t631.json',
 'experiments/scalar_lift/williamson_inputs.json',
 'experiments/scalar_lift/near_williamson_inputs.json',
 'experiments/scalar_lift/bh13_6_seed.json',
 'certificates/scalar-lift/williamson-spence-family.json')
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
    require(manifest["schema"] == "butson-towers-additional-source-reproducibility-v1", "unexpected manifest schema")
    require(manifest["lean_toolchain"] == TOOLCHAIN, "manifest toolchain mismatch")
    require(manifest["mathlib_revision"] == MATHLIB_REVISION, "manifest mathlib mismatch")
    require(manifest["expected_theorems"] == list(THEOREMS), "manifest theorem set mismatch")
    require(manifest["expected_audits"] == {name: list(values) for name, values in AUDITS.items()},
            "manifest complete audit set mismatch")
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
    for filename, names in AUDITS.items():
        audit = (formal / filename).read_text()
        require(tuple(re.findall(r"^#print axioms (\S+)\s*$", audit, re.MULTILINE)) == names,
                f"axiom audit does not name exactly the expected theorems: {filename}")
    expected_files = set(PAYLOAD_FILES) | {"MANIFEST.json"}
    actual_files = set()
    for directory, subdirs, files in os.walk(root, followlinks=False):
        relative_dir = Path(directory).relative_to(root)
        subdirs[:] = [name for name in subdirs
                      if name != "__pycache__" and not (relative_dir == Path("formal") and name == ".lake")]
        require(not any((Path(directory) / name).is_symlink() for name in subdirs),
                "symlinked payload directory")
        actual_files.update(str((relative_dir / name)) for name in files)
    require(actual_files == expected_files, "payload file set differs from explicit whitelist")
    # Every canonical local import, including the three audits, must be present.
    modules = {".".join(Path(path).relative_to("formal").with_suffix("").parts)
               for path in SOURCE_FILES if path.startswith("formal/") and path.endswith(".lean")}
    for relative in SOURCE_FILES:
        if relative.endswith(".lean"):
            imports = re.findall(r"^import\s+([A-Za-z0-9_.]+)", (root / relative).read_text(), re.MULTILINE)
            require(all(not name.startswith("ButsonTowers") or name in modules for name in imports),
                    f"missing canonical local import in {relative}")
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


def parse_axioms(output, expected_theorems=THEOREMS):
    reports = {}
    pattern = re.compile(
        r"'([^']+)'\s+(?:depends\s+on\s+axioms:\s*\[([^\]]*)\]"
        r"|does\s+not\s+depend\s+on\s+any\s+axioms)")
    position = 0
    for match in pattern.finditer(output):
        require(not output[position:match.start()].strip(), "unexpected axiom-audit output")
        theorem = match.group(1)
        require(theorem in expected_theorems and theorem not in reports, "unexpected or duplicate audited theorem")
        axioms = set(filter(None, (s.strip() for s in (match.group(2) or "").split(","))))
        require(axioms <= ALLOWED_AXIOMS, f"unsupported axiom in {theorem}")
        reports[theorem] = sorted(axioms)
        position = match.end()
    require(not output[position:].strip(), "unexpected trailing axiom-audit output")
    require(set(reports) == set(expected_theorems), "missing theorem axiom reports")
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
    parser.add_argument("--check-only", action="store_true", help="check hashes, import closure and literal bindings without compiling Lean")
    parser.add_argument("--self-test", action="store_true", help="test the exact axiom-output parser only")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    started = time.monotonic()
    root = Path(__file__).resolve().parent
    packages = verify_payload(root)
    formal = root / "formal"
    for checker in ("check_lean_literal_data.py", "check_lean_additional_sources.py"):
        run([sys.executable, root / "scripts" / checker], root)
    if args.check_only:
        pins = check_dependencies(args.cache_dir.resolve(), packages) if args.cache_dir else None
        report = {"schema": "checked-additional-source-payload-v1", "success": True,
                  "kernel_checked": False, "scope": "payload hashes, exact source/import set and literal data bindings only",
                  "expected_additional_theorems": len(THEOREMS), "dependency_pins": pins}
        if args.report:
            args.report.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
        return
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
    audits = {filename: parse_axioms(run([lake, "env", "lean", filename], formal), names)
              for filename, names in AUDITS.items()}
    audited = audits["AuditAdditionalSources.lean"]
    verify_payload(root)  # The build must not rewrite the portable lock or proof sources.
    check_dependencies(formal / ".lake" / "packages", packages)
    report = {"schema": "verified-butson-towers-additional-source-kernel-v1", "success": True,
              "payload_manifest_sha256": digest((root / "MANIFEST.json").read_bytes()),
              "scope": "complete canonical proof closure; sixteen additional finite orders, six additional Spence recurrences and the shifted109 family",
              "lean_toolchain": TOOLCHAIN, "lake_version": version,
              "mathlib_revision": MATHLIB_REVISION, "dependency_pins": dependency_report,
              "proof_theorems_checked": sum(map(len, audits.values())),
              "additional_theorems_checked": len(audited), "axioms_by_audit": audits,
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
