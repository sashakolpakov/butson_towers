#!/usr/bin/env python3
"""Verify the pinned, partial Butson Towers Lean construction kernel."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

TOOLCHAIN = "leanprover/lean4:v4.32.0"
MATHLIB_REVISION = "81a5d257c8e410db227a6665ed08f64fea08e997"
ALLOWED_AXIOMS = frozenset({"propext", "Classical.choice", "Quot.sound"})
THEOREMS = (
    "ButsonTowers.inverse_metric_adjoint", "ButsonTowers.framed_zero",
    "ButsonTowers.frame_pair_transport", "ButsonTowers.frame_sum_transport",
    "ButsonTowers.weighted_pair_transport", "ButsonTowers.weighted_sum_transport",
    "ButsonTowers.weighted_gram_block", "ButsonTowers.frame_gram_block",
    "ButsonTowers.framed_matrix_flat", "ButsonTowers.framed_matrix_gram",
    "ButsonTowers.framed_matrix_hadamard",
    "ButsonTowers.Spence.doubled_minus_norm", "ButsonTowers.Spence.doubled_plus_norm",
    "ButsonTowers.Spence.sum_squares_from_strong_norms",
    "ButsonTowers.Spence.strong_norm_induction", "ButsonTowers.Spence.entry_step_sign",
)
SOURCE_FILES = (
    "formal/.gitignore", "formal/lean-toolchain", "formal/lakefile.toml",
    "formal/lake-manifest.json", "formal/ButsonTowers.lean", "formal/AuditAxioms.lean",
    "formal/ButsonTowers/ScalarFrame.lean", "formal/ButsonTowers/Spence.lean",
    "docs/FORMALIZATION.md",
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
    audit = (formal / "AuditAxioms.lean").read_text()
    require(tuple(re.findall(r"^#print axioms (\S+)\s*$", audit, re.MULTILINE)) == THEOREMS,
            "axiom audit does not name exactly the expected sixteen theorems")
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
    for line in output.splitlines():
        if not line.strip():
            continue
        match = re.fullmatch(r"'([^']+)' depends on axioms: \[([^\]]*)\]", line.strip())
        empty = re.fullmatch(r"'([^']+)' does not depend on any axioms", line.strip())
        require(match is not None or empty is not None, "unexpected axiom-audit output")
        theorem = (match or empty).group(1)
        require(theorem in THEOREMS and theorem not in reports, "unexpected or duplicate audited theorem")
        axioms = set(filter(None, (s.strip() for s in match.group(2).split(",")))) if match else set()
        require(axioms <= ALLOWED_AXIOMS, f"unsupported axiom in {theorem}")
        reports[theorem] = sorted(axioms)
    require(set(reports) == set(THEOREMS), "missing theorem axiom reports")
    return reports


def self_test():
    valid = "\n".join(f"'{name}' depends on axioms: [propext, Quot.sound]" for name in THEOREMS)
    require(len(parse_axioms(valid)) == 16, "valid parser fixture failed")
    for invalid in ("\n".join(valid.splitlines()[:-1]), valid + "\n" + valid.splitlines()[0],
                    valid.replace("[propext, Quot.sound]", "[sorryAx]", 1),
                    valid.replace(THEOREMS[0], "Unrelated.theorem", 1)):
        try:
            parse_axioms(invalid)
        except VerificationError:
            continue
        raise VerificationError("invalid parser fixture accepted")
    print("Axiom-parser positive, missing, duplicate, foreign-axiom and foreign-theorem checks passed.")


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
    build_command = [lake] + (["--no-cache"] if args.from_source or args.cache_dir else []) + ["build"]
    run(build_command, formal)
    audited = parse_axioms(run([lake, "env", "lean", "AuditAxioms.lean"], formal))
    verify_payload(root)  # The build must not rewrite the portable lock or proof sources.
    check_dependencies(formal / ".lake" / "packages", packages)
    report = {"schema": "verified-butson-towers-lean-kernel-v1", "success": True,
              "payload_manifest_sha256": digest((root / "MANIFEST.json").read_bytes()),
              "scope": "sixteen general construction theorems; concrete matrices and all endpoint claims are not yet formally instantiated",
              "lean_toolchain": TOOLCHAIN, "lake_version": version,
              "mathlib_revision": MATHLIB_REVISION, "dependency_pins": dependency_report,
              "proof_theorems_checked": len(audited), "axioms_by_theorem": audited,
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
