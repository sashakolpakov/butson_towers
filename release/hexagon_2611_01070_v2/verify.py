#!/usr/bin/env python3
"""Check the v2 submission archive and replay only its extracted public files.

The Python replays prove finite identities. The pinned Lean report and source
packet establish a separate reproducible proof scope; this command does not
pretend to rerun the Lean kernel. Use --with-scalar-algebra for the optional
python-flint field-factor/CRT reconstruction audit.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import subprocess
import sys
import sysconfig
import tempfile
import time
import zipfile

HERE = Path(__file__).resolve().parent
MAX_MEMBER = 4*1024*1024
MAX_TOTAL = 100*1024*1024


def require(value, message):
    if not value:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def safe_member(name):
    path = PurePosixPath(name)
    require(isinstance(name, str) and str(path) == name and not path.is_absolute()
            and ".." not in path.parts and "\\" not in name and ":" not in name,
            "unsafe archive path")
    require(not any(part in {".git", ".lake", ".venv", "__pycache__", ".local"}
                    or part.casefold().endswith(".local") for part in path.parts),
            "private or generated archive path")
    require(name == "main.tex" or (name.startswith("anc/") and path.suffix in {".py", ".json"}),
            "unsupported ancillary path or extension")
    return path


def expected_members(payload):
    require(payload["schema"] == "butson-hexagon-v2-payload-v1", "unknown payload schema")
    repository = payload["repository_files"]
    require(len(repository) == len(set(repository)), "duplicate repository payload path")
    names = {"main.tex", "anc/manifest.json"}
    names.update("anc/"+name for name in repository)
    names.update(payload["release_files"].values())
    require(len(names) == 2+len(repository)+len(payload["release_files"]), "payload destination collision")
    for name in names:
        safe_member(name)
    return names


def verify_archive(directory, target):
    """Validate the exact whitelist before extracting a single member."""
    directory, target = Path(directory), Path(target)
    payload = read_json(directory / "payload.json")
    manifest = read_json(directory / "submission-manifest.json")
    archive = directory / "source.zip"
    require(digest(archive.read_bytes()) == manifest["source_archive_sha256"], "source ZIP hash mismatch")
    require(digest((directory / "announcement.pdf").read_bytes()) == manifest["reference_pdf_sha256"],
            "reference PDF hash mismatch")
    require(digest((directory / "metadata.template.json").read_bytes()) == manifest["metadata_sha256"],
            "metadata hash mismatch")
    expected = expected_members(payload)
    records = manifest["source_archive_members"]
    require(len(records) == len(expected) and {r["path"] for r in records} == expected,
            "submission manifest differs from the explicit payload whitelist")
    records = {r["path"]: r for r in records}
    with zipfile.ZipFile(archive) as saved:
        infos = saved.infolist()
        require(len(infos) == len(expected) and {x.filename for x in infos} == expected,
                "ZIP members differ from the explicit payload whitelist")
        require(sum(x.file_size for x in infos) <= MAX_TOTAL, "archive expands beyond the size limit")
        for info in infos:
            safe_member(info.filename)
            mode = stat.S_IFMT(info.external_attr >> 16)
            require(mode in {0, stat.S_IFREG} and not info.is_dir(), "nonregular ZIP member")
            require(info.file_size <= MAX_MEMBER, "archive member exceeds the size limit")
            record = records[info.filename]
            require(info.file_size == record["bytes"], "archive member size mismatch")
            data = saved.read(info)
            require(digest(data) == record["sha256"], "archive member hash mismatch")
        require(saved.read("anc/payload.json") == (directory / "payload.json").read_bytes(),
                "archived payload whitelist differs")
        require(saved.read("anc/verify_release.py") == (directory / "verify.py").read_bytes(),
                "archived verifier differs")
        for info in infos:
            destination = target.joinpath(*PurePosixPath(info.filename).parts)
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(saved.read(info))
    ancillary = read_json(target / "anc/manifest.json")
    for name, expected_hash in ancillary.get("included_files", {}).items():
        require(name in expected and digest((target / name).read_bytes()) == expected_hash,
                "ancillary manifest hash mismatch")
    source = (target / "main.tex").read_text()
    metadata = read_json(directory / "metadata.template.json")
    title = re.search(r"\\title\{([^{}]*)\}", source)
    abstract = source.split(r"\begin{abstract}", 1)[1].split(r"\end{abstract}", 1)[0]
    abstract = re.sub(r"(?<!\\)%[^\n]*", "", abstract).replace("$", "")
    abstract = re.sub(r"\s+", " ", abstract).strip()
    require(title is not None and title.group(1) == metadata["title"], "metadata title mismatch")
    require(abstract == metadata["abstract"], "metadata abstract mismatch")
    pointer = read_json(target / "anc/proof-sources.json")
    lean = read_json(target / "anc/lean-verification.json")
    require(lean["success"] is True and lean["proof_theorems_checked"] == 143
            and lean["additional_theorems_checked"] == 48
            and lean["project_axioms_or_placeholders_used"] is False,
            "current complete Lean verification report is missing")
    require(pointer["source_package_sha256"] == lean["archive_sha256"], "Lean packet/report hash mismatch")
    require(re.fullmatch(r"[0-9a-f]{40}", pointer["proof_commit"]) is not None
            and pointer["proof_commit"] in pointer["source_package_url"]
            and pointer["source_package_url"].startswith("https://github.com/sashakolpakov/butson_towers/"),
            "Lean source URL is not pinned to an immutable public revision")
    return payload, manifest


def run_recipe(name, root):
    """Small additional checks, also executable directly inside anc/."""
    root = Path(root)
    if name == "completed-inventory":
        ledger = read_json(root / "certificates/scalar-lift/completed-order-inventory.json")
        records = ledger["orders"]
        require(len(records) == len({x["order"] for x in records}) == 65,
                "completed order inventory is not a distinct65-order union")
        new = [x["order"] for x in records if x["status"] == "new_to_best_of_authors_knowledge"]
        known = [x["order"] for x in records if x["status"] == "known"]
        require(len(new) == 34 and len(known) == 31 and new == ledger["completed_new_real_orders"]
                and known == ledger["completed_known_real_orders"], "completed order labels differ")
        for relative, expected in ledger["source_certificate_sha256"].items():
            require(digest((root / relative).read_bytes()) == expected, "completed inventory source hash differs")
        print(json.dumps({"new_orders": len(new), "known_orders": len(known), "distinct_orders": len(records)}))
    elif name == "base-sixth-root-sources":
        import numpy as np
        from experiments.scalar_lift.reuse_hadamard_7572_frame import verify_source_gram
        saved = read_json(root / "certificates/scalar-lift/hadamard-7572-frame-reuse.json")
        orders = []
        for item in saved["outputs"]:
            exponents = np.array(item["exponents"], dtype=np.int64)
            proof = verify_source_gram(exponents)
            require(set(map(int, np.unique(exponents))) <= {0, 1, 3, 4}, "base source left the admitted palette")
            require(proof["order"] == item["source_order"]
                    and 2524*proof["order"] == item["output_order"], "base source order mismatch")
            orders.append(proof["order"])
        source3 = np.ones((3, 3), dtype=np.int64)
        np.fill_diagonal(source3, 3)
        verify_source_gram(source3)
        require(sorted(orders) == [7, 9, 13], "base source list differs")
        print(json.dumps({"all_complete_source_Grams_verified": True, "source_orders": [3, *orders]}))
    else:
        raise ValueError("unknown extracted recipe")


def replay(target, payload, *, scalar_algebra=False):
    root = Path(target) / "anc"
    (root / "_replay").mkdir()
    env = os.environ.copy()
    # -S below prevents editable-install .pth files from adding a checkout
    # outside the extraction. Expose installed third-party packages directly
    # without executing their site initialization or editable import hooks.
    library_paths = dict.fromkeys(sysconfig.get_paths()[name] for name in ("purelib", "platlib"))
    env["PYTHONPATH"] = os.pathsep.join((str(root), str(root / "src"), *library_paths))
    env.update(PYTHONNOUSERSITE="1", PYTHONDONTWRITEBYTECODE="1", OPENBLAS_NUM_THREADS="1",
               OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
    results = []
    for recipe in payload["replays"]:
        if recipe.get("optional_dependency") and not scalar_algebra:
            continue
        require(recipe["cwd"] == "anc", "unexpected replay working directory")
        started = time.monotonic()
        command = [sys.executable, "-S", *recipe["argv"]]
        completed = subprocess.run(command, cwd=root, env=env, text=True,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        require(completed.returncode == 0, f"{recipe['name']} replay failed:\n{completed.stdout}")
        comparison = recipe.get("compare_json")
        if comparison:
            require(read_json(root / comparison["output"]) == read_json(root / comparison["certificate"]),
                    f"{recipe['name']} regenerated certificate differs")
        results.append({"name": recipe["name"], "success": True,
                        "elapsed_seconds": round(time.monotonic()-started, 2)})
        print("Extracted " + recipe["name"] + " replay passed", flush=True)
    return results


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", type=Path, default=HERE)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--with-scalar-algebra", action="store_true")
    parser.add_argument("--recipe", choices=("completed-inventory", "base-sixth-root-sources"))
    args = parser.parse_args(argv)
    if args.recipe:
        run_recipe(args.recipe, Path.cwd())
        return
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="butson-hexagon-v2-replay-") as temp:
        payload, manifest = verify_archive(args.release_dir, Path(temp))
        results = [] if args.check_only else replay(Path(temp), payload, scalar_algebra=args.with_scalar_algebra)
    report = {"schema": "butson-hexagon-v2-extracted-replay-v1", "success": True,
              "source_archive_sha256": manifest["source_archive_sha256"],
              "all_explicit_archive_members_and_hashes_checked": True,
              "title_and_abstract_synchronized": True,
              "source_packet_pinned_to_complete_143_report_Lean_audit": True,
              "Lean_kernel_rerun_by_this_command": False,
              "rank_comparison_Lean_checked": False,
              "only_extracted_recipe_files_used": True,
              "inherited_PYTHONPATH_replaced_with_extracted_source_paths": True,
              "editable_install_pth_execution_disabled": True,
              "check_only": args.check_only, "replays": results,
              "optional_full_scalar_algebra_replay_requested": args.with_scalar_algebra,
              "remote_requests_performed": False,
              "elapsed_seconds": round(time.monotonic()-started, 2)}
    if args.report:
        args.report.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
