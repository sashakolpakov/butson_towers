import hashlib
import importlib.util
import json
from pathlib import Path
import zipfile
from types import SimpleNamespace

import pytest


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "release/hexagon_2611_01070_v2/verify.py"
spec = importlib.util.spec_from_file_location("hexagon_v2_verifier", SOURCE)
verifier = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verifier)


def write_release(directory):
    directory.mkdir()
    payload = {"schema": "butson-hexagon-v2-payload-v1", "repository_files": [],
               "release_files": {"proof-sources.json": "anc/proof-sources.json",
                                 "lean-verification.json": "anc/lean-verification.json",
                                 "payload.json": "anc/payload.json", "verify.py": "anc/verify_release.py"}}
    commit, packet = "a"*40, "b"*64
    pointer = {"proof_commit": commit, "source_package_sha256": packet,
               "source_package_url": f"https://github.com/sashakolpakov/butson_towers/raw/{commit}/packet.zip"}
    lean = {"success": True, "proof_theorems_checked": 143, "additional_theorems_checked": 48,
            "project_axioms_or_placeholders_used": False, "archive_sha256": packet}
    def encoded(value):
        return (json.dumps(value, indent=2)+"\n").encode()
    files = {"main.tex": b"\\title{Test}\n\\begin{abstract}Exact $1$ signs.\\end{abstract}\n",
             "anc/payload.json": encoded(payload), "anc/manifest.json": encoded({"included_files": {}}),
             "anc/proof-sources.json": encoded(pointer), "anc/lean-verification.json": encoded(lean),
             "anc/verify_release.py": SOURCE.read_bytes()}
    (directory / "payload.json").write_bytes(files["anc/payload.json"])
    (directory / "verify.py").write_bytes(SOURCE.read_bytes())
    (directory / "announcement.pdf").write_bytes(b"test reference PDF bytes")
    (directory / "metadata.template.json").write_bytes(encoded({"title": "Test", "abstract": "Exact 1 signs."}))
    archive = directory / "source.zip"
    with zipfile.ZipFile(archive, "w") as out:
        for name, value in files.items():
            out.writestr(name, value)
    manifest = {"source_archive_sha256": verifier.digest(archive.read_bytes()),
                "reference_pdf_sha256": verifier.digest((directory / "announcement.pdf").read_bytes()),
                "metadata_sha256": verifier.digest((directory / "metadata.template.json").read_bytes()),
                "source_archive_members": [{"path": p, "bytes": len(v), "sha256": verifier.digest(v)}
                                           for p, v in files.items()]}
    (directory / "submission-manifest.json").write_bytes(encoded(manifest))
    return manifest


def test_minimal_archive_checks_all_hashes_and_proof_pointer(tmp_path):
    directory = tmp_path / "release"
    manifest = write_release(directory)
    destination = tmp_path / "extracted"
    payload, saved = verifier.verify_archive(directory, destination)
    assert saved == manifest
    assert {str(p.relative_to(destination)) for p in destination.rglob("*") if p.is_file()} == (
        verifier.expected_members(payload)
    )
    damaged = json.loads((directory / "submission-manifest.json").read_text())
    damaged["source_archive_members"][0]["sha256"] = "0"*64
    (directory / "submission-manifest.json").write_text(json.dumps(damaged))
    with pytest.raises(ValueError, match="member hash mismatch"):
        verifier.verify_archive(directory, tmp_path / "bad")
    assert not (tmp_path / "bad").exists()


@pytest.mark.parametrize("name", ["../outside.py", "anc/../outside.py", "anc\\escape.py",
                                 "anc/.local/key.json", "anc/hexagon_key.local", "anc/formal.lean",
                                 "anc/nested.zip", "anc/./same.py", "/anc/absolute.py"])
def test_private_unsafe_and_unaccepted_archive_paths_are_rejected(name):
    with pytest.raises(ValueError):
        verifier.safe_member(name)


def test_manifest_cannot_sneak_extra_file_into_whitelist(tmp_path):
    directory = tmp_path / "release"
    manifest = write_release(directory)
    manifest["source_archive_members"].append({"path": "anc/extra.py", "bytes": 0, "sha256": hashlib.sha256(b"").hexdigest()})
    (directory / "submission-manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="explicit payload whitelist"):
        verifier.verify_archive(directory, tmp_path / "bad")


def test_child_replay_cannot_use_editable_install_paths(tmp_path, monkeypatch):
    (tmp_path / "anc").mkdir()
    monkeypatch.setenv("PYTHONPATH", "/a/forbidden/editable/checkout")
    seen = []
    def fake_run(command, **kwargs):
        seen.append((command, kwargs))
        return SimpleNamespace(returncode=0, stdout="")
    monkeypatch.setattr(verifier.subprocess, "run", fake_run)
    payload = {"replays": [{"name": "test", "argv": ["a.py"], "cwd": "anc"}]}
    assert verifier.replay(tmp_path, payload)[0]["success"]
    command, kwargs = seen[0]
    assert command[1:] == ["-S", "a.py"]
    paths = kwargs["env"]["PYTHONPATH"].split(verifier.os.pathsep)
    assert paths[:2] == [str(tmp_path / "anc"), str(tmp_path / "anc/src")]
    assert "/a/forbidden/editable/checkout" not in paths
