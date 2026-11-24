#!/usr/bin/env python3
"""Read-only privacy audit of a Git index, reachable history and releases.

Usage: python3 check_private_files.py --root /path/to/repository
Exit 0 means all checks passed; 1 means a finding or an incomplete check.
Only the two named root private files are read as privacy inputs. Credential
bytes stay in memory. Administrative preference leaf values are not scanned.
No credentials, credential hashes, Git object IDs or exception text are output.
"""
from __future__ import annotations

import argparse
import io
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import sys
import zipfile

PRIVATE_FILES = ("hexagon_key.local", "hexagon_parameters.local")
SKIP_DIRS = {".git", ".venv", ".lake"}
MAX_BYTES = 256 * 1024 * 1024
MAX_ARCHIVE_BYTES = 512 * 1024 * 1024
MAX_ARCHIVE_DEPTH = 4
OID = re.compile(rb"(?:[0-9a-f]{40}|[0-9a-f]{64})\Z")


class AuditError(Exception):
    """A deliberately content-free operational error."""


def is_private_path(name: str) -> bool:
    return any(part.casefold() == ".local" or part.casefold().endswith(".local")
               for part in name.replace("\\", "/").split("/"))


def credential_needles(raw: bytes) -> set[bytes]:
    """Support raw keys, KEY=value with optional quotes, and Bearer values."""
    raw = raw.strip()
    if not raw:
        return set()
    values = {raw}
    match = re.fullmatch(
        rb"(?:export\s+)?(?:HEXAGON_API_KEY|HEXAGON_TOKEN)\s*=\s*(.+)",
        raw, flags=re.DOTALL)
    value = match.group(1).strip() if match else raw
    if len(value) >= 2 and value[:1] == value[-1:] and value[:1] in (b"'", b'"'):
        value = value[1:-1].strip()
    values.add(value)
    if value[:7].lower() == b"bearer ":
        values.add(value[7:].strip())
    return {value for value in values if value}


class Audit:
    def __init__(self, root: Path):
        self.root = root
        self.needles: list[tuple[str, bytes]] = []
        self.findings: set[tuple[str, str]] = set()
        self.errors: set[tuple[str, str]] = set()
        self.counts = {
            "private_files_checked": 0, "tracked_paths_checked": 0,
            "git_index_blob_entries_checked": 0,
            "git_objects_checked": 0, "git_blobs_checked": 0,
            "git_index_only_blobs_checked": 0,
            "reachable_git_objects_checked": 0, "reachable_git_trees_checked": 0,
            "reachable_git_blobs_checked": 0, "git_blob_bytes_checked": 0,
            "release_archives_checked": 0, "zip_archives_checked": 0,
            "zip_members_checked": 0, "zip_expanded_bytes_checked": 0,
        }
        self.checks = {
            "root_private_files_present_regular_and_mode_0600": True,
            "credential_loaded": False, "parameters_loaded": False,
            "git_index_checked": False, "all_reachable_git_objects_checked": False,
            "all_git_index_blobs_checked": False,
            "release_zip_walk_completed": False,
        }

    def redact(self, path: str) -> str:
        raw = path.encode("utf-8", "surrogateescape")
        for _, needle in sorted(self.needles, key=lambda pair: -len(pair[1])):
            raw = raw.replace(needle, b"<redacted>")
        return raw.decode("utf-8", "replace")

    def fail(self, code: str, path: str, *, operational: bool = False):
        (self.errors if operational else self.findings).add((code, path))

    def load_private_inputs(self):
        for name in PRIVATE_FILES:
            path = self.root / name
            try:
                mode = path.lstat().st_mode
                if not stat.S_ISREG(mode) or stat.S_IMODE(mode) != 0o600:
                    self.checks["root_private_files_present_regular_and_mode_0600"] = False
                    self.fail("private_file_not_regular_mode_0600", name)
                if not stat.S_ISREG(mode):
                    continue
                if path.stat().st_size > MAX_BYTES:
                    raise AuditError()
                raw = path.read_bytes()
                self.counts["private_files_checked"] += 1
                if name == "hexagon_key.local":
                    needles = credential_needles(raw)
                    if not needles:
                        raise AuditError()
                    self.needles.extend(("credential_value", value) for value in needles)
                    self.checks["credential_loaded"] = True
                else:
                    if not raw.strip():
                        raise AuditError()
                    # Check whole-file copies only: public administrative leaf
                    # values legitimately appear in article metadata.
                    self.needles.append(("private_parameters_file_copy", raw))
                    self.checks["parameters_loaded"] = True
            except (OSError, AuditError):
                self.checks["root_private_files_present_regular_and_mode_0600"] = False
                self.fail("private_input_unavailable", name, operational=True)

    def inspect_content(self, data: bytes, location: str):
        for code, needle in self.needles:
            if needle in data:
                self.fail(code, location)

    def inspect_zip(self, data: bytes, location: str, depth: int = 0):
        stream = io.BytesIO(data)
        if not zipfile.is_zipfile(stream):
            self.fail("invalid_zip", location, operational=True)
            return
        if depth > MAX_ARCHIVE_DEPTH:
            self.fail("archive_depth_limit", location, operational=True)
            return
        self.counts["zip_archives_checked"] += 1
        expanded = 0
        try:
            with zipfile.ZipFile(stream) as archive:
                for info in archive.infolist():
                    self.counts["zip_members_checked"] += 1
                    member = location + "!" + info.filename
                    if is_private_path(info.filename):
                        self.fail("private_zip_member_path", member)
                    self.inspect_content(info.filename.encode("utf-8"), member)
                    if info.is_dir():
                        continue
                    expanded += info.file_size
                    if info.file_size > MAX_BYTES or expanded > MAX_ARCHIVE_BYTES:
                        self.fail("archive_size_limit", member, operational=True)
                        continue
                    # Bounded read; encrypted/corrupt/unsupported archives fail closed.
                    with archive.open(info) as entry:
                        content = entry.read(MAX_BYTES + 1)
                    if len(content) > MAX_BYTES:
                        self.fail("archive_size_limit", member, operational=True)
                        continue
                    self.counts["zip_expanded_bytes_checked"] += len(content)
                    self.inspect_content(content, member)
                    if zipfile.is_zipfile(io.BytesIO(content)):
                        self.inspect_zip(content, member, depth + 1)
        except (OSError, ValueError, RuntimeError, NotImplementedError,
                zipfile.BadZipFile, EOFError):
            self.fail("zip_read_failed", location, operational=True)

    def git(self, *args: str) -> bytes:
        try:
            result = subprocess.run(
                ["git", "--no-optional-locks", "-C", str(self.root), *args],
                stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=False)
        except OSError:
            raise AuditError() from None
        if result.returncode:
            raise AuditError()
        return result.stdout

    def inspect_git(self):
        try:
            indexed_blobs = []
            for entry in self.git("ls-files", "--stage", "-z").split(b"\0"):
                if not entry:
                    continue
                metadata, separator, raw = entry.partition(b"\t")
                fields = metadata.split()
                if not separator or len(fields) != 3 or not OID.fullmatch(fields[1]):
                    raise AuditError()
                mode, oid, _stage = fields
                path = raw.decode("utf-8", "surrogateescape")
                self.counts["tracked_paths_checked"] += 1
                if is_private_path(path):
                    self.fail("private_tracked_path", path)
                self.inspect_content(raw, path)
                # Gitlinks name commits in a separate repository, not file
                # contents. Regular files and symlinks have index blob IDs.
                if mode != b"160000":
                    indexed_blobs.append((oid, path))
                    self.counts["git_index_blob_entries_checked"] += 1
            self.checks["git_index_checked"] = True
            # rev-list supplies an example public path for each unique object.
            # Every tree is also inspected below, so aliases cannot hide a
            # private filename when its blob was already listed at a public path.
            objects = []
            for line in self.git("rev-list", "--objects", "--all").splitlines():
                oid, _, raw_path = line.partition(b" ")
                if not OID.fullmatch(oid):
                    raise AuditError()
                objects.append((oid, raw_path.decode("utf-8", "surrogateescape")))
            # Newly staged public files can contain private bytes before any
            # commit makes them reachable from a ref. Include all index blobs,
            # retaining every reachable tree for historical alias checks.
            reachable = {oid for oid, _path in objects}
            seen = set(reachable)
            for oid, path in indexed_blobs:
                if oid not in seen:
                    objects.append((oid, path))
                    seen.add(oid)
            process = subprocess.Popen(
                ["git", "--no-optional-locks", "-C", str(self.root), "cat-file", "--batch"],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
            try:
                assert process.stdin is not None and process.stdout is not None
                for oid, path in objects:
                    process.stdin.write(oid + b"\n")
                    process.stdin.flush()
                    header = process.stdout.readline().split()
                    if len(header) != 3 or header[0] != oid:
                        raise AuditError()
                    kind, size = header[1], int(header[2])
                    if size > MAX_BYTES:
                        raise AuditError()
                    data = process.stdout.read(size)
                    if len(data) != size or process.stdout.read(1) != b"\n":
                        raise AuditError()
                    self.counts["git_objects_checked"] += 1
                    self.counts["reachable_git_objects_checked"] += int(oid in reachable)
                    if kind == b"tree":
                        self.counts["reachable_git_trees_checked"] += 1
                        offset, oid_bytes = 0, len(oid) // 2
                        while offset < len(data):
                            zero = data.index(b"\0", offset)
                            raw_name = data[offset:zero].split(b" ", 1)[1]
                            name = raw_name.decode("utf-8", "surrogateescape")
                            location = (path + "/" if path else "") + name
                            if is_private_path(name):
                                self.fail("private_historical_path", location)
                            self.inspect_content(raw_name, location)
                            offset = zero + 1 + oid_bytes
                        if offset != len(data):
                            raise AuditError()
                    elif kind == b"blob":
                        self.counts["git_blobs_checked"] += 1
                        self.counts["reachable_git_blobs_checked"] += int(oid in reachable)
                        self.counts["git_index_only_blobs_checked"] += int(oid not in reachable)
                        self.counts["git_blob_bytes_checked"] += len(data)
                        location = path or "<reachable-git-blob>"
                        self.inspect_content(data, location)
                        if zipfile.is_zipfile(io.BytesIO(data)):
                            self.inspect_zip(data, location)
                process.stdin.close()
                if process.wait() != 0:
                    raise AuditError()
                self.checks["all_reachable_git_objects_checked"] = True
                self.checks["all_git_index_blobs_checked"] = True
            finally:
                if process.poll() is None:
                    process.kill()
                    process.wait()
                if process.stdin and not process.stdin.closed:
                    process.stdin.close()
                if process.stdout:
                    process.stdout.close()
        except (OSError, ValueError, IndexError, AuditError):
            self.fail("git_history_check_incomplete", "<git-history>", operational=True)

    def inspect_releases(self):
        release = self.root / "release"
        try:
            def walk_error(_error):
                raise AuditError()
            if release.is_symlink():
                raise AuditError()
            for directory, dirs, files in os.walk(release, onerror=walk_error, followlinks=False):
                retained = []
                for name in dirs:
                    if name in SKIP_DIRS or is_private_path(name):
                        continue
                    path = Path(directory) / name
                    if path.is_symlink():
                        self.fail("release_directory_symlink_not_scanned",
                                  path.relative_to(self.root).as_posix(), operational=True)
                    else:
                        retained.append(name)
                dirs[:] = retained
                for name in files:
                    if not name.casefold().endswith(".zip"):
                        continue
                    path = Path(directory) / name
                    relative = path.relative_to(self.root).as_posix()
                    self.counts["release_archives_checked"] += 1
                    if is_private_path(relative):
                        self.fail("private_release_archive_path", relative)
                    if path.is_symlink() or path.stat().st_size > MAX_BYTES:
                        self.fail("release_archive_unreadable", relative, operational=True)
                        continue
                    self.inspect_zip(path.read_bytes(), relative)
            self.checks["release_zip_walk_completed"] = True
        except (OSError, ValueError, AuditError):
            self.fail("release_walk_incomplete", "release", operational=True)

    def run(self) -> dict:
        self.load_private_inputs()
        self.inspect_git()
        self.inspect_releases()
        return {
            "ok": not self.findings and not self.errors and all(self.checks.values()),
            "checks": self.checks,
            "counts": {**self.counts, "findings": len(self.findings), "errors": len(self.errors)},
            "findings": [{"check": code, "path": self.redact(path)}
                         for code, path in sorted(self.findings)],
            "errors": [{"check": code, "path": self.redact(path)}
                       for code, path in sorted(self.errors)],
        }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        report = Audit(args.root.resolve()).run()
    except Exception:
        # Exception messages can contain file contents or credential-bearing
        # names. Deliberately emit neither exception text nor a traceback.
        report = {"ok": False, "counts": {"errors": 1},
                  "errors": [{"check": "audit_incomplete", "path": "<audit>"}]}
    print(json.dumps(report, sort_keys=True))
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
