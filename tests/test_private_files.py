"""Privacy checker controls use fake credentials and temporary Git repositories."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
import zipfile

CHECKER_PATH = Path(__file__).with_name("check_private_files.py")
if not CHECKER_PATH.exists():
    CHECKER_PATH = Path(__file__).resolve().parents[1] / "scripts" / "check_private_files.py"
SPEC = importlib.util.spec_from_file_location("privacy_check", CHECKER_PATH)
CHECKER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECKER)
FAKE_KEY = "synthetic-test-credential-for-privacy-check-0123456789"


class PrivacyChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.git("init", "-q")
        self.git("config", "user.email", "fake@example.invalid")
        self.git("config", "user.name", "Synthetic privacy test")
        (self.root / ".gitignore").write_text("*.local\n.local/\n")
        (self.root / "README.md").write_text("Public references to hexagon_key.local are permitted.\n")
        self.key('HEXAGON_API_KEY="Bearer ' + FAKE_KEY + '"\n')
        preferences = self.root / "hexagon_parameters.local"
        preferences.write_text('{"publicAdministrativeChoice":"permitted-in-metadata","privateFileMarker":"synthetic-root-preferences"}\n')
        preferences.chmod(0o600)
        self.commit()

    def tearDown(self):
        self.temp.cleanup()

    def git(self, *args):
        return subprocess.run(["git", "-C", str(self.root), *args], check=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout

    def commit(self):
        self.git("add", "-A")
        self.git("commit", "-qm", "Synthetic fixture")

    def key(self, value):
        path = self.root / "hexagon_key.local"
        path.write_text(value)
        path.chmod(0o600)

    def archive(self, members):
        path = self.root / "release" / "fixture" / "bundle.zip"
        path.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as output:
            for name, data in members.items():
                output.writestr(name, data)
        return path

    def audit(self):
        result = CHECKER.Audit(self.root).run()
        self.assertNotIn(FAKE_KEY, json.dumps(result))
        return result

    def test_safe_repository_and_zip(self):
        self.archive({"metadata.json": '{"choice":"permitted-in-metadata"}',
                      "README.md": "Configuration filename hexagon_parameters.local is documented."})
        result = self.audit()
        self.assertTrue(result["ok"], result)
        self.assertEqual(result["counts"]["release_archives_checked"], 1)

    def test_reject_private_zip_member(self):
        self.archive({"anc/.local/receipt.json": "fake", "nested/private.local/file": "fake"})
        result = self.audit()
        self.assertFalse(result["ok"])
        self.assertEqual(sum(x["check"] == "private_zip_member_path" for x in result["findings"]), 2)

    def test_reject_fake_key_in_zip(self):
        self.archive({"innocent.txt": "prefix " + FAKE_KEY + " suffix"})
        self.assertFalse(self.audit()["ok"])

    def test_reject_staged_public_secret_before_commit(self):
        path = self.root / "article-metadata.json"
        path.write_text(json.dumps({"accidentalCredential": FAKE_KEY}))
        self.git("add", "article-metadata.json")
        result = self.audit()
        self.assertFalse(result["ok"])
        self.assertTrue(any(x["check"] == "credential_value" and x["path"] == "article-metadata.json"
                            for x in result["findings"]))

    def test_reject_staged_compressed_zip_outside_release(self):
        path = self.root / "paper-source.zip"
        with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as output:
            output.writestr("innocent.txt", FAKE_KEY)
        self.git("add", "paper-source.zip")
        result = self.audit()
        self.assertFalse(result["ok"])
        self.assertTrue(any(x["check"] == "credential_value" and x["path"] == "paper-source.zip!innocent.txt"
                            for x in result["findings"]))

    def test_reject_deleted_historical_secret(self):
        path = self.root / "public.txt"
        path.write_text(FAKE_KEY)
        self.commit()
        path.unlink()
        self.commit()
        result = self.audit()
        self.assertFalse(result["ok"])
        self.assertTrue(any(x["check"] == "credential_value" for x in result["findings"]))

    def test_private_historical_alias_of_public_blob(self):
        path = self.root / "alias.local"
        path.write_bytes((self.root / "README.md").read_bytes())
        self.git("add", "-f", "alias.local")
        self.git("commit", "-qm", "Synthetic alias")
        path.unlink()
        self.commit()
        result = self.audit()
        self.assertFalse(result["ok"])
        self.assertTrue(any(x["check"] == "private_historical_path" for x in result["findings"]))

    def test_deleted_historical_compressed_zip_secret(self):
        path = self.archive({"innocent.txt": FAKE_KEY})
        self.commit()
        path.unlink()
        self.commit()
        result = self.audit()
        self.assertFalse(result["ok"])
        self.assertTrue(any(x["check"] == "credential_value" and "!innocent.txt" in x["path"]
                            for x in result["findings"]))

    def test_release_symlink_fails_closed(self):
        external = self.root / "external"
        external.mkdir()
        release = self.root / "release"
        release.mkdir()
        (release / "linked").symlink_to(external, target_is_directory=True)
        result = self.audit()
        self.assertFalse(result["ok"])
        self.assertTrue(any(x["check"] == "release_directory_symlink_not_scanned"
                            for x in result["errors"]))

    def test_reject_tracked_private_directory(self):
        path = self.root / "anything.local" / "nested.txt"
        path.parent.mkdir()
        path.write_text("fake")
        self.git("add", "-f", "anything.local/nested.txt")
        result = self.audit()
        self.assertFalse(result["ok"])
        self.assertTrue(any(x["check"] == "private_tracked_path" for x in result["findings"]))

    def test_key_formats_and_mode(self):
        for key in [FAKE_KEY, "Bearer " + FAKE_KEY,
                    "HEXAGON_TOKEN='" + FAKE_KEY + "'",
                    "HEXAGON_API_KEY=Bearer " + FAKE_KEY]:
            with self.subTest(format=key.split("=", 1)[0] if "=" in key else "raw-or-bearer"):
                self.key(key)
                self.archive({"member": FAKE_KEY})
                self.assertFalse(self.audit()["ok"])
        (self.root / "hexagon_key.local").chmod(0o644)
        self.assertFalse(self.audit()["ok"])

    def test_reject_parameters_whole_copy_but_not_leaf(self):
        self.archive({"metadata.json": (self.root / "hexagon_parameters.local").read_bytes()})
        result = self.audit()
        self.assertFalse(result["ok"])
        self.assertTrue(any(x["check"] == "private_parameters_file_copy" for x in result["findings"]))

    def test_cli_redacts_credential_bearing_public_path(self):
        self.archive({FAKE_KEY + ".txt": FAKE_KEY})
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            code = CHECKER.main(["--root", str(self.root)])
        self.assertEqual(code, 1)
        self.assertNotIn(FAKE_KEY, output.getvalue())
        self.assertIn("redacted", output.getvalue())


if __name__ == "__main__":
    unittest.main()
