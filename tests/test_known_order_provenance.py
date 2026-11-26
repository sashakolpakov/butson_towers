"""Replay and corruption controls for the known-order provenance certificate."""

import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/audit_known_order_provenance.py"
SPEC = importlib.util.spec_from_file_location("known_order_provenance", SCRIPT)
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


@pytest.fixture
def copied_inputs(tmp_path):
    report = audit.build_report(ROOT)
    for relative in [*report["source_sha256"], audit.CERTIFICATE]:
        target = tmp_path / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / relative, target)
    return tmp_path


def rewrite(path, change):
    value = json.loads(path.read_text())
    change(value)
    path.write_text(json.dumps(value))


def test_public_certificate_replays_exactly():
    report = audit.build_report(ROOT)
    assert (ROOT / audit.CERTIFICATE).read_bytes() == audit.canonical_bytes(report)
    assert report["known_order_count"] == 31
    assert report["endpoint_selectors_checked"] == 32
    assert report["ledger_status_selectors_checked"] == 32
    assert report["detailed_prior_proofs_checked"] == 32
    assert sum(report["construction_counts"].values()) == 31


@pytest.mark.parametrize("raw,order", [
    ({"family": "Paley I", "parameters": {"q": 1521}}, 1522),
    ({"family": "Paley I", "parameters": {"q": 13}}, 14),
    ({"family": "Paley II", "parameters": {"q": 13}}, 30),
    ({"family": "Paley sixteen-block product", "q": 19, "s": 5, "t": 21}, 420),
])
def test_false_positive_construction_hypotheses_are_rejected(raw, order):
    with pytest.raises(audit.AuditError, match="failed integer hypotheses"):
        audit.normalize_prior(order, raw)


def test_strict_selector_rejects_skipped_syntax_and_missing_keys():
    assert audit.select({"items": [{"order": 1524}]}, "items[0].order") == 1524
    for selector in ("items[0]..order", "items[0]junk", "items[-1]", "items[9]", ""):
        with pytest.raises(audit.AuditError):
            audit.select({"items": [{"order": 1524}]}, selector)


def test_changed_endpoint_order_is_rejected(copied_inputs):
    path = copied_inputs / "certificates/scalar-lift/hadamard-frame-p127.json"
    rewrite(path, lambda value: value["known_order_positive_repeat"]["endpoint"].update(order=1528))
    with pytest.raises(audit.AuditError, match="endpoint selector has a different order"):
        audit.build_report(copied_inputs)


def test_changed_ledger_status_source_is_rejected(copied_inputs):
    path = copied_inputs / "certificates/scalar-lift/higher-order-records.json"
    rewrite(path, lambda value: value["completed_order_statuses"][0].update(status="new"))
    with pytest.raises(audit.AuditError, match="ledger status selector disagrees"):
        audit.build_report(copied_inputs)


def test_check_rejects_input_hash_drift_without_a_semantic_change(copied_inputs, capsys):
    path = copied_inputs / "experiments/scalar_lift/williamson_order_status.json"
    path.write_bytes(path.read_bytes() + b"\n")
    assert audit.main(["--repo", str(copied_inputs), "--check"]) == 1
    assert "differs from current inputs" in capsys.readouterr().err
    assert audit.main(["--repo", str(copied_inputs)]) == 0
    assert audit.main(["--repo", str(copied_inputs), "--check"]) == 0


def test_check_rejects_changed_certificate(copied_inputs, capsys):
    rewrite(copied_inputs / audit.CERTIFICATE,
            lambda value: value.update(known_order_count=30))
    assert audit.main(["--repo", str(copied_inputs), "--check"]) == 1
    assert "differs from current inputs" in capsys.readouterr().err


def test_invalid_hypotheses_remain_rejected_under_python_optimization(copied_inputs):
    path = copied_inputs / "experiments/scalar_lift/williamson_order_status.json"
    rewrite(path, lambda value: value["631:5"]["prior_construction"]["parameters"].update(q=12615))
    result = subprocess.run([sys.executable, "-O", str(SCRIPT), "--repo", str(copied_inputs)],
                            capture_output=True, text=True, check=False)
    assert result.returncode == 1
    assert "failed integer hypotheses" in result.stderr
