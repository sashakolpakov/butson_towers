import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest


@pytest.fixture(scope="module")
def binary(tmp_path_factory):
    compiler = shutil.which("clang++")
    if compiler is None:
        pytest.skip("the optional congruence helper requires C++17")
    target = tmp_path_factory.mktemp("congruence")/"helper"
    source = Path(__file__).parents[1]/"experiments/scalar_lift/sextic_row_congruence.cpp"
    subprocess.run([compiler, "-O3", "-std=c++17", str(source), "-o", str(target)], check=True)
    return target


def invoke(binary, path, bits, target, selected):
    n, m = bits.shape
    path.write_text(f"{n} {m} {len(selected)}\n"+" ".join(map(str, target))+"\n"
                    +"\n".join(" ".join(map(str, row)) for row in bits)+"\n"
                    +" ".join(map(str, selected))+"\n")
    result = json.loads(subprocess.check_output([str(binary), str(path), ".02", "717"], text=True))
    actual = np.asarray(result["selected"])
    residual = bits[actual].sum(axis=0)-target
    assert len(set(map(int, actual))) == len(selected)
    assert np.array_equal(residual, result["residual"])
    assert int(residual @ residual) == result["best_energy"]
    return result


def test_exact_syndrome_matching_finds_a_complete_binary_witness(binary, tmp_path):
    bits = np.eye(10, dtype=np.int64)
    target = np.asarray((1, 1, 1, 1, 0, 0, 0, 0, 0, 0))
    result = invoke(binary, tmp_path/"exact.case", bits, target, (4, 5, 6, 3))
    assert not result["initial_mod2_sector"]
    assert result["best_energy"] == 0
    assert result["all_mod2_passed"]
    assert result["selected"] == [0, 1, 2, 3]


def test_sector_claim_is_not_made_for_an_impossible_syndrome(binary, tmp_path):
    result = invoke(binary, tmp_path/"impossible.case", np.zeros((8, 2), dtype=np.int64),
                    np.asarray((1, 0)), (0, 1, 2))
    assert not result["all_mod2_passed"]
    assert result["best_energy"] == 1


def test_neutral_moves_preserve_the_passed_sector_without_claiming_a_zero(binary, tmp_path):
    result = invoke(binary, tmp_path/"even.case", np.zeros((8, 2), dtype=np.int64),
                    np.asarray((2, 0)), (0, 1, 2))
    assert result["all_mod2_passed"]
    assert result["passed_modulus"] == 2
    assert result["best_energy"] == 4
    assert result["syndrome_matches"] > 1
