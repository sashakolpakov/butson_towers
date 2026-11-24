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
        pytest.skip("the optional macro helper requires C++17")
    target = tmp_path_factory.mktemp("mod4")/"helper"
    source = Path(__file__).parents[1]/"experiments/scalar_lift/sextic_row_mod4.cpp"
    subprocess.run([compiler, "-O3", "-std=c++17", str(source), "-o", str(target)], check=True)
    return target


def invoke(binary, case, bits, target, selected):
    n, m = bits.shape
    case.write_text(f"{n} {m} {len(selected)}\n"+" ".join(map(str, target))+"\n"
                    +"\n".join(" ".join(map(str, row)) for row in bits)+"\n"
                    +" ".join(map(str, selected))+"\n")
    result = json.loads(subprocess.check_output([str(binary), str(case), ".05", "171"], text=True))
    actual = np.asarray(result["selected"])
    residual = bits[actual].sum(axis=0)-target
    assert len(set(map(int, actual))) == len(selected)
    assert np.array_equal(residual, result["residual"])
    assert int(residual @ residual) == result["best_energy"]
    assert np.all(residual % result["passed_modulus"] == 0)
    assert Path(str(case)+".kernel.tsv").exists()
    return result


def test_three_disjoint_kernel_moves_compose_to_an_exact_planted_lift(binary, tmp_path):
    rng = np.random.default_rng(2026112421)
    bits = rng.integers(0, 2, size=(24, 16), dtype=np.int64)
    for group in range(3):
        out = list(range(3*group, 3*group+3))
        incoming = list(range(9+3*group, 12+3*group))
        bits[incoming[2]] = np.bitwise_xor.reduce(bits[out+incoming[:2]], axis=0)
    selected = np.arange(9)
    target = bits[np.arange(9, 18)].sum(axis=0)
    result = invoke(binary, tmp_path/"planted.case", bits, target, selected)
    assert result["best_energy"] == 0
    assert result["best_macro_moves"] == 3
    assert len(result["chosen_macro"]) == 3
    removed, inserted = [], []
    for move in result["chosen_macro"]:
        removed.extend(move["out"])
        inserted.extend(move["in"])
        delta = bits[move["in"]].sum(axis=0)-bits[move["out"]].sum(axis=0)
        assert np.all(delta % 2 == 0)
        assert sum(1 << q for q, x in enumerate(delta) if x > 0) == move["positive_half_mask"]
        assert sum(1 << q for q, x in enumerate(delta) if x < 0) == move["negative_half_mask"]
    assert len(set(removed)) == len(set(inserted)) == 9


def test_missing_catalogue_span_does_not_claim_a_mod4_lift(binary, tmp_path):
    result = invoke(binary, tmp_path/"no-span.case", np.zeros((8, 2), dtype=np.int64),
                    np.asarray((2, 0)), np.arange(3))
    assert result["half_syndrome_rank"] == 0
    assert not result["target_in_catalogue_syndrome_span"]
    assert not result["mod4_lift_found"]
    assert result["passed_modulus"] == 2


def test_continuing_a_mod4_state_cannot_drop_to_mod2(binary, tmp_path):
    result = invoke(binary, tmp_path/"preserve.case", np.zeros((8, 2), dtype=np.int64),
                    np.asarray((4, 0)), np.arange(3))
    assert result["best_energy"] == 16
    assert result["passed_modulus"] == 4
    assert result["mod4_lift_found"]
