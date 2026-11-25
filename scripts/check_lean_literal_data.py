#!/usr/bin/env python3
"""Regenerate/check public Lean literal data using only the standard library.

The JSON inputs and these comparisons are provenance, not Lean proof premises.
Lean checks the emitted literal source properties directly in its kernel.
Default: check. --write updates literal declarations without changing proofs.
"""
import argparse
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "formal/ButsonTowers/Data"


def load(path):
    return json.loads((ROOT / path).read_text())


def tree(values):
    if len(values) == 1:
        assert type(values[0]) is int
        return f"(.leaf ({values[0]}))"
    split = len(values) // 2
    return f"(.node {split} {tree(values[:split])} {tree(values[split:])})"


def signs(words):
    assert all(set(word) <= {"+", "-"} for word in words)
    return [[1 if symbol == "+" else -1 for symbol in word] for word in words]


def expected_trees():
    parts = ["import ButsonTowers.Data.LiteralTree\n\nnamespace ButsonTowers.Data\n"]
    for p, path in [(127, "artifacts/t-matrices-127/t127.json"),
                    (631, "artifacts/hadamard-7572/t631.json")]:
        rows = load(path)["b_sequences"]
        assert len(rows) == 4 and all(len(row) == p for row in rows)
        codes = []
        for column in zip(*rows):
            assert all(type(value) is int and value in (-1, 0, 1) for value in column)
            assert sum(abs(value) for value in column) == 1
            channel = next(i for i, value in enumerate(column) if value)
            codes.append(2 * channel + int(column[channel] < 0))
        parts.append(f"def t{p}Tree : LiteralTree ℕ :=\n  {tree(codes)}\n"
                     f"def t{p}CodesFast : ℕ → ℕ := treeAt t{p}Tree\n")
    for item in load("experiments/scalar_lift/williamson_inputs.json")["entries"]:
        rows = signs(item["first_rows_pm"])
        parts.append(f'def williamson{item["order"]}Trees : Fin 4 → LiteralTree ℤ :=\n  !['
                     + ",\n    ".join(tree(row) for row in rows) + "]\n")
    family = load("certificates/scalar-lift/williamson-spence-family.json")
    for item in family["starters"]:
        rows = signs(item["literal_first_rows_pm"])
        parts.append(f'def spence{item["q"]}Trees : Fin 4 → LiteralTree ℤ :=\n  !['
                     + ",\n    ".join(tree(row) for row in rows) + "]\n")
    for item in family["Paley_parameters"]:
        parts.append(f'def paley{item["q"]}Tree : LiteralTree ℤ :=\n  '
                     + tree(item["literal_quadratic_character_first_row"]) + "\n")
    parts.append("end ButsonTowers.Data\n")
    return "\n".join(parts)


def expected_unreal_trees():
    seed = load("experiments/scalar_lift/bh13_6_seed.json")
    e = seed["exponents"]
    r, c = seed["unreal_row_phase_shifts"], seed["unreal_column_phase_shifts"]
    assert len(e) == len(r) == len(c) == 13 and all(len(row) == 13 for row in e)
    assert all(type(value) is int and 0 <= value < 6 for row in e for value in row)
    assert all(type(value) is int for value in r + c)
    phases = [(e[i][j] + r[i] + c[j]) % 6 for i in range(13) for j in range(13)]
    assert set(phases) <= {1, 2, 4, 5}
    a = [{1: 1, 2: -1, 4: -1, 5: 1}[phase] for phase in phases]
    b = [{1: 1, 2: 1, 4: -1, 5: -1}[phase] for phase in phases]
    return {"unreal13ATree": tree(a), "unreal13BTree": tree(b)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    path = DATA / "TreeLiterals.lean"
    expected = expected_trees()
    changed = 0
    if path.read_text() != expected:
        if not args.write:
            raise SystemExit("Lean tree literal data differ from the public source records")
        path.write_text(expected)
        changed += 1
    path = DATA / "Unreal13.lean"
    content = path.read_text()
    for name, value in expected_unreal_trees().items():
        line = f"def {name} : LiteralTree ℤ := {value}"
        pattern = rf"^def {name} : LiteralTree ℤ := .*?$"
        matches = re.findall(pattern, content, flags=re.MULTILINE)
        assert len(matches) == 1
        if matches[0] != line:
            if not args.write:
                raise SystemExit("Lean unreal source literals differ from the public source records")
            content = re.sub(pattern, lambda _: line, content, flags=re.MULTILINE)
    if content != path.read_text():
        path.write_text(content)
        changed += 1
    print(json.dumps({"ok": True, "literal_modules_checked": 2, "modules_updated": changed}))


if __name__ == "__main__":
    main()
