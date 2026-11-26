#!/usr/bin/env python3
"""Bind the additional Lean literals to the public source records.

This is a provenance check. The Lean kernel independently proves the
sign, symmetry, complete correlation, strong-starter and Paley identities.
Only the Python standard library is required.
"""
from pathlib import Path
import json
import re

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "formal/ButsonTowers/Data"
WIDTHS = (47, 53, 59, 67, 69, 73, 79, 93, 97, 103, 109, 127)
STARTER_SIGNS = {
    5: (1, -1, 1, 1), 53: (1, -1, 1, 1), 97: (1, 1, 1, -1),
    109: (1, 1, 1, 1), 113: (1, 1, 1, 1), 137: (1, -1, 1, 1),
}


def tree(values):
    if len(values) == 1:
        assert type(values[0]) is int
        return f"(.leaf ({values[0]}))"
    split = len(values) // 2
    return f"(.node {split} {tree(values[:split])} {tree(values[split:])})"


def rows(record):
    words = record["first_rows_pm"]
    assert len(words) == 4 and all(set(word) <= {"+", "-"} for word in words)
    return [[1 if c == "+" else -1 for c in word] for word in words]


def check_quad(path, name, values):
    expected = (f"def {name} : Fin 4 → LiteralTree ℤ :=\n  !["
                + ",\n    ".join(tree(row) for row in values) + "]\n")
    matches = re.findall(rf"^def {name} : Fin 4 → LiteralTree ℤ :=\n  !\[.*?\]\n",
                         path.read_text(), re.MULTILINE | re.DOTALL)
    if matches != [expected]:
        raise ValueError(f"Lean literal differs from public rows: {path}:{name}")


def main():
    old = json.loads((ROOT / "experiments/scalar_lift/williamson_inputs.json").read_text())
    near = json.loads((ROOT / "experiments/scalar_lift/near_williamson_inputs.json").read_text())
    records = near.get("entries", near.get("records", []))
    by_width = {record.get("width", record.get("order")): record for record in records}
    old_by_width = {record["order"]: record for record in old["entries"]}
    for width in WIDTHS:
        record = by_width[width]
        values = rows(record)
        assert all(len(row) == width for row in values)
        path = DATA / f"NearWilliamson/W{width}.lean"
        check_quad(path, f"nearWilliamson{width}Trees", values)
        constructor = "nearMatrices" if record["first_matrix_reversed"] else "rowMatrices"
        definition = (f"def nearWilliamson{width} : Fin 4 → Matrix (ZMod {width}) "
                      f"(ZMod {width}) ℚ :=\n  {constructor} {width} nearWilliamson{width}Row\n")
        if definition not in path.read_text():
            raise ValueError(f"Lean reversal convention differs from the source record: {path}")
    for q, signs in STARTER_SIGNS.items():
        width = (q + 1) // 2
        record = by_width[width] if q == 137 else old_by_width[width]
        values = [[sign * value for value in row] for sign, row in zip(signs, rows(record))]
        path = DATA / f"AdditionalSpence/Q{q}.lean"
        check_quad(path, f"additionalSpence{q}Trees", values)
        paley = [0 if x == 0 else (1 if pow(x, (q - 1) // 2, q) == 1 else -1)
                 for x in range(q)]
        expected = f"def additionalPaley{q}Tree : LiteralTree ℤ :=\n  {tree(paley)}\n"
        if expected not in path.read_text():
            raise ValueError(f"Lean literal Paley row differs: {path}")
    print(json.dumps({"ok": True, "near_Williamson_widths": list(WIDTHS),
                      "additional_Spence_parameters": list(STARTER_SIGNS),
                      "literal_modules_checked": len(WIDTHS) + len(STARTER_SIGNS),
                      "data_domains": ["near-Williamson rows and reversal flags",
                                       "signed strong-starter rows", "Paley rows"]}))


if __name__ == "__main__":
    main()
