#!/usr/bin/env python3
"""Typeset and prepare the versioned announcement locally, without network access."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CANONICAL = ROOT / "manuscript/hadamard_7572_announcement.tex"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def prose(text):
    """Normalize the announcement's small TeX prose/math subset for API text."""
    text = re.sub(r"(?<!\\)%[^\n]*", "", text)
    text = re.sub(r"\\(?:begin|end)\{(?:gathered|aligned|array)\}", " ", text)
    text = re.sub(r"\\(?:\(|\)|\[|\])", " ", text)
    text = re.sub(r"\\(?:,|!|;|quad\b|qquad\b|allowbreak\b)", " ", text)
    text = text.replace(r"\\", " ").replace("$", "").replace("&", " ")
    text = re.sub(r"\\(?:text|textrm|mathrm|mathbf|emph|textit|textbf)\{([^{}]*)\}", r"\1", text)
    text = re.sub(r"\\geq?(?![A-Za-z])", ">=", text)
    text = re.sub(r"\\(?:ldots|cdots|dots)\b", "...", text)
    text = text.replace(r"\times", "*").replace(r"\zeta", "zeta")
    text = text.replace("{", "").replace("}", "").replace("~", " ")
    if "\\" in text:
        raise ValueError("Unsupported TeX command in article metadata; extend prose() explicitly")
    text = re.sub(r"\s+([,.;:])", r"\1", text)
    return re.sub(r"\s+", " ", text).strip()


def title_argument(source):
    start = source.index(r"\title{") + len(r"\title{")
    depth = 1
    for index in range(start, len(source)):
        if source[index] == "{":
            depth += 1
        elif source[index] == "}":
            depth -= 1
            if depth == 0:
                return source[start:index]
    raise ValueError("Unclosed title")


def main():
    source = CANONICAL.read_bytes()
    text = source.decode()
    title = prose(title_argument(text))
    assert title == "Butson towers and high-rank Hadamard matrices"
    abstract = prose(text.split(r"\begin{abstract}", 1)[1].split(r"\end{abstract}", 1)[0])
    (HERE / "main.tex").write_bytes(source)
    command = ["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", "-outdir=build", "main.tex"]
    result = subprocess.run(command, cwd=HERE, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    build_dir = HERE / "build"
    (build_dir / "latexmk-output.txt").write_text(result.stdout)
    if result.returncode:
        print(result.stdout)
        raise SystemExit(result.returncode)
    log = (build_dir / "main.log").read_text()
    warnings = [line for line in log.splitlines()
                if "Warning:" in line or "Overfull" in line or "Underfull" in line]
    if warnings:
        raise ValueError("Resolve typesetting warnings: " + "; ".join(warnings))
    page_match = re.search(r"Output written on .*?\((\d+) pages?", log, re.DOTALL)
    if not page_match:
        raise ValueError("Page count missing from TeX log")
    pages = int(page_match.group(1))
    pdf = (build_dir / "main.pdf").read_bytes()
    (HERE / "announcement.pdf").write_bytes(pdf)
    build = {
        "source": "main.tex", "source_sha256": digest(source),
        "reference_pdf": "announcement.pdf", "pdf_sha256": digest(pdf), "pages": pages,
        "command": " ".join(command), "warnings": warnings,
        "source_and_reference_pdf_checked": True,
        "all_pages_visually_reviewed": False,
        "visual_review_note": "Set true only after inspecting all pages of this exact PDF.",
    }
    (HERE / "build.json").write_text(json.dumps(build, indent=2) + "\n")
    metadata_path = HERE / "metadata.template.json"
    metadata = json.loads(metadata_path.read_text())
    metadata.update(title=title, abstract=abstract)
    metadata["comments"] = (
        f"{pages} pages. Eighteen new real Hadamard ranks to the best of the author's knowledge; "
        "six new infinite output families; new periodic T-matrices of order 127. "
        "All announced ranks and six all-iteration families have complete Lean existence proofs "
        "from literal finite witnesses. Includes standalone Python/NumPy reconstruction, exact "
        "source and endpoint Gram proofs, an optional SymPy trace audit, and a pinned link to "
        "the separate source-only Lean reproducibility package. Known ranks are identified "
        "individually in the ancillary status records."
    )
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    subprocess.run([sys.executable, str(HERE / "prepare.py")], cwd=ROOT, check=True)
    print(json.dumps({"title": title, "pages": pages, "abstract_words": len(abstract.split()),
                      "warnings": len(warnings), "remote_requests_performed": False}))


if __name__ == "__main__":
    main()
