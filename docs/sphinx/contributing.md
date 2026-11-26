(contributing)=
# Build and extend the manual

The documentation source is MyST Markdown read by Sphinx. Chapters follow the
construction and its verification steps. Add an implementation discussion at
the mathematical step it implements, with its required invariants, input
representation and replay evidence nearby.

## Reproduce the documentation build

From the repository root:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[docs]' -c docs/requirements.txt
./scripts/build_docs.sh html
./scripts/build_docs.sh doctest
```

HTML is written to `docs/_build/html/index.html`. The build treats warnings,
missing references and invalid source links as errors. The doctest builder
executes the explicitly marked small examples. Ordinary displayed shell
commands are documentation; they do not start long searches during a build.

The frontend versions tested for this manual are pinned in
{repo-file}`docs/requirements.txt <docs/requirements.txt>`. Runtime arithmetic
is supplied by the package's base dependencies. Add `.[search]` to execute
CP-SAT examples, and install the `cryptominisat5` executable for native SAT.

For the printable manual, install XeLaTeX, `latexmk`, the Latin Modern fonts
and DejaVu Sans Mono. Then run:

```sh
./scripts/build_docs.sh pdf
```

The resulting book is
`docs/_build/latex/butson-towers-manual.pdf`. `./scripts/build_docs.sh all`
builds HTML, runs doctests and creates this PDF. A separate interpreter can
be selected with `BUTSON_DOCS_PYTHON`; this allows documentation builds to
use their own environment while research jobs retain frozen dependencies.

```sh
BUTSON_DOCS_PYTHON=/path/to/docs-env/bin/python ./scripts/build_docs.sh html
```

The optional `linkcheck` target checks external URLs. Its network results are
separate from the deterministic local source/reference and doctest checks.

## Publish on GitHub Pages

The public manual is at
[sashakolpakov.github.io/butson_towers](https://sashakolpakov.github.io/butson_towers/),
with a [downloadable PDF](https://sashakolpakov.github.io/butson_towers/butson-towers-manual.pdf).
The {repo-file}`documentation workflow <.github/workflows/docs.yml>` runs after
relevant pushes to `main` and can also be started from the Actions tab. It
installs the pinned frontend, builds HTML, executes the small examples,
builds the PDF and publishes only `docs/_build/html`. The PDF is copied into
that directory for download. Sphinx's build caches stay outside the site.

Deployment uses the repository's `github-pages` environment and GitHub's
short-lived workflow token. It needs no repository API secret or local
submission configuration. The repository's Pages source is GitHub Actions;
the workflow's official actions are pinned to release commits.

## Add a mathematical construction

State the input's algebra, involution, norm scale and actual alphabet.
Introduce every new symbol and define the indexing convention before using
it. Derive the norm and entry properties separately. When the final order
has grown, show why its ordinary Gram scale equals that order.

Beside a code example, identify which objects it materializes and which stay
implicit. Name the exact arithmetic type, the returned data contract, the
scope of exceptions/statuses, and the corresponding proof or test. A few
checked small examples are more informative than an expensive stochastic
run whose result depends on a timeout.

A completed new input needs a literal witness and independent full checks.
It then needs endpoint data, the dated literature comparison and a new Lean
literal proof if formal existence is claimed. Update {ref}`results`, the
source ledger and relevant API examples together. Preserve the submitted
manuscript and frozen release hashes as historical records.

## Keep excerpts and references attached to actual code

Use a literal include when a short implementation is itself informative:

````markdown
```{literalinclude} ../../experiments/scalar_lift/williamson_scalar_endpoints.py
:pyobject: ScalarEndpoint
:language: python
```
````

Use the local `repo-file` role for public files that should remain links
rather than copied text:

```text
{repo-file}`endpoint verifier <experiments/scalar_lift/williamson_scalar_endpoints.py>`
```

The role checks that a target exists within the repository and rejects local
credentials, generated research stores and build products. The output link
points to the public repository. API directives import selected modules
without invoking their CLI entrypoints. Optional native solver calls remain
inside functions, so an ordinary book build does not start a solver.

Keep `*.local`, `.local/`, active ledgers and credentials outside the source
publication. Sphinx reads only its source directory and explicitly selected
public excerpts; generated output is ignored by Git. The new manual does not
repackage private submission configuration or modify frozen construction
archives.

## Validation commands and their scopes

| Check | What it establishes |
|---|---|
| `./scripts/build_docs.sh html` | Source parses, all local cross-references and selected source targets resolve |
| `./scripts/build_docs.sh doctest` | Marked small Python examples execute with the documented results |
| `python docs/examples/verify_small_endpoint.py` | Literal T7 input, polynomial design, width-three source and complete dense H84 Gram |
| `./scripts/run_tests.sh` | The repository's executable arithmetic/construction/search regression suite |
| `./scripts/check_formal.sh` | Lean build and the selected axiom audits, when their pinned dependencies are present |
