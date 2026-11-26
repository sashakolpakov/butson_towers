(verification)=
# Verification, finite witnesses and Lean proofs

The existence proof starts from literal finite data, not from confidence in
a search. For a periodic T-input, the essential finite equality is

```{math}
:label: eq-verification-T
\sum_{j=0}^3\sum_{x\in\mathbb Z/p\mathbb Z}b_j(x)b_j(x+s)
=p\delta_{s,0}.
```

Together with the signed coordinate alphabet, this feeds the orthogonal
design and source-substitution proofs. The repository checks these
obligations in several ways. Each has a specific object and scope.

## Begin with the literal data

From the repository root, the following checks do not repeat a search:

```sh
python3 scripts/verify_manuscript_data.py
python3 scripts/check_lean_literal_data.py
```

The first uses the Python standard library and carries the complete literal
data in its own {repo-file}`source <scripts/verify_manuscript_data.py>`.
It verifies the T127/T631 functions, the thirteen printed Williamson
quadruples, the five strong starters, the finite character matrices and
the order-13 pair. All indices and inputs in those finite equalities are
specified. The second checks agreement between the public construction
data and the Lean literals. Its name does not mean that it invokes Lean.

The detailed proof correspondence is in
{repo-file}`MANUSCRIPT_LEAN_AUDIT.md <docs/MANUSCRIPT_LEAN_AUDIT.md>`.
It records the checked manuscript tables and the exact conventions used
by the formal construction.

## Different exact checks prove different claims

| Check | What is actually verified | Scope of the conclusion |
|---|---|---|
| Alphabet and support check | Each column of $b$ is one signed standard basis vector | The four coordinate supports partition every position |
| Integer correlation check | Every shift in {eq}`eq-verification-T` | The complete periodic T-input |
| Polynomial Gram check | All four norm coefficients and six mixed coefficients of $D(x)D(x)^{\mathsf T}$ | The full orthogonal-design identity for arbitrary substituted scalars satisfying the required relations |
| Source check | Signs, amicability and summed Gram of the actual source matrices | The particular source used by the substitution |
| Dense endpoint check | Every diagonal norm and every distinct row-pair inner product of a materialized matrix | The particular stored or reconstructed representative |
| Selected-row control | Entry/row API agreement and some full row products | The evaluated interfaces; the full Gram still follows from the construction proof |
| Lean kernel proof | The stated theorem, its literal hypotheses and its dependencies | Exactly that formal statement, including all parameters quantified by it |

For the design coefficients $X_j$ and source quadruple $A_j$, the full
algebraic endpoint proof is

```{math}
:label: eq-verification-substitution
\begin{aligned}
H&=\sum_jX_j\otimes A_j,\\
HH^{\mathsf T}
&=pI_{4p}\otimes\sum_jA_jA_j^{\mathsf T}
=4pwI_{4pw}.
\end{aligned}
```

The mixed terms vanish by the proved coefficient and source identities.
A separate disjoint-support argument proves every entry is a sign. An
exact polynomial Gram proof establishes every endpoint row product even
when the dense matrix is too large to allocate.

The implementation is
{repo-file}`williamson_scalar_endpoints.py <experiments/scalar_lift/williamson_scalar_endpoints.py>`.
Its `verify_orthogonal_design`, `verify_williamson` and `row_controls`
functions correspond to different rows of the table above. A returned
row-control record is not presented as a dense audit.

### The saved order-7572 representative

The standalone
{repo-file}`H7572 recipe <artifacts/hadamard-7572/recipe.py>` reconstructs
the fixed matrix and checks all $\binom{7572}{2}=28\,663\,806$ distinct
row pairs with exact Hamming distances. For two sign rows of length $N$,

```{math}
:label: eq-verification-hamming
\langle u,v\rangle=N-2\operatorname{Ham}(u,v).
```

Thus distance $N/2$ proves orthogonality exactly. The diagonal norm is
$N$ because every entry is a sign. This audit supplements the source and
frame proofs for that representative.

```sh
python3 artifacts/hadamard-7572/recipe.py
```

The bundle also contains independent sequence, scalar-frame and packed
matrix audits. See
{repo-file}`the bundle README <artifacts/hadamard-7572/README.md>` and
{repo-file}`independent-full-npz-audit.json <artifacts/hadamard-7572/independent-full-npz-audit.json>`.
The large family endpoints instead have exact entry recipes and complete
algebraic Gram proofs; the documentation does not claim dense audits for
every one of them.

### Exact finite transform checks

The first expanded Spence sources are checked coefficient by coefficient
using an integer number-theoretic transform. This is not a floating-point
Fourier residual. The retained modulus is $998244353$ and the checker
establishes a strict integer residual bound $|r|\leq8w<998244353$.
A zero residue therefore forces the integer residual to be zero.
Multivariate padding and folding implement the declared finite group,
including $C_{11}\times C_{11}$ at character order 121.

The algorithm and its bounds are in
{repo-file}`exact_group_norms <experiments/scalar_lift/williamson_spence_family.py>`;
the recorded finite checks are in
{repo-file}`williamson-spence-family.json <certificates/scalar-lift/williamson-spence-family.json>`.
The general family theorem comes from induction, not from a finite
transform experiment.

## What the Lean project proves

`ButsonTowers.HadamardExists n` asserts a positive order $n$ and the
existence of a rational matrix on `Fin n` with entries $1,-1$ and Gram
$nI_n$. The named endpoint theorems and all six family statements are in
{repo-file}`Endpoints.lean <formal/ButsonTowers/Endpoints.lean>`.
In particular the family parameters are universally quantified natural
numbers.

The dependency chain connects:

1. Literal finite T and source data to sound integer checkers and their
   rational matrix meanings.
2. Circulants, the Goethals–Seidel Gram identity and all orthogonal-design
   coefficient identities.
3. Disjoint signed supports, amicable substitution and equality of final
   order and Gram scale.
4. The five strong Spence recurrences and the unreal cubic recurrence,
   including their signs, norms and index cardinalities for every iteration.

The finite certificates use kernel reduction (`decide +kernel`). No
external search result, source-file hash, native evaluator or numerical
residual is a mathematical hypothesis of these endpoint proofs.

The regular group-algebra action, its faithfulness, involution and
normalized trace metric are proved in
{repo-file}`FiniteGroupAlgebra.lean <formal/ButsonTowers/FiniteGroupAlgebra.lean>`.
The endpoint proofs use concrete tensor substitution; they do not depend
on a general cyclotomic decomposition theorem. The trace and inverse-metric
identities have their own scope in
{repo-file}`TraceTowers.lean <formal/ButsonTowers/TraceTowers.lean>` and
{repo-file}`ScalarFrame.lean <formal/ButsonTowers/ScalarFrame.lean>`.

### The 95 selected axiom reports

Run the current repository audit with

```sh
./scripts/check_formal.sh
```

The script builds the project and checks two exact lists: 85 statements
in {repo-file}`AuditAxioms.lean <formal/AuditAxioms.lean>` and ten additional
rational sixth-root matrix identities in
{repo-file}`AuditSixthRootPair.lean <formal/AuditSixthRootPair.lean>`.
These are **95 selected axiom reports**, not a count of every theorem in
the project. The lists include the eighteen announced orders and all six
all-iteration families. Only `propext`, `Classical.choice` and `Quot.sound`
are permitted; missing, repeated or unexpected reports are rejected.

The project pins Lean 4.32.0, mathlib revision
`81a5d257c8e410db227a6665ed08f64fea08e997` and the nine transitive dependency
revisions in the Lake manifest. Reproduction retains that lock. A matching
precompiled dependency cache can be reused after its revisions and tracked
files are checked; the package also documents a dependency source-build
route. These checks do not rebuild the Lean compiler itself.

### Formal scope and representative conventions

The paper's character identities have a general finite-field proof. Lean
checks the five required finite matrices from literals. Its endpoint
proofs therefore need no general quadratic-character theorem.

The core formal GS reversal uses $u+v=0$; Python uses $u+v=-1$.
{repo-file}`ReflectionConvention.lean <formal/ButsonTowers/ReflectionConvention.lean>`
proves the reindexing relation at periods 127 and 631. The formal existence
statements do not identify their constructed matrix with the saved packed
H7572 representative.

The ten-statement
{repo-file}`SixthRootPair.lean <formal/ButsonTowers/SixthRootPair.lean>`
supplement derives the rational action identities from the paired-matrix
hypotheses. It does not instantiate a cyclotomic field or the manuscript's
particular Hermitian form. The full transfer of an arbitrary incoming
Hermitian matrix is proved in the manuscript; the existing formal trace
module has the narrower named trace-metric and tower statements described
in the correspondence audit. None of these distinctions leaves an
assumption missing from an announced endpoint theorem.

## Releases, checksums and reproducible runs

The frozen submission has a standalone executable source package and a
separate source-only Lean package. They serve different purposes:

| Package or record | Binding and replay |
|---|---|
| {repo-file}`Submitted source archive <release/butson_towers_v3/source.zip>` | Literal data, reconstruction recipes, finite source/trace replays and the manuscript |
| {repo-file}`Submission manifest <release/butson_towers_v3/submission-manifest.json>` | Archive-member sizes and hashes |
| {repo-file}`Build record <release/butson_towers_v3/build.json>` | The reference PDF and exact TeX source |
| {repo-file}`Complete Lean source packet <release/lean_verification/README.md>` | Pinned formal sources, literal data, dependency lock and the original 85-report verifier |
| {repo-file}`Sixth-root supplement <release/lean_sixth_root_pair/README.md>` | The additional ten matrix identities and their separate verifier |
| {repo-file}`Manuscript/Lean audit <docs/MANUSCRIPT_LEAN_AUDIT.md>` | Current interpretation of the exact theorem lists and their mathematical scope |

The older Lean packet's narrative says “five Spence families”; its exact
audited list also includes `unreal_cubic_family`. The certified total is
six output families. Frozen older prose is preserved, and the current
scope is stated explicitly here and in the correspondence audit.

```sh
python3 release/butson_towers_v3/verify.py
```

This extracts a fresh source archive, verifies its bindings and runs the
seven standalone executable replays. It checks the separate Lean packet's
hash; it does not itself rerun the Lean kernel. For a fresh kernel replay,
extract the appropriate Lean packet and run its `verify.py`, following
its README and pinned-cache instructions.

A source digest identifies bytes. The H7572 packed-bit digest identifies
the ordered sign data, whereas the NPZ digest also depends on its container.
The Lean records include `payload_manifest_sha256` and the exact
`axioms_by_theorem` map. These bindings make a result traceable to its
inputs; a digest alone does not establish a Gram identity or a theorem.

For a new local audit, retain a transcript separately from the structured
verification report:

```sh
set -o pipefail
./scripts/check_formal.sh 2>&1 | tee /tmp/butson-formal.log
```

After the command succeeds, its transcript can be bound by a digest:

```python
from hashlib import sha256
from pathlib import Path

log = Path('/tmp/butson-formal.log')
print(sha256(log.read_bytes()).hexdigest())
```

That hash identifies the captured run output; the successful build and
strict theorem-report parser establish what was checked. Current search
logs and elapsed counters are not substituted for either proof.

The complete executable regression suite runs with
`./scripts/run_tests.sh`. New T167/T179/T223 search inputs remain separate
from the frozen existence proofs until a literal witness and its required
new formal bridge are added and actually checked.
