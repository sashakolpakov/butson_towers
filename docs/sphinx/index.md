(manual)=
# Butson Towers

This manual follows the construction from a finite signed-coordinate input to
an exactly verified Hadamard matrix. Definitions, proofs, data layouts and
executable checks appear together at the step where they are needed.

[Read the manuscript (PDF)](https://sashakolpakov.github.io/butson_towers/butson-towers-paper.pdf) ·
[Download the printable manual (PDF)](https://sashakolpakov.github.io/butson_towers/butson-towers-manual.pdf).

Begin with {ref}`getting-started` to reconstruct a small complete example.
Then follow {ref}`scalar-restriction` and {ref}`construction` for the
identities behind the order-7572 construction and its later source reuse.
The subsequent chapters explain how inputs are found, how families remain
implicit at large orders, and how the executable and Lean checks establish
precise claims.

| Stage | Mathematical obligation | Executable object |
|---|---|---|
| Finite input | Disjoint signed basis values and zero nontrivial periodic correlations | Literal states, coordinate sequences and integer verifier |
| Disjoint design | Polynomial Gram identity and entrywise support partition | Four coefficient matrices and circulant/block recipes |
| Scalar source | Exact norm, involution, action and trace form | Butson exponents, group polynomials or explicit rational operators |
| Restriction and frame | Transport the norm and test the used entry alphabet | Exact rational action, sign blocks or regular representation |
| Completed endpoint | Sign entries and Gram scale equal to expanded order | Literal entry recipe, certificates and formal existence proof |

The project name describes chains of scalar restrictions. Its
{ref}`references` record the earlier morphism constructions and the immediate
starting reference. The current paper's completed results and the live
search targets have separate entries in {ref}`results`.

```{toctree}
:maxdepth: 2

getting-started
scalar-restriction
construction
families
exact-search
sign-solvers
campaigns
verification
results
reference
contributing
references
```
