# Formal proof status and remaining kernel

The Lean project in [`formal`](../formal) contains checked general
construction theorems. It does **not yet formally certify the concrete
T631/T127 witnesses, every source table, or the completed endpoint list**.
Those objects retain their separate exact executable certificates.

## Checked now

[`ScalarFrame.lean`](../formal/ButsonTowers/ScalarFrame.lean) proves over a
commutative rational star algebra, which may have zero divisors:

* The inverse-metric adjoint identity from a metric and its two-sided inverse.
* Weighted pair, sum, and Gram transport for an algebra representation.
* Flat-frame pair and sum transport, and the **full expanded matrix Gram**.
* Preservation of zero entries and a separate used-entry flatness implication.
* A positive-order Hadamard conclusion that explicitly requires
  `frame scale * incoming Gram scale = expanded matrix order`.

In the Euclidean case, for a rational algebra homomorphism `rho` with
`rho(star a)=transpose(rho a)` and a frame `F F^T=cI`, the proved identity is

\[
\sum_k(\rho(a_k)F)(\rho(b_k)F)^T
=c\,\rho\left(\sum_k a_k\overline{b_k}\right).
\]

The theorem `framed_matrix_gram` expands this identity to all individual
rows indexed by the product of the source and coordinate index types.
`framed_matrix_hadamard` additionally requires literal sign images for the
actual source entries and the scale/order equality. Neither property is
inferred from the source Gram alone. The source need not be Butson.

[`Spence.lean`](../formal/ButsonTowers/Spence.lean) proves:

* Both doubled Spence norm recurrences in an arbitrary **commutative ring**.
  Doubling removes any implicit field or division-by-two assumption.
* The sum-of-four-squares identity from the two strong norms.
* Strong-norm induction, assuming the stated recurrence laws.
* The exact rational pointwise sign test for the actual half-sum recurrence.

These are sixteen proved theorems. No project axiom, `sorry`, `admit`, or
native-evaluation certificate is used. The exact Spence matrix/group-ring
realization and the concrete starter checks are not yet formalized.

## Build and trust record

The portable configuration pins Lean **4.32.0** and mathlib commit
`81a5d257c8e410db227a6665ed08f64fea08e997`. The checked build reused an
already installed matching mathlib cache; its local symlinks are ignored
build state and are absent from the portable dependency manifest.

```sh
cd formal
lake build
lake env lean AuditAxioms.lean
```

The second command prints the dependencies of all sixteen theorems.
The reported axioms are subsets of Lean's standard
`propext`, `Classical.choice`, and `Quot.sound`; there is no `sorryAx` or
project-specific axiom. This is a partial mathematical formalization,
not a claim that the Python verifiers or all recorded matrices have been
translated into Lean.

## Minimal remaining theorems

The following order keeps the trusted construction kernel small.

1. **Finite-group regular action.** For a finite abelian group `G`, define
   convolution and inversion on rational coefficient functions and prove
   `rho(f)[r,s]=f(s-r)` preserves addition, multiplication, the identity,
   and the involution. Prove the identity-coefficient trace gives metric
   `I` in the inverse-group basis. This directly instantiates the checked
   frame theorems for the Williamson scalar restrictions.
2. **T tuple to orthogonal design.** For four circulant coefficients prove
   the Goethals--Seidel Gram identity using their commutativity and the
   specified reversal matrix. Prove the literal four-variable rotation
   satisfies `L(x)^T L(x)=(sum x_j^2)I`. A periodic T norm then gives all ten
   OD coefficient identities; its one-hot support gives disjoint complete
   coefficient supports.
3. **OD substitution and flatness.** Prove substitution of self-adjoint
   commuting group-ring elements into the OD preserves the displayed
   weighted Gram. Prove that disjoint OD coefficient supports and sign
   coefficient polynomials give only sign entries after regular
   restriction. This connects the previous two steps to the existing
   expanded Hadamard theorem.
4. **Spence realization.** Define the tensor/group extension and identify
   the actual half-sum outputs with the doubled algebraic expressions
   already proved. Establish symmetry, commutativity and the two norm
   recurrences from `Q^2=qI-J`, `J^2=qJ`, `JQ=0`. The pointwise sign lemma
   and strong-norm induction then give every family member. The dimension
   changes at each step must be represented by the actual group extension
   or explicit algebra embeddings, rather than suppressed in the types.
5. **Concrete witnesses.** Define the T631/T127 data and Spence starter
   rows as literal Lean objects; prove their finite checkers sound, then
   prove each checker returns true. Search may supply witnesses but must
   not become an axiom. A stored hash binds files but is not a mathematical
   proof of their Gram. The NTT modulus/bound argument likewise needs a
   proved soundness theorem before its external result can replace a
   formal correlation check.
6. **General trace restriction and towers.** Given a finite free star
   algebra, a star-compatible nondegenerate trace functional and a basis,
   prove the multiplication adjoint relation for its trace Gram. Use the
   checked inverse-metric lemma to obtain row-Gram transport. Prove trace
   transitivity and basis congruence for successive restrictions. Extend
   rational transpose to conjugate transpose over a lower star field.
7. **Finite etale interpretation.** Prove the cyclic and finite-abelian
   group algebra decompositions into cyclotomic fields, including repeated
   factors, and identify the normalized trace and factor metrics. These
   results justify the field-factor descriptions; they are not needed to
   search again for any orthogonality identity.

## Semantic constraints to keep explicit

* **Algebra versus field.** `Q[G]` is generally a product of fields, not a
  single field. No cancellation of a possibly zero group-algebra element
  is permitted. `q=121` uses the additive group `C11 x C11`, not `C121`.
* **Metric direction.** With `rho(a)^T G=G rho(star a)`, row-Gram transport
  uses `G^-1`: `A(I tensor G^-1)A^T=w(I tensor G^-1)`. The column convention
  uses `G` and an incoming column Gram. These are distinct statements.
* **Positivity.** A finite etale algebra with an arbitrary involution need
  not have a positive trace form: swapping the factors of `Q x Q` gives
  an indefinite form. The Gram transport only needs the advertised
  nondegeneracy/adjoint identities. Positivity requires the correct CM or
  Euclidean form and a separate proof.
* **Trace normalization.** The group metric `I` uses normalized regular
  trace. Unnormalized regular trace gives `|G|I`. In cyclotomic factor
  power bases the normalized blocks are Ramanujan trace matrices divided
  by `|G|`; their multiplicities must be retained.
* **Dimension and scale.** The target action dimension is not automatically
  a field degree. The sextic action of dimension2524 has field degree2 and
  module rank1262. The incoming Gram scale is not multiplied by a field
  degree merely because coordinates are expanded. The final equality
  with expanded order is an explicit hypothesis of the checked theorem.
* **Flatness and zeros.** A weighted intermediate need not have root
  entries. Every used entry image needs its own exact alphabet test;
  zeros remain zero under a linear scalar restriction.
* **Noncommutative substitutions.** The four OD coefficient identities
  alone do not identify a Hamilton quaternion action. The present
  substitution proof uses commuting self-adjoint scalar coefficients.
* **Novelty.** Priority and author-knowledge status belong to the dated
  literature ledger, not to a mathematical existence theorem.
