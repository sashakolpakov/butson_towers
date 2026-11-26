(families)=
# Exact sources and infinite output families

The reusable finite ingredient is the T631 orthogonal design. Its four
coefficient matrices $X_0,\ldots,X_3$ have order $d=2524$, disjoint signed
supports covering every entry, and satisfy

```{math}
:label: eq-families-design
X_jX_j^{\mathsf T}=631I_d,\qquad
X_iX_j^{\mathsf T}+X_jX_i^{\mathsf T}=0\quad(i\ne j).
```

A source quadruple consists of sign matrices $A_j$ of a common order $w$
with

```{math}
:label: eq-families-source
A_iA_j^{\mathsf T}=A_jA_i^{\mathsf T},\qquad
\sum_{j=0}^3A_jA_j^{\mathsf T}=4wI_w.
```

These conditions are called pairwise amicability and the summed Gram
identity. Substitution gives

```{math}
:label: eq-families-endpoint
H=\sum_{j=0}^3X_j\otimes A_j,\qquad
HH^{\mathsf T}=2524wI_{2524w}.
```

For distinct $i,j$, the two cross terms combine to
$(X_iX_j^{\mathsf T}+X_jX_i^{\mathsf T})\otimes A_iA_j^{\mathsf T}=0$.
The diagonal terms give $631I_d\otimes4wI_w$. Disjointness of the $X_j$
proves the entry condition: exactly one signed source entry contributes at
each position, so $H$ is a sign matrix. These are separate arguments for
orthogonality and flatness.

The executable checks are in
{repo-file}`williamson_scalar_endpoints.py <experiments/scalar_lift/williamson_scalar_endpoints.py>`;
the general proof is in
{repo-file}`OrthogonalSubstitution.lean <formal/ButsonTowers/OrthogonalSubstitution.lean>`.

## Literal Williamson sources

The retained input table contains 27 symmetric circulant Williamson
quadruples. For first rows $a_j:\mathbb Z/w\mathbb Z\to\{1,-1\}$, the
verifier checks $a_j(-x)=a_j(x)$ and every coefficient of

```{math}
:label: eq-families-williamson
\sum_{j=0}^3\sum_{x\in\mathbb Z/w\mathbb Z}a_j(x)a_j(x+s)
=4w\delta_{s,0}.
```

Their circulants commute; symmetry then gives amicability. Consequently
these finite checks establish all hypotheses of
{eq}`eq-families-source`, rather than only some selected row products.
The literal rows, provenance and source hashes are in
{repo-file}`williamson_inputs.json <experiments/scalar_lift/williamson_inputs.json>`.
The paper prints thirteen of these cyclic quadruples and the twelve further
compatible inputs described below; the repository formally checks all 27
original cyclic inputs. A malformed entry labelled order 45 in the
pinned upstream table is excluded because its literal words have length 43.

From the repository root:

```python
from experiments.scalar_lift.williamson_scalar_endpoints import (
    verified_williamson, endpoint,
)

A, source_proof = verified_williamson(17)
H = endpoint(631, 17)
assert H.order == 42908
assert H.entry(0, 100) in (-1, 1)
```

`verified_williamson` checks the source identities and input binding.
`endpoint` also verifies the finite T-design coefficients. The resulting
object evaluates entries and rows without allocating the complete endpoint.
Its `.row(i)` method does materialize the requested row.

## Further compatible inputs and their scalar algebras

The {repo-file}`additional literal table <experiments/scalar_lift/near_williamson_inputs.json>`
contains 19 published inputs through width 181. Let $C(a)_{i,j}=a(j-i)$
and let $R_{i,j}=1$ when $i+j=w-1$, with zero otherwise. Except at width
69, the represented source is
$A_0=C(a_0)R$, $A_j=C(a_j)$ for $j=1,2,3$; at 69 all four are ordinary
circulants. The last three rows are symmetric. Thus $C(a_0)R$ is symmetric,
commutes with the other three matrices and has the same Gram contribution
as $C(a_0)$. Every source norm and pairwise identity is checked in full.

```python
from experiments.scalar_lift.near_williamson_scalar_endpoints import (
    verified_source, endpoint,
)

A, proof = verified_source(47)
H = endpoint(47)
assert H.order == 118628
assert H.entry(0, 100) in (-1, 1)
```

These matrices generate a commutative reduced algebra
$\mathcal A=\mathbb Q[A_0,A_1,A_2,A_3]$ on $V=\mathbb Q^w$, with physical
metric $I_w$. They are not assumed to belong to a cyclic group algebra.
The {repo-file}`exact algebra audit <experiments/scalar_lift/commutative_scalar_algebra.py>`
expresses every $A_j$ as a polynomial in $T=A_0+3A_1+5A_2+7A_3$, factors its
minimal polynomial, and reconstructs all field modules by rational spectral
idempotents. If $\mathcal A=\prod_iK_i$, it retains
$V=\bigoplus_iK_i^{m_i}$ and $\sum_im_i[K_i:\mathbb Q]=w$.
The normalized trace is $\sum_i m_i\operatorname{Tr}_{K_i/\mathbb Q}/w$.
For example, the width-73 input has algebra dimension 25 and module
dimension 73: its degree-24 field occurs three times. The audit preserves
all three copies. Its basis change $B$ transports the metric to
$G=B^{\mathsf T}B$, and all represented scalar actions are checked to be
self-adjoint for this $G$.

This is how nonflat field components return to literal sign blocks after
the recorded rational reconstruction. The
{repo-file}`algebra certificate <certificates/scalar-lift/near-williamson-scalar-algebras.json>`
stores compact reconstruction recipes and exact hashes, rather than large
rational matrices. Its optional factor audit requires `python-flint`:

```sh
python -m experiments.scalar_lift.commutative_scalar_algebra \
  --check certificates/scalar-lift/near-williamson-scalar-algebras.json
```

The {repo-file}`source inventory <docs/SOURCE_INVENTORY.md>` records which
twelve outputs add new orders and which outputs were already known.

## Source recurrences of Spence

We use Spence's classical recurrence, with the five parameters
$q\in\{13,17,61,101,121\}$ and $w_0=(q+1)/2$. The sources and their
recurrence are classical; the T631 construction supplies the new output
multiplier. The retained source attribution and precise starter data are
in {repo-file}`SPENCE_SCALAR_FAMILIES.md <docs/SPENCE_SCALAR_FAMILIES.md>`.

The separate
{repo-file}`additional executor <experiments/scalar_lift/additional_spence_families.py>`
also verifies prime parameters $5,53,97,109,113,137$ without changing the
earlier five-source records. The current paper states the three additional
families $7572\cdot5^r$, $174156\cdot137^r$ and $15131380\cdot109^r$.
The last starts at recurrence stage $r+1$; the other three prime recurrences
remain construction controls without a new-family claim. Each added starter
satisfies the same two strong identities below, and the same induction
proves every parameter value.

Write the four starting matrices as $D_0,E_0,F_0,G_0$. Besides being
symmetric commuting sign matrices, they must satisfy **both** identities

```{math}
:label: eq-families-strong-starter
(D_0-E_0)^2+(F_0-G_0)^2=4I_{w_0},\qquad
(D_0+E_0)^2+(F_0+G_0)^2=4qI_{w_0}.
```

An arbitrary Williamson quadruple has the summed norm, but need not have
these two stronger norms. The implementation's `starter(q)` checks signs,
symmetry, commutation and both equations. For $q=13,17$ it constructs the
starter from the recorded conference halves. For the other three values
it applies specified signed permutations to retained literal quadruples.

Let $Q$ be the symmetric quadratic-character matrix indexed by
$\mathbb F_q$, with $Q_{x,y}=\chi(y-x)$. Here $\chi(0)=0$, and on nonzero
elements $\chi$ is $+1$ for squares and $-1$ for nonsquares. Let $J$ be the
all-one matrix of order $q$. The required
finite identities are

```{math}
:label: eq-families-paley
Q^{\mathsf T}=Q,\quad Q^2=qI-J,\quad QJ=JQ=0,\quad J^2=qJ.
```

The diagonal of $Q$ is zero and its off-diagonal entries are signs. At
$q=121$ the field is $\mathbb F_{11}[u]/(u^2-2)$, with additive group
$C_{11}\times C_{11}$; it is not a cyclic group of order 121. Python checks
the complete finite matrix identities. Lean checks the corresponding
literal matrices in
{repo-file}`Data/Paley.lean <formal/ButsonTowers/Data/Paley.lean>` and its
finite-check modules.

Apply the following operation to each pair $(A,B)=(D_r,E_r)$ and
$(F_r,G_r)$:

```{math}
:label: eq-families-spence-step
\begin{aligned}
A'&=\frac{A-B}{2}\otimes J+\frac{A+B}{2}\otimes(Q+I),\\
B'&=\frac{A-B}{2}\otimes J+\frac{A+B}{2}\otimes(Q-I).
\end{aligned}
```

For paired sign entries $a,b$, exactly one of $(a-b)/2,(a+b)/2$ is a
sign and the other is zero. In the second case, $Q\pm I$ has sign entries
because $Q$ has zero diagonal. Thus the recurrence preserves signs.
Symmetry and commutativity follow from the corresponding properties of
the old matrices and of $I,J,Q$.

To prove the norm for every iteration, define

```{math}
:label: eq-families-strong-induction
\begin{aligned}
M_r&=(D_r-E_r)^2+(F_r-G_r)^2,\\
P_r&=(D_r+E_r)^2+(F_r+G_r)^2.
\end{aligned}
```

The difference of a new pair is $(A+B)\otimes I$; its sum is
$(A-B)\otimes J+(A+B)\otimes Q$. Squaring and using
{eq}`eq-families-paley` gives

```{math}
:label: eq-families-strong-recurrence
M_{r+1}=P_r\otimes I,\qquad
P_{r+1}=qP_r\otimes I+(qM_r-P_r)\otimes J.
```

Starting from {eq}`eq-families-strong-starter`, induction yields
$M_r=4q^rI$ and $P_r=4q^{r+1}I$. Hence the component order and norm are

```{math}
:label: eq-families-spence-size
w_r=w_0q^r,\qquad
D_r^2+E_r^2+F_r^2+G_r^2=\tfrac12(M_r+P_r)=4w_rI.
```

This is the invariant implemented by
{repo-file}`SpenceSource <experiments/scalar_lift/williamson_spence_family.py>`
and proved by `strong_step`, `family_strong` and `family_card` in
{repo-file}`SpenceMatrix.lean <formal/ButsonTowers/SpenceMatrix.lean>`.
The entry method recursively decodes the indices; one entry costs $O(r)$
arithmetic operations rather than a dense matrix construction.

## The scalar algebra and its dimension

For $q=\ell^e$, these components are group-circulant on
$G_r=C_{w_0}\times(C_\ell)^{er}$, of order $w_r$. Let
$a_0,\ldots,a_3\in\mathbb Q[G_r]$ denote their first-row group
polynomials. They are self-adjoint, have sign coefficients and satisfy
$\sum_j a_j^2=4w_rc_0$, where $c_0$ is the group identity. Therefore
$S_r=\sum_jX_ja_j$ is a complete weighted source of order 2524 over this
algebra and Gram scale $2524w_r$. Its entries are signed copies of the
$a_j$; they are not asserted to be roots of unity.

With $m_r=(q^r-1)/(\ell-1)$, the algebra decomposes as

```{math}
:label: eq-families-field-factors
\mathbb Q[G_r]\cong
\prod_{t\mid w_0}\mathbb Q(\zeta_t)
\;\times\;
\prod_{t\mid w_0}\mathbb Q(\zeta_{\ell t})^{m_r}.
```

To count the factors, put $E=(C_\ell)^{er}$. It has $q^r$ complex
characters. The trivial character gives a rational factor. Each nontrivial
character takes values in the $\ell$th roots, and the $\ell-1$ nonzero
power maps act freely on it. Thus the nontrivial rational character orbits
give $(q^r-1)/(\ell-1)=m_r$ copies of $\mathbb Q(\zeta_\ell)$.
Meanwhile $\mathbb Q[C_{w_0}]$ has one $\mathbb Q(\zeta_t)$ for each
$t\mid w_0$. Since $\gcd(\ell,w_0)=1$, tensoring each such factor with
$\mathbb Q(\zeta_\ell)$ gives $\mathbb Q(\zeta_{\ell t})$: their
compositum contains both roots and has the product of the two field degrees.
This proves the displayed multiplicities.

The total rational dimension is
$w_0[1+m_r(\ell-1)]=w_0q^r=w_r$. In particular $w_r$ is an algebra
dimension, not the degree of a single field. The involution sends
$c_g$ to $c_{-g}$. Use the normalized regular trace

```{math}
:label: eq-families-trace
\tau(a)=[c_0]a=\frac{1}{w_r}\sum_{\text{field factors }K}
\operatorname{Tr}_{K/\mathbb Q}(a_K).
```

Each repeated factor occurs in this sum. In the inverse-group basis
$\beta_g=c_{-g}$, one has
$\tau(\beta_u^*\beta_v)=\delta_{u,v}$ and
$\rho(a)_{u,v}=a_{v-u}$. Multiplication therefore sends each $a_j$ to its
literal sign group-circulant matrix. Restriction retains the source Gram
scale $2524w_r$ and expands its order from 2524 to $2524w_r$.
Disjointness of the $X_j$ proves the final signs. The complete reconstruction,
including every factor and multiplicity, is essential to this entry test.

| $q$ | Source component order $w_r$ | Hadamard order for every $r\geq0$ |
|---:|---:|---:|
| 13 | $7\cdot13^r$ | $17668\cdot13^r$ |
| 17 | $9\cdot17^r$ | $22716\cdot17^r$ |
| 61 | $31\cdot61^r$ | $78244\cdot61^r$ |
| 101 | $51\cdot101^r$ | $128724\cdot101^r$ |
| 121 | $61\cdot121^r$ | $153964\cdot121^r$ |

An implicit endpoint can be evaluated directly:

```python
from experiments.scalar_lift.williamson_spence_family import SpenceSource
from experiments.scalar_lift.williamson_scalar_endpoints import (
    verified_t_design, ImplicitScalarEndpoint,
)

X, design_proof = verified_t_design(631)
source = SpenceSource(13, 2)
H = ImplicitScalarEndpoint(X, source)
assert H.order == 17668 * 13**2
assert H.entry(0, 100) in (-1, 1)
```

`.entry` stays implicit; `.row` and `.coefficient_rows()` allocate the
requested rows. The finite `certificate()` replay checks the first
expanded source for each $q$, using an exact integer NTT with a proved
residual bound. The all-$r$ conclusion uses the induction above and its
Lean proof, rather than extrapolation from those five finite replays.

## The sixth family: unreal sixth-root tripling

Write a current unreal source of order $n$ as $U=(A+\sqrt{-3}B)/2$,
where $A,B$ are sign matrices of order $n$ satisfying the two identities
below. The retained initial source has $n=13$.

```{math}
:label: eq-families-unreal-pair
AA^{\mathsf T}+3BB^{\mathsf T}=4nI,\qquad
AB^{\mathsf T}=BA^{\mathsf T}.
```

These identities imply $UU^*=nI$. Put

```{math}
:label: eq-families-tripling
C=\begin{pmatrix}0&1&-1\\-1&0&1\\1&-1&0\end{pmatrix},\quad
A'=B\otimes J_3,\quad B'=A\otimes I_3+B\otimes C.
```

Here $J_3$ is the all-one matrix of order three, $I_3$ is the identity
matrix of that order, $C^{\mathsf T}=-C$, $CC^{\mathsf T}=3I-J_3$ and
$CJ_3=J_3C=0$. The diagonal/off-diagonal supports in $B'$ are disjoint,
so both new matrices have sign entries. Direct multiplication gives

```{math}
:label: eq-families-tripling-norm
A'A'^{\mathsf T}+3B'B'^{\mathsf T}=12nI,\qquad
A'B'^{\mathsf T}-B'A'^{\mathsf T}
=(BA^{\mathsf T}-AB^{\mathsf T})\otimes J_3=0.
```

Thus $n_r=13\cdot3^r$ is a valid source order for every $r\geq0$.
The source entries are the unreal sixth roots with exponents
$\{1,2,4,5\}$. A common factor $\zeta_6^{-1}$ changes that palette to
$\{0,1,3,4\}$. For the T631 action/frame, precisely these entries have
images $P,Q,-P,-Q$. Their sign entries and frame scale 2524 give the sixth
output family, $32812\cdot3^r$.

The Python recurrence follows the formula literally:

```{literalinclude} ../../experiments/scalar_lift/tripled_unreal_scalar_endpoints.py
:pyobject: triple
:language: python
```

`verified_source(n)` in that module fully materializes and verifies the
retained source orders $13,39,117,351$; `endpoint(39)` gives literal access
to order 98436. The finite executor deliberately has those supported
orders. The proof for arbitrary $r$ is in
{repo-file}`UnrealTripling.lean <formal/ButsonTowers/UnrealTripling.lean>`;
`unreal_cubic_family` in
{repo-file}`Endpoints.lean <formal/ButsonTowers/Endpoints.lean>` uses the
amicable quadruple $(A_r,B_r,B_r,B_r)$ for the final substitution.

## What is new

The source recurrences are established constructions; the six displayed
output families use the new T631 ingredient. Each has a new initial order
to the best of the author's knowledge. Individual members can have prior
constructions: for example 885924 in the cubic family and 18629644 in the
$q=121$ family are known orders. The exact finite labels and prior
constructions are recorded in
{repo-file}`HADAMARD_HIGHER_ORDER_RECORDS.md <docs/HADAMARD_HIGHER_ORDER_RECORDS.md>`.
No equivalence-class conclusion follows merely from two recipes having
the same order.
