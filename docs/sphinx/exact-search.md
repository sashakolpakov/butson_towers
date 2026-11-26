(exact-search)=
# Exact search for the finite input

The search supplies four disjoint periodic sequences. The construction then
uses those sequences in exact matrix identities and in restriction of scalars.
These are distinct steps: a search program supplies no new Hadamard matrix
until it returns a literal witness that passes the integer verifier and the
endpoint checks. The original successful length-631 input and its resulting
order-7572 matrix belong to the retained construction results. The searches at
lengths 223, 167 and 179 have not yet supplied verified new inputs.

The method used successfully at length 631 is retained:

```text
unsigned coordinates satisfying every parity equation
    → affine sign equations modulo 4
    → quadratic sign equations modulo 8
    → exact coordinate sums and integer correlations
    → independent expansion at every cyclic shift
```

The newer native, SAT and CP-SAT implementations change how the sign equations
are solved. Their common mathematical constraints are described below; their
encodings and return statuses are explained in {ref}`sign-solvers`.

(exact-search-object)=
## The object being sought

Let $p$ be an odd prime and let $e_0,e_1,e_2,e_3$ denote the standard basis of
$\mathbb Z^4$. A candidate is a function

```{math}
b:\mathbb Z/p\mathbb Z\longrightarrow
\{e_0,-e_0,e_1,-e_1,e_2,-e_2,e_3,-e_3\}.
```

Write $b_j(x)$ for its coordinate $j$. Exactly one of the four integers
$b_j(x)$ is nonzero at each residue $x$. Define the periodic correlation

```{math}
C(s)=\sum_{x\in\mathbb Z/p\mathbb Z} b(x)^{\mathsf T}b(x+s).
```

The required identity is $C(0)=p$ and $C(s)=0$ for every nonzero $s$. To see
why this gives periodic T-matrices, set
$T_j(u,v)=b_j(v-u)$, with indices reduced modulo $p$. Then

```{math}
\left(\sum_{j=0}^3T_jT_j^{\mathsf T}\right)_{u,v}
=\sum_{j=0}^3\sum_x b_j(x)b_j(x+u-v)
=C(u-v)=p\delta_{u,v}.
```

The supports of the $T_j$ are disjoint and cover every position because the
same statement holds for the four coordinates of $b(x)$. Thus both support
and orthogonality follow from the stated signed-basis alphabet and correlation
identity. The subsequent disjoint orthogonal design and scalar restriction
still require their own exact identities and final entry check; see the
{ref}`construction`.

### Sites, multiplication orbits and gauges

The unrestricted chart has one unknown signed basis vector at every residue.
The function
{repo-file}`periodic_chart <experiments/scalar_lift/unrestricted_parity_map.py>`
records the actual site pairs for shifts $1,\ldots,(p-1)/2$. These shifts
suffice during solving because $C(s)=C(-s)$, by substituting $x-s$ for $x$.
The final verifier nevertheless checks all $p$ shifts.

The successful 631 chart imposes an additional symmetry. Let $\Pi$ fix
$e_0$ and cycle $e_1,e_2,e_3$. For an element $\mu$ of order $m$ modulo $p$,
the chart requires

```{math}
b(\mu x)=\Pi b(x).
```

Here $3\mid m$. A choice at one representative determines the values on its
multiplication orbit. For $(p,m,\mu)=(631,21,458)$ there are 31 orbits,
including $\{0\}$, and 15 nonzero shift classes under multiplication by
$\mu$ and negation. For $(223,3,183)$ these numbers are 75 and 37. The
{repo-file}`chart and expansion code <experiments/scalar_lift/hadamard_frame_search.py>`
stores each residue's orbit index and its power of the coordinate cycle.
Changing an orbit representative therefore requires transporting that cycle
power, as implemented by `large_parity_map.normalize_head`.

Both chart types fix $b(0)=e_0$. The symmetric sign-lift interface additionally
requires the orbit represented by 1 to use coordinate zero and fixes its
sign. At length 223 the fixed sign is negative; at the retained 127 and 631
controls it is positive. These are recorded restrictions of the solver's
input. The unrestricted model fixes no second site unless `fixed_signs` is
supplied explicitly.

(exact-search-parity)=
## Unsigned correlations as a determinant equation

Encode a coordinate label by $a=u+v\omega\in\mathbb F_4$, where
$u,v\in\mathbb F_2$ and $\omega^2+\omega+1=0$. For
$z=u+v\omega$,

```{math}
z^3=u+v+uv\in\mathbb F_2.
```

The value is zero precisely when $z=0$. Consequently the equality indicator
of two labels is

```{math}
[a_i=a_j]=1+(a_i+a_j)^3.
```

It is a Boolean polynomial of degree at most two in their binary coordinates.
In a symmetric chart one label is first multiplied by the appropriate power
of $\omega$; this is a binary linear operation and preserves that degree
bound. Signs disappear modulo two, so the unsigned correlation parities form
a quadratic map
$F:\mathbb F_2^{4d}\to\mathbb F_2^d$ for the charts considered here.

The coordinate recovery constructs invertible binary maps $T$ and $L$ with

```{math}
F\bigl(T(P,Q,U,V)\bigr)=L(PV+QU+1),
\qquad P,Q,U,V\in\mathbb F_{2^d}.
```

The implementation recovers the scalar field from ratios of polar forms and
constructs two hyperbolic planes. It then checks the resulting map directly
against the original site-pair correlations. Field reconstruction alone is
not accepted as validation.

Since $L$ is invertible, all parity conditions hold exactly when
$PV+QU=1$. This is the determinant-one condition in characteristic two, so
the entire unsigned parity locus is parametrized by
$\operatorname{SL}_2(\mathbb F_{2^d})$. It imposes no sign or exact-count
condition yet.

### Why the coefficient check proves the whole identity

A Boolean quadratic has a unique squarefree expression consisting of its
constant, linear and pairwise quadratic coefficients. Values on zero, every
unit vector and every sum of two distinct unit vectors determine all those
coefficients. For example, the coefficient of $z_i z_j$ is
$F(e_i+e_j)+F(e_i)+F(e_j)+F(0)$ in characteristic two.

Thus checking
$1+4d+\binom{4d}{2}$ vector coefficients proves the displayed identity on
every input. The rank checks also prove that both coordinate maps are
invertible. This is a finite polynomial-identity check, rather than evidence
from random samples. The unrestricted implementation makes the probes and
rank conditions explicit:

```{literalinclude} ../../experiments/scalar_lift/unrestricted_parity_map.py
:language: python
:pyobject: validate_map
```

A separate `validate_expanded_geometry` reconstructs each probe as a
length-$p$ coordinate word and compares the parity at every nonzero shift.
It checks that the compressed pair geometry agrees with literal cyclic
correlations. The retained certificates record the following checks:

| Chart | Field degree $d$ | Binary coordinate variables | Coefficient probes | Expanded nonzero-shift comparisons |
|---|---:|---:|---:|---:|
| Length 223, multiplier order 3 | 37 | 148 | 11,027 | 2,447,994 |
| Unrestricted length 167 | 83 | 332 | 55,279 | 9,176,314 |
| Unrestricted length 179 | 89 | 356 | 63,547 | 11,311,366 |

The source data are the
{repo-file}`223 map <certificates/scalar-lift/t223-parity-map.json>`,
{repo-file}`167 map <certificates/scalar-lift/unrestricted-parity-p167.json>`
and {repo-file}`179 map <certificates/scalar-lift/unrestricted-parity-p179.json>`.
The original 631 map has its own complete coefficient check in
`hadamard_frame_search.validate_parity_map`.

### Sampling the proven locus

Choose $(P,Q)\ne(0,0)$ uniformly and choose $W$ uniformly in the same field.
If $P\ne0$, set $U=W$ and $V=(1+QU)/P$. If $P=0$, set $U=Q^{-1}$ and
$V=W$. Each determinant-one matrix has exactly one such description. The
unconditioned procedure therefore samples uniformly among the
$2^d(2^{2d}-1)$ matrices in this finite group.

Filtering by support sizes or normalizing a representative changes the
sampling distribution being reported. In particular, representative
normalization is not a claim of uniform sampling in a quotient. The Python
implementation documents this distinction:

```{literalinclude} ../../experiments/scalar_lift/large_parity_map.py
:language: python
:pyobject: DeterminantSampler.draw_packed
```

A validated parametrization describes all parity solutions. A sequence of
draws from it is still a sampled search, even if every sign assignment is
checked for each drawn coordinate pattern.

(exact-search-modular-lifting)=
## Affine signs modulo four

Fix the occupied coordinates and write the remaining sign at index $i$ as
$\sigma_i=1-2x_i$, with $x_i\in\{0,1\}$. An index is an actual site in
the unrestricted chart and a multiplication orbit in a symmetric chart.
For shift class $\ell$, let $D_\ell$ be its sign-independent diagonal
contribution and let $w_{ij,\ell}$ count its equal-coordinate pairs between
distinct indices. The exact correlation is

```{math}
C_\ell=D_\ell+\sum_{i<j}w_{ij,\ell}\sigma_i\sigma_j
=A_\ell-2\sum_{i<j}w_{ij,\ell}(x_i\mathbin{\oplus}x_j),
\qquad
A_\ell=D_\ell+\sum_{i<j}w_{ij,\ell}.
```

The parity stage ensures that $A_\ell$ is even. Since
$x_i\oplus x_j\equiv x_i+x_j\pmod2$, the requirement
$C_\ell\equiv0\pmod4$ is the binary linear equation

```{math}
\sum_i\left(\sum_{j\ne i}w_{ij,\ell}\bmod2\right)x_i
\equiv A_\ell/2\pmod2.
```

The model appends the fixed signs and the parities of the exact coordinate
counts. Binary Gaussian elimination then gives

```{math}
x=x^{(0)}\mathbin{\oplus}By,\qquad y\in\mathbb F_2^k.
```

`particular` represents $x^{(0)}$; `basis` contains the columns of $B$;
`forms[i]` is the packed binary linear form giving the change in $x_i$.
These masks are Python integers, so the model construction has no
machine-word dimension limit. The shared elimination is implemented by
`quadratic_sign_lift._eliminate`; the unrestricted model builds its input
rows as follows:

```{literalinclude} ../../experiments/scalar_lift/unrestricted_t_sign_model.py
:language: python
:start-at:     rows = []
:end-at:     particular, basis = _eliminate(rows,n)
```

An inconsistent linear system excludes the fixed coordinate pattern with its
specified counts and fixed signs. It does not exclude other coordinate
patterns, multiplier charts or Hadamard constructions.

## Quadratic signs modulo eight

For a pair $i,j$, put
$\varepsilon_{ij}=(-1)^{x_i^{(0)}+x_j^{(0)}}$ and let
$L_{ij}(y)=\bigoplus_u a_{ij,u}y_u$ be the XOR of their affine forms.
As an integer-valued Boolean function,

```{math}
L_{ij}(y)\equiv\sum_u a_{ij,u}y_u
-2\sum_{u<v}a_{ij,u}a_{ij,v}y_uy_v\pmod4.
```

Terms of degree three and above have coefficients divisible by four. Using
$\sigma_i\sigma_j=\varepsilon_{ij}(1-2L_{ij})$ therefore gives

```{math}
C_\ell\equiv c_\ell-2\sum_u S_{\ell,u}y_u
+4\sum_{u<v}T_{\ell,uv}y_uy_v\pmod8,
```

where
$c_\ell=D_\ell+\sum_{i<j}\varepsilon_{ij}w_{ij,\ell}$,
$S_{\ell,u}=\sum_{i<j}\varepsilon_{ij}w_{ij,\ell}a_{ij,u}$ and
$T_{\ell,uv}=\sum_{i<j}\varepsilon_{ij}w_{ij,\ell}a_{ij,u}a_{ij,v}$.
The affine modulo-four conditions ensure $4\mid c_\ell$ and
$2\mid S_{\ell,u}$. Division by four is consequently exact and yields the
Boolean quadratic equation

```{math}
q_\ell(y)=c_\ell/4-\sum_u(S_{\ell,u}/2)y_u
+\sum_{u<v}T_{\ell,uv}y_uy_v=0\quad\text{in }\mathbb F_2.
```

This is precisely $C_\ell/4\pmod2$ on the affine space, not an approximation
to the integer correlation. The code checks the divisibilities before
exporting any equation:

```{literalinclude} ../../experiments/scalar_lift/quadratic_sign_lift.py
:language: python
:pyobject: quadratic_mod8
```

(exact-search-counts)=
## Exact sums and final verification

Let $r_j=\sum_x b_j(x)$ be the four coordinate sums. Summing all
correlations and changing the order of summation gives the necessary identity

```{math}
\sum_s C(s)=\sum_j\left(\sum_x b_j(x)\right)^2
=\sum_j r_j^2=p.
```

If coordinate $j$ has support size $N_j$, its required number of negative
entries is exactly

```{math}
n_{j,-}=(N_j-r_j)/2.
```

This number must be an integer between zero and $N_j$. Fixed signs may make
the count impossible or force every remaining sign in that coordinate.
These consequences are applied before elimination. The modulo-four model
uses count parities; a final witness must satisfy the complete integer
counts. The unrestricted builder retains four separate groups and four
targets in `groups` and `negative_counts`.

In the symmetric 223 chart, the sums are $(-14,3,3,3)$. Indeed symmetry makes
them $(a,t,t,t)$ with $a\equiv1\pmod3$, while
$a^2+3t^2=223$ forces $a=-14$ and $t=\pm3$; changing the signs of the last
three coordinates permits $t=3$. If $h$ nonzero orbits use coordinate zero,
their signs sum to $-5$, and the other $74-h$ orbit signs sum to 3. Hence
$h$ is odd, $5\le h\le71$, and

```{math}
(n_{0,+},n_{0,-},n_{\ne0,+},n_{\ne0,-})
=\left(\frac{h-5}{2},\frac{h+5}{2},
        \frac{77-h}{2},\frac{71-h}{2}\right).
```

At $h=5$, every coordinate-zero orbit sign is forced negative. A negative
such orbit always exists; rescaling residues can place it at 1. At $h=7$
and $h=9$, after that fixed negative sign, the complete positive-sign
branches number $\binom61=6$ and $\binom82=28$. Keeping all head signs
within the affine system is another complete way to solve the same fixed
coordinate pattern; it avoids explicit branch multiplication.

The unrestricted 167 and 179 models instead use the four individual site
counts. Representative prepared models retain 80 and 86 affine bits. The
631 retained control has 13 free bits; the current symmetric 223 models
typically have 32--35. These are dimensions of particular reduced models,
not new existence results or a universal rank claim for every input.

No modular survivor is accepted directly. The final verifier reconstructs a
fresh integer $4\times p$ array from the signed coordinate word, checks its
alphabet and origin, checks every one of the $p$ correlations, and checks the
requested sums. It does not reuse the sign solver's pair tables:

```{literalinclude} ../../experiments/scalar_lift/unrestricted_t_sign_model.py
:language: python
:pyobject: verify_states
```

## Choose a multiplier before choosing the matrix order

A large multiplier can make a large prime much cheaper than a smaller
unrestricted input. In this chart the multiplier has odd order $m=3k$,
where $k$ is odd. Write $N=(p-1)/m$ for the number of nonzero multiplication
orbits. There are $N+1$ orbit-sign variables and $N/2$ nonzero correlation
classes. If the modulo-four correlation rows have rank $N/2$ and the fixed
signs and count parities add three independent conditions, the nullity is
$(N+1)-N/2-3=N/2-2$. Exact count constraints can force more signs; the actual
dimension must be measured for each prepared pattern.

Let $a$ and $b$ be the signed sums of the nonzero coordinate-zero orbits
and the remaining orbits, respectively. Changing the last three coordinate
signs permits $b>0$. The four expanded coordinate sums are
$(1+ma,kb,kb,kb)$. Their squared norm gives

```{math}
N=2a+k(3a^2+b^2),\qquad p=3kN+1,\qquad a\equiv b\pmod2.
```

This is a finite arithmetic way to select charts. It also gives the bound
$p\leq1+3N(N+2)/4$, proved in the
{repo-file}`compact scan note <docs/HADAMARD_FRAME_TARGET_SCAN.md>`.
The {repo-file}`replayable scanner <experiments/scalar_lift/multiplier_target_scan.py>`
agrees with an independent prime-first enumeration. At $N\le160$ it finds
87 parameter choices on 79 primes, recorded in the
{repo-file}`finite parameter certificate <certificates/scalar-lift/multiplier-target-scan-20261126.json>`.
Each surviving chart still proceeds through sign lifting and the complete
integer correlation checks.

```{list-table}
:header-rows: 1
:widths: 10 15 15 16 44

* - $p$
  - Multiplier order $m$
  - Sign variables
  - Free bits
  - State
* - 631
  - 21
  - 31
  - 13
  - Completed input
* - 919
  - 27
  - 35
  - 15
  - Verified map; no completed input
* - 1453
  - 33
  - 45
  - 20 or 21
  - Verified map, including the zero-head case; no completed input
* - 1951
  - 39
  - 51
  - 23 in checked patterns
  - Verified map; active target
* - 3571
  - 51
  - 71
  - 33 in checked patterns
  - Verified map; no completed input
* - 4219
  - 57
  - 75
  - 35 in checked patterns
  - Verified map; no completed input
```

For $p=1951$, the
{repo-file}`coordinate map <certificates/scalar-lift/hadamard-frame-p1951.json>`
uses $\mu=1127$ of order 39 and the field $\mathbb F_{2^{25}}$. All 5,051
constant, linear and quadratic coefficients, and all 9,849,450 independently
expanded parity correlations, were checked. The requested sums are
$(-38,13,13,13)$. A complete input would give the candidate real orders
23412, 54628 and 101452 through the retained source orders 3, 7 and 13.
Source order 9 would give 70236, already a Paley order.

The larger maps at
{repo-file}`3571 <certificates/scalar-lift/hadamard-frame-p3571.json>` and
{repo-file}`4219 <certificates/scalar-lift/hadamard-frame-p4219.json>` also
pass every quadratic coefficient probe at every expanded nonzero shift:
35,239,470 and 46,511,886 parity comparisons. Four prepared patterns per
chart leave 33 and 35 free sign bits. Their 71 and 75 orbit-sign variables
require the modern arbitrary-size masks. Larger multiplier orders have
therefore not made these two searches cheaper than 631.

There is also a finite limitation on finding another equally small chart
with this coordinate action. The scanner finds exactly 16 charts with at
most 30 nonzero orbits, the number used at 631. Every other prime in that
list already has Williamson-type matrices. Combining those with classical
T-inputs at our 27 small source widths recovers known Hadamard orders.
The {repo-file}`target-selection note <docs/HADAMARD_FRAME_TARGET_SCAN.md>`
gives the complete list and primary references. Broader coordinate actions
can supply different small searches without altering the subsequent scalar
restriction construction.

Even multiplier orders cannot provide a shortcut in this prime-length
signed-permutation model. The
{repo-file}`inversion obstruction <docs/EVEN_MULTIPLIER_OBSTRUCTION.md>`
proves this for every prime greater than three, including actions that
exchange two coordinates at inversion. The identity action remains a
separate possibility from the three-cycle action used at 631.

The compact $p=397$ chart is an engine control. Its formerly listed target
14292 is known: Paley gives a symmetric conference matrix of order 198;
Seberry–Yamada's Corollary 29(ii) then gives symmetric Williamson-type
matrices of order 397, which Cooper--Wallis combines with T9. The
{repo-file}`dated records update <docs/HADAMARD_HIGHER_ORDER_RECORDS.md>`
gives the primary references and the other known-order gates used in target
selection.


(exact-search-replay)=
## Small deterministic replay

Run the following from the repository root in the project's Python
environment. It constructs and proves the complete small-field parity map;
it does not launch a stochastic search or a frontier computation.

```sh
PYTHONPATH=src:. .venv/bin/python - <<'PY'
from experiments.scalar_lift import unrestricted_parity_map as maps

config = maps.recover_map(7, seed=7)
identity = maps.validate_map(config)
geometry = maps.validate_expanded_geometry(config)
assert identity["coefficients_checked"] == 79
assert geometry["exact_correlations_compared"] == 474
print("GF(8): complete quadratic identity and all probe shifts verified")
PY
```

The {repo-file}`unrestricted model tests <tests/test_unrestricted_t_sign_model.py>`
also compare the entire length-7 pipeline with independent enumeration of
all $8^6$ signed words having positive origin. They recover 180 solutions
with sums $(2,1,1,1)$, from 504 parity coordinate patterns. The
{repo-file}`original finite controls <tests/test_quadratic_sign_lift.py>`
compare all 90 stated symmetric T19 coordinate patterns and their 150 exact
sign solutions. These are complete finite controls in their declared
alphabets and gauges; neither is a classification at lengths 167, 179 or 223.
