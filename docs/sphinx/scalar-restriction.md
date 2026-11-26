(scalar-restriction)=
# Restrict scalars while retaining the norm

Restriction of scalars replaces a scalar by its multiplication matrix. The
source norm then has a transported matrix identity. A final Hadamard matrix
also needs the required entries and a Gram scale equal to its order. This
page derives the transport identities, specifies their conventions, and
connects them to the exact field and endpoint implementations.

The concrete T-matrix construction in {ref}`construction` supplies two
examples: a large rational action of the sixth-root field, and an action
of a rational group algebra. Neither requires intermediate multiplication
matrices to have equal entry magnitudes.

(restriction-conventions)=
## Declare the source, involution and dimensions

Let $K\subseteq L\subseteq\mathbb C$ be finite extensions stable under
complex conjugation. Write $\bar z$ for conjugation and
$X^*=\overline X^{\mathsf T}$ for the conjugate transpose of a matrix.
For real matrices this is ordinary transpose.

A complete source in this page is a specified square matrix
$N\in M_n(L)$ with an exact identity

```{math}
:label: eq-restriction-source
NN^*=\omega I_n,\qquad \omega\in K\cap\mathbb R_{>0}.
```

Its order is $n$ and its Gram scale is $\omega$. A Butson matrix
$N\in\operatorname{BH}(n,q)$ has all entries in
$\mu_q=\{e^{2\pi i k/q}:0\leq k<q\}$ and $\omega=n$. A complex weighing
matrix permits zero entries and has its separately declared weight
$\omega$. A weighted algebra-valued source can have another entry
alphabet; {ref}`construction-group-algebra` specifies one exactly.
The assumption $\omega\in K$ ensures $\rho(\omega)=\omega I_d$ and permits
pulling the scale through the relative trace. Without it, the right side of
the row identity is $I_n\otimes(\rho(\omega)G^{-1})$, rather than a scalar
multiple of the metric. All implemented construction scales are positive
integers.

Choose a basis $\beta_0,\ldots,\beta_{d-1}$ of $L$ over $K$, where
$d=[L:K]$. Coordinates are columns:

```{math}
u=\sum_a x_a\beta_a,\qquad
z\beta_b=\sum_a\rho(z)_{ab}\beta_a.
```

Thus $\rho(z)$ acts on the coordinate column of $u$. Applying both sides
to each basis vector proves

```{math}
\rho(z+t)=\rho(z)+\rho(t),\qquad
\rho(zt)=\rho(z)\rho(t),\qquad
\rho(k)=kI_d\quad(k\in K).
```

The action is faithful: multiplication by $z$ sends $1$ to $z$.
Its block extension to matrices is denoted by the same symbol:
$X=\rho(N)$ has $(i,j)$ block $\rho(N_{ij})$ and order $nd$.

These are the regular multiplication coordinates. A separate
representation can have dimension larger than $[L:K]$; that distinction
is treated in {ref}`restriction-representation-dimension`.

(restriction-trace-adjoint)=
## Compute the trace form and the adjoint

Define the relative trace metric

```{math}
:label: eq-restriction-trace-metric
G_{ab}=\operatorname{Tr}_{L/K}(\bar\beta_a\beta_b).
```

For coordinate columns $x,y$ of $u,v$, expansion gives

```{math}
x^*Gy=\operatorname{Tr}_{L/K}(\bar u v).
```

Conjugation commutes with the relative trace: conjugating a basis changes
the multiplication matrix of $z$ into a matrix for $\bar z$ with
conjugated entries, and matrix trace is independent of the basis. Hence
$G^*=G$. The trace pairing for a finite characteristic-zero field extension
is nondegenerate. One proof uses a primitive element and its distinct
conjugates: the embedding matrix in the power basis is a nonsingular
Vandermonde matrix, and the trace-pairing matrix is its transpose times
itself. Any other basis changes that matrix by an invertible congruence.
Composing one argument with the invertible map $u\mapsto\bar u$ preserves
nondegeneracy. Thus $G$ is invertible.
For a more general trace functional on an algebra, invertibility is an
explicit hypothesis and is checked by the implementation.

For any $z\in L$,

```{math}
\begin{aligned}
x^*\rho(z)^*Gy
 &=\operatorname{Tr}_{L/K}(\overline{zu}\,v)\\
 &=\operatorname{Tr}_{L/K}(\bar u\bar z v)
 =x^*G\rho(\bar z)y.
\end{aligned}
```

Taking coordinate basis vectors for $x,y$ proves

```{math}
:label: eq-restriction-adjoint
\rho(z)^*G=G\rho(\bar z),\qquad
\rho(\bar z)=G^{-1}\rho(z)^*G,\qquad
\rho(\bar z)G^{-1}=G^{-1}\rho(z)^*.
```

These identities specify the adjoint in the chosen coordinates. In
particular, transposing a multiplication block without its metric is
justified only when the declared representation has that property.

The code {repo-file}`relative_scalar_restriction.py <experiments/scalar_lift/relative_scalar_restriction.py>`
builds the relative basis, its inverse coordinate change, every upper-root
multiplication block, $G$, and $G^{-1}$. `RelativeRestriction` checks the
complete root table by exact basis conjugacy, the cyclotomic polynomial,
the inverse-root relation, and the displayed adjoint identity.
`edge(upper_phase, lower_phase, basis_exponents)` returns this checked
object. `action(exponent)` returns a multiplication block; `proof()` and
`metadata()` retain the definitions and checks.

Coefficients are stored as integer numerators with a common denominator.
The calculations use exact rational arithmetic, and fixed-width expansion
steps have overflow guards. Acceptance is not based on a floating-point
tolerance.

```{literalinclude} ../../experiments/scalar_lift/relative_scalar_restriction.py
:language: python
:pyobject: RationalArray.exact
```

(restriction-row-column)=
## Use the correct form for row and column identities

For the row Gram in {eq}`eq-restriction-source`, the transported identity
contains the **inverse** metric:

```{math}
:label: eq-restriction-row-gram
X(I_n\otimes G^{-1})X^*
  =\omega(I_n\otimes G^{-1}),\qquad X=\rho(N).
```

Indeed, its $(i,j)$ block is

```{math}
\begin{aligned}
\sum_k\rho(N_{ik})G^{-1}\rho(N_{jk})^*
 &=\sum_k\rho(N_{ik})\rho(\overline{N_{jk}})G^{-1}\\
 &=\rho\!\left(\sum_kN_{ik}\overline{N_{jk}}\right)G^{-1}
 =\omega\delta_{ij}G^{-1}.
\end{aligned}
```

The first equality is exactly {eq}`eq-restriction-adjoint`.
The Gram scale remains $\omega$, while the order becomes $nd$.

There is an equivalent column convention. Since $N$ is square and
$NN^*=\omega I_n$, its right inverse $N^*/\omega$ is also its inverse.
Therefore $N^*N=\omega I_n$, which transports to

```{math}
:label: eq-restriction-column-gram
X^*(I_n\otimes G)X=\omega(I_n\otimes G).
```

Its $(i,j)$ block is
$G\rho(\sum_k\overline{N_{ki}}N_{kj})=\omega\delta_{ij}G$.
Thus the row formula uses $G^{-1}$ and the column formula uses $G$.
The two are not competing normalizations.

`RelativeRestriction.proof()` records the row identity. The state object
`ExactNormState` in {repo-file}`exact_norm_descent.py <experiments/scalar_lift/exact_norm_descent.py>`
uses the column convention, storing an operator $M$ and a Hermitian form
$E$ with $M^*EM=nE$. Its `verify()` checks that exact identity:

```{literalinclude} ../../experiments/scalar_lift/exact_norm_descent.py
:language: python
:pyobject: ExactNormState.verify
```

Here the implementation's `source_order` is the original Butson order
$n$, which is also that source's norm scale. The generic symbol
$\omega$ above also covers weighing and other weighted sources.

(restriction-towers)=
## Transfer an incoming form through every stage

At a later stage the form need not be an identity matrix or have constant
entries. Let $E=(E_{ij})\in M_r(L)$ be Hermitian, and suppose
$T^*ET=\omega E$. Its restricted form is

```{math}
:label: eq-restriction-incoming-form
\widehat E_{(i,a),(j,b)}
  =\operatorname{Tr}_{L/K}(\bar\beta_a E_{ij}\beta_b).
```

If $u_i=\sum_a x_{ia}\beta_a$ and $v_j=\sum_b y_{jb}\beta_b$, expansion
proves $x^*\widehat E y=\operatorname{Tr}_{L/K}(u^*Ev)$. It follows that
$\widehat E$ is Hermitian and

```{math}
\rho(T)^*\widehat E\rho(T)=\omega\widehat E.
```

For computation, the $(i,j)$ block of the form is
$G\rho(E_{ij})$, since its $(a,b)$ entry is
$\operatorname{Tr}_{L/K}(\bar\beta_a E_{ij}\beta_b)$. Therefore

```{math}
:label: eq-restriction-form-product
\widehat E=(I_r\otimes G)\rho(E).
```

The factor $\rho(E)$ retains all incoming field coefficients. Replacing
it by an identity would change the form being transferred.

Now let $F\subseteq K\subseteq L$ be a tower of fields stable under
conjugation, with $\omega\in F\cap\mathbb R_{>0}$ when the same scalar
scale is retained at every stage. Let $(\gamma_b)$ be an $L/K$ basis and $(\beta_a)$ a $K/F$
basis. The product basis, ordered by $(b,a)$, is $(\gamma_b\beta_a)$.
If the $L/K$ multiplication matrix of $z$ has entries $c_{bd}(z)$,
its $F$-matrix has block $\rho_{K/F}(c_{bd}(z))$. Taking the diagonal
matrix traces gives

```{math}
\operatorname{Tr}_{L/F}(z)
 =\operatorname{Tr}_{K/F}\bigl(\operatorname{Tr}_{L/K}(z)\bigr).
```

Applying this equality to the entries of
{eq}`eq-restriction-incoming-form` gives

```{math}
\begin{aligned}
&\operatorname{Tr}_{K/F}\!\left(
 \bar\beta_a\beta_c\,
 \operatorname{Tr}_{L/K}(\bar\gamma_b E_{ij}\gamma_d)\right)\\
&\qquad=\operatorname{Tr}_{L/F}\!\left(
 \overline{\gamma_b\beta_a}E_{ij}\gamma_d\beta_c\right).
\end{aligned}
```

Thus repeated restriction and direct restriction have exactly the same
operator and form in the product basis. If another basis is used, the
actual basis change must be applied. For a basis-change matrix $S$, the
new multiplication action is $S^{-1}\rho(z)S$ and the new metric is
$S^*GS$.

`ExactNormState.restrict` implements {eq}`eq-restriction-form-product`,
checks the child norm, and records the parent and restriction hashes:

```{literalinclude} ../../experiments/scalar_lift/exact_norm_descent.py
:language: python
:pyobject: ExactNormState.restrict
```

The test
{repo-file}`test_nonconstant_incoming_metric_transfers_by_trace_transitivity <tests/test_exact_norm_descent.py>`
uses the path $\mathbb Q(\zeta_{30})\supset\mathbb Q(\zeta_6)\supset\mathbb Q$.
Its middle form has nonconstant coefficients. The test compares both the
operator and metric with direct restriction by their explicit basis change.

(restriction-products)=
## Retain every factor of a product algebra

The same adjoint proof works for a finite-dimensional commutative algebra
$A$ with involution and a specified trace functional $\tau$, provided
$\tau(\bar z)=\overline{\tau(z)}$ and
$G_{ab}=\tau(\bar\beta_a\beta_b)$ is invertible.
For a product $A=\prod_\alpha L_\alpha$ of fields whose factors are
preserved by conjugation, one can declare

```{math}
\tau((z_\alpha)_\alpha)
 =\sum_\alpha c_\alpha\operatorname{Tr}_{L_\alpha/K}(z_\alpha),
```

with specified nonzero weights $c_\alpha$ fixed by conjugation. In a
concatenated factor basis the metric is block diagonal, with blocks
$c_\alpha G_\alpha$. Repeated factors occur separately. The operator,
involution, trace weights and reconstruction from factor coordinates
remain part of the data.

For the group-algebra sources, the direct basis has a particularly simple
metric. Let $G_0$ be a finite abelian group of order $w$. With
$c_gc_h=c_{g+h}$, $c_g^*=c_{-g}$ and
$\tau(\sum_g a_gc_g)=a_0$, the inverse-group basis $\beta_g=c_{-g}$ gives

```{math}
\rho(a)_{u,v}=a_{v-u},\qquad
\tau(\beta_u^*\beta_v)=\delta_{u,v}.
```

This is the identity metric used in {ref}`construction-group-algebra`.
The regular matrix trace is $wa_0$, so the normalization of $\tau$ is
fixed rather than implicit.

For a cyclic group,

```{math}
\mathbb Q[C_w]\cong\mathbb Q[t]/(t^w-1)
 \cong\prod_{m\mid w}\mathbb Q(\zeta_m),\qquad
\sum_{m\mid w}\varphi(m)=w.
```

Here $\varphi$ is Euler's totient function. The factorization
$t^w-1=\prod_{m\mid w}\Phi_m(t)$ has pairwise coprime irreducible
cyclotomic factors over $\mathbb Q$, so the Chinese remainder theorem
gives this product. In those factor coordinates,

```{math}
\tau(a)=\frac1w\sum_{m\mid w}
             \operatorname{Tr}_{\mathbb Q(\zeta_m)/\mathbb Q}(a_m).
```

For example, $w=15$ gives field degrees $1,2,4,8$, summing to $15$;
the scalar algebra is not one degree-15 number field.
For a noncyclic group $G_0=\prod_i C_{n_i}$, the factors can have
multiplicity. The number of group elements of order dividing $m$ is
$\prod_i\gcd(n_i,m)$. Subtracting the counts at proper divisors gives
the number of elements of exact order $m$, and division by $\varphi(m)$
gives the number of corresponding rational cyclotomic factors. Indeed,
evaluation at the complex characters diagonalizes the group algebra.
Characters of exact order $m$ have values generating
$\mathbb Q(\zeta_m)$, and their rational Galois orbits have size
$\varphi(m)$. Grouping these characters into their orbits gives the
stated rational factors and multiplicities.

`algebra_description(w)` and `group_algebra_description(dimensions)` in
{repo-file}`williamson_scalar_endpoints.py <experiments/scalar_lift/williamson_scalar_endpoints.py>`
record these factors, dimensions and trace normalization and check their
dimension sum. The concrete regular action and its sign entries are
checked by the source and endpoint verifiers. The general cyclotomic
`ExactNormState` executor and these product-algebra constructions have
separate APIs; the common transport proof is the one above.

Individual factor matrices need not be flat. The final entry test is made
after reconstructing the full regular-action block.

(restriction-flatness)=
## Supply the frame, test the used palette, and match the scale

Suppose a unital rational representation
$\rho:A\to M_D(\mathbb Q)$ has the ordinary transpose adjoint
$\rho(\bar z)=\rho(z)^{\mathsf T}$, and the complete source has positive
rational Gram scale $\omega$. Let a rational matrix $F$ satisfy
$FF^{\mathsf T}=cI_D$ with $c>0$, and form blocks

```{math}
H_{ij}=\rho(N_{ij})F.
```

For every block row pair,

```{math}
\begin{aligned}
\sum_kH_{ik}H_{jk}^{\mathsf T}
 &=c\sum_k\rho(N_{ik})\rho(\overline{N_{jk}})\\
 &=c\rho\!\left(\sum_kN_{ik}\overline{N_{jk}}\right)
 =c\omega\delta_{ij}I_D.
\end{aligned}
```

Thus the output order is $nD$ and its Gram scale is $c\omega$.
For a real Hadamard endpoint, two further checks are required:

```{math}
:label: eq-restriction-flatness
(\rho(z)F)_{ab}^2=1
\quad\text{for every used scalar }z\text{ and every }a,b,
\qquad c\omega=nD.
```

A rational entry satisfies $t^2=1$ exactly when it is $1$ or $-1$.
These are finite polynomial entry conditions on the actual scalar action
and frame. Equal nonzero magnitudes can first be expressed with a common
scale; the displayed equations are their sign-normalized form.
Restriction and reconstruction can pass through nonflat matrices before
these final equations hold.

For a dense Butson source $\omega=n$, the scale condition is $c=D$.
For a weighted source it is the separate equation $c\omega=nD$.
A zero source entry always gives the zero block $\rho(0)F=0$, so a
weighing or conference source retains its zero positions under this
scalar-linear substitution. A verified weighing endpoint is described by
its weight, not assigned the expanded order as its Gram scale.

The set of scalars actually used must be stated. In
{ref}`construction-sixth-root`, the frame $F=P$ certifies
$\{1,\zeta_6,-1,-\zeta_6\}$. The remaining sixth-root image
$\rho(\zeta_6^2)P=2V$ is not a sign block. The source for order 7572
uses only the certified scalars $-1,\zeta_6$.

Returning to the regular $d$-dimensional action and metric $G$ of
{ref}`restriction-trace-adjoint`, one sufficient normalization is to supply
matrices $C,D_0\in M_d(K)$ with
$C^*C=G$ and $D_0D_0^*=cG^{-1}$. Then the blocks
$C\rho(z)D_0$ give an ordinary row Gram $c\omega I$ by
{eq}`eq-restriction-row-gram`. Their entries still require the palette
check. Such rational factors are additional data; a numerical square
root of $G$ is not a certificate that they exist over the required field.

There is a related whole-matrix interface at a rational endpoint. Starting
from $M^{\mathsf T}EM=nE$, supply nonsingular rational input and output
bases $S,T$, and put

```{math}
B=T^{-1}MS,\qquad E_{\rm in}=S^{\mathsf T}ES,\qquad
E_{\rm out}=T^{\mathsf T}ET.
```

Substitution gives $B^{\mathsf T}E_{\rm out}B=nE_{\rm in}$.
The sufficient ordinary-norm conditions
$E_{\rm out}=cI$ and $E_{\rm in}=acI$, with $c>0$, give
$B^{\mathsf T}B=naI$. `ExactNormState.endpoint` checks these supplied
bases, the weighted identity, the claimed amplification $a$, and every
entry. Its `sign_Hadamard_verified` result requires both the norm
normalization and zero entry defects. In the Butson restriction workflow,
$a$ is the accumulated relative degree, so $na$ is the expanded order.

(restriction-representation-dimension)=
## Separate the field degree from the representation dimension

For $K=\mathbb Q(\zeta_6)$ the regular representation has dimension two.
The representation constructed from a T-design instead acts on
$\mathbb Q^d$, with $d=4p$, by
$\rho(\zeta_6)=R=QP^{\mathsf T}/d$.
Since $K$ is a field and the representation is unital and faithful,
$\mathbb Q^d$ is a vector space over $K$ and

```{math}
d=[K:\mathbb Q]\dim_K\mathbb Q^d
  =2\dim_K\mathbb Q^d.
```

For order 7572, $p=631$, $d=2524$ and this module dimension is $1262$.
It is not the field extension degree.

The actual Euclidean trace form in those coordinates is explicit. Put
$S_0=2R-I_d$ and $\xi=\sqrt{-3}$. The frame identities imply
$S_0^{\mathsf T}=-S_0$ and $S_0^2=-3I_d$. Define

```{math}
h(u,v)=\frac{u^{\mathsf T}v}{2}
       -\frac{u^{\mathsf T}S_0v}{6}\,\xi.
```

Skew symmetry gives $h(v,u)=\overline{h(u,v)}$. If
$s=u^{\mathsf T}v$ and $t=u^{\mathsf T}S_0v$, then
$h(u,S_0v)=t/2+s\xi/2=\xi h(u,v)$. Together with rational bilinearity,
this proves $K$-linearity in the second argument and conjugate linearity
in the first. Finally,

```{math}
\operatorname{Tr}_{K/\mathbb Q}h(u,v)=u^{\mathsf T}v,
\qquad h(u,u)=\tfrac12u^{\mathsf T}u>0\quad(u\ne0).
```

Thus the trace metric of this module in its actual rational coordinates
is $I_d$. The two-dimensional power-basis metric below describes a
different coordinate realization of the same scalar field.

(restriction-eisenstein)=
## Check the two-dimensional sixth-root example

Let $\zeta=\zeta_6=(1+\sqrt{-3})/2$. In the ordered basis $(1,\zeta)$,

```{math}
A=\rho(\zeta)=\begin{pmatrix}0&-1\\1&1\end{pmatrix},\qquad
G=\begin{pmatrix}2&1\\1&2\end{pmatrix},\qquad
G^{-1}=\frac13\begin{pmatrix}2&-1\\-1&2\end{pmatrix}.
```

The first matrix follows from $\zeta\cdot1=\zeta$ and
$\zeta^2=\zeta-1$. The metric follows from
$\operatorname{Tr}(1)=2$, $\operatorname{Tr}(\zeta)=1$ and
$\bar\zeta\zeta=1$. It satisfies
$A^{\mathsf T}GA=G$ and $AG^{-1}A^{\mathsf T}=G^{-1}$.
It does not satisfy $AA^{\mathsf T}=I_2$.

These exact values can be inspected through the same API used in a tower:

```python
from sympy import Matrix, Rational
from experiments.scalar_lift.relative_scalar_restriction import edge, regular_matrix

restriction = edge(6, 2)
A = regular_matrix(restriction.action(1), 2)
G = regular_matrix(restriction.metric, 2)
assert A == Matrix([[0, -1], [1, 1]])
assert G == Matrix([[2, 1], [1, 2]])
assert regular_matrix(restriction.inverse_metric, 2) == G.inv()
assert A.T * G * A == G
```

{repo-file}`The relative-restriction tests <tests/test_relative_scalar_restriction.py>`
also substitute a complete order-three sixth-root source. They check its
six-dimensional weighted row identity exactly and check that the unframed
matrix is not a sign Hadamard matrix of order six.

The `phase` arguments are declared root labels. The labels one and two
both have rational coefficient field, but designate different root
alphabets. Consequently `edge(6,1)` and `edge(6,2)` retain distinct lower
root labels even though their rational multiplication and metric arrays
coincide. A field conductor and an allowed root alphabet are separate data.

(restriction-example-tower)=
## Replay a complete tower and a nonflat norm control

`exact_norm_descent.positive_control()` starts with the Fourier source of
order four, written as $N_{ij}=\zeta_8^{2ij}$. This is a complete
$\operatorname{BH}(4,8)$ source whose entries already lie in the fourth
roots. The declared field is $\mathbb Q(\zeta_8)$, and the tested path is

```{math}
\mathbb Q(\zeta_8)\supset\mathbb Q(i)\supset\mathbb Q.
```

The relative degrees are two and two. The path basis has exponent order
$(0,2,1,3)$, while the direct power basis has order $(0,1,2,3)$.
The code checks both action conjugacy and metric congruence under the
explicit permutation. It then supplies the input basis
$I_4\otimes(H_2\otimes H_2)$, with
$H_2=\left(\begin{smallmatrix}1&1\\1&-1\end{smallmatrix}\right)$,
and identity output basis, and verifies a known sign Hadamard matrix of
order sixteen.

```{literalinclude} ../../experiments/scalar_lift/exact_norm_descent.py
:language: python
:pyobject: positive_control
```

A separate test shows why the entry check is necessary even after ordinary
norm normalization. Starting from the order-four Fourier source over
$\mathbb Q(i)$, it uses $H_2$ as the input frame and the rational orthogonal
output change

```{math}
Q_0=\frac15\begin{pmatrix}3&-4\\4&3\end{pmatrix}.
```

The ordinary Gram scale is still eight. But a block corresponding to the
source scalar one is

```{math}
Q_0^{\mathsf T}H_2
 =\frac15\begin{pmatrix}7&-1\\-1&-7\end{pmatrix},
```

which is not a sign matrix. The exact norm is preserved while terminal
flatness fails. Both the successful tower and this negative entry control
are in {repo-file}`test_exact_norm_descent.py <tests/test_exact_norm_descent.py>`.

```sh
PYTHONPATH=src:. .venv/bin/python -m pytest -q \
  tests/test_relative_scalar_restriction.py tests/test_exact_norm_descent.py
```

These controls establish the implementation's transport and endpoint
checks. They make no claim to new Hadamard orders.
