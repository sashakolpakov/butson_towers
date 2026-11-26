(construction)=
# From a periodic T-design to a Hadamard matrix

A finite signed-coordinate function supplies the input to this construction.
Its exact correlation identity gives a four-variable orthogonal design. Two
substitutions then produce real Hadamard matrices: the sixth-root source
used for order 7572, and verified commuting quadruples interpreted in a
rational group algebra. This page proves the identities and identifies the
code that checks their hypotheses.

The general rule for transporting a norm and checking the final entries is
proved in {ref}`scalar-restriction`. Here every matrix is real until the
sixth-root field or group algebra is introduced. Transpose is denoted by
$^{\mathsf T}$; $I_r$ is the identity matrix of order $r$.

(construction-t-witness)=
## Verify the finite input

Let $p$ be a positive integer and write $G_p=\mathbb Z/p\mathbb Z$ additively.
Let $e_0,e_1,e_2,e_3$ be the standard basis of $\mathbb Z^4$. The input is a
function

```{math}
b:G_p\longrightarrow\{\pm e_0,\pm e_1,\pm e_2,\pm e_3\}
```

whose summed periodic correlation is

```{math}
:label: eq-construction-t
C_b(s):=\sum_{u\in G_p}b(u)^{\mathsf T}b(u+s)
       =p\delta_{s,0}\qquad(s\in G_p).
```

Here $\delta_{s,0}$ is one when $s=0$ and zero otherwise. Equivalently, the
four coordinate functions $b_j:G_p\to\{0,1,-1\}$ have exactly one nonzero
entry at each residue and satisfy the displayed identity. These are the
periodic T-sequences used in the repository. The corresponding T-matrices
are defined below.

The input has to satisfy every shift. Reducing the number of variables by a
multiplier symmetry is a search device. For example, the retained
length-631 witness uses multiplication by $458$, of order $21$ modulo $631$,
together with the permutation fixing coordinate zero and cycling the other
three coordinates. The reconstruction expands this description to four
literal length-631 integer sequences before checking the construction.

The source records are {repo-file}`the length-631 witness <artifacts/hadamard-7572/t631.json>`
and {repo-file}`the length-127 witness <artifacts/t-matrices-127/t127.json>`.
`expand_t` in {repo-file}`the scalar endpoint module <experiments/scalar_lift/williamson_scalar_endpoints.py>`
checks the orbit cover and reconstructs the entries. `verify_t_rows` accepts
expanded integer sequences and checks their alphabet, the partition of
positions, and all $p$ correlations:

```{literalinclude} ../../experiments/scalar_lift/williamson_scalar_endpoints.py
:language: python
:pyobject: verify_t_rows
```

Thus the subsequent construction does not assume that a sampler, a modular
solver, or a numerical search found a genuine witness. It checks the finite
identity that the proof uses.

(construction-circulants)=
## Turn correlations into a matrix identity

For a real function $f:G_p\to\mathbb R$, define

```{math}
L_f(u,v)=f(v-u),\qquad f^-(u)=f(-u),\qquad
(f*g)(z)=\sum_{u\in G_p}f(u)g(z-u).
```

Directly multiplying matrix entries gives

```{math}
L_f^{\mathsf T}=L_{f^-},\qquad
L_fL_g=L_{f*g}=L_gL_f.
```

Indeed, the $(u,v)$ entry of the product is
$\sum_x f(x)g(v-u-x)$. Replacing $x$ by $v-u-x$ interchanges $f$ and $g$.
In particular, all these circulants and their transposes commute.
Their row Gram entries are

```{math}
:label: eq-construction-circulant-gram
(L_fL_f^{\mathsf T})_{u,v}
 =\sum_x f(x)f(x+v-u).
```

Consequently $T_j=L_{b_j}$ satisfies

```{math}
\sum_{j=0}^3T_jT_j^{\mathsf T}=pI_p.
```

At position $(u,v)$ the four entries $T_j(u,v)$ are the coordinates of the
same signed basis vector $b(v-u)$. Their supports are therefore disjoint
and cover all matrix positions. This proves the T-matrix assertions.

Define the reversal matrix by

```{math}
J_{u,v}=\delta_{u+v,p-1}.
```

Indices are residues modulo $p$. Its permutation $u\mapsto p-1-u$ is an
involution, so $J^2=I_p$ and $J^{\mathsf T}=J$. Moreover,

```{math}
(JL_fJ)_{u,v}=f((p-1-v)-(p-1-u))=f(u-v),
\qquad JL_fJ=L_f^{\mathsf T}.
```

For four functions $f=(f_0,f_1,f_2,f_3)$ put
$A=L_{f_0}$, $B=L_{f_1}$, $C=L_{f_2}$ and $D=L_{f_3}$. The exact
Goethals--Seidel array convention used here is

```{math}
:label: eq-construction-gs
\operatorname{GS}(f)=
\begin{pmatrix}
 A&BJ&CJ&DJ\\
 -BJ&A&-D^{\mathsf T}J&C^{\mathsf T}J\\
 -CJ&D^{\mathsf T}J&A&-B^{\mathsf T}J\\
 -DJ&-C^{\mathsf T}J&B^{\mathsf T}J&A
\end{pmatrix}.
```

Each diagonal Gram block is
$N=AA^{\mathsf T}+BB^{\mathsf T}+CC^{\mathsf T}+DD^{\mathsf T}$.
Using $JXJ=X^{\mathsf T}$, the six upper off-diagonal blocks simplify to

```{math}
\begin{array}{c|l}
(0,1)&(BA-AB)J+DC-CD\\
(0,2)&(CA-AC)J+BD-DB\\
(0,3)&(DA-AD)J+CB-BC\\
(1,2)&(AD^{\mathsf T}-D^{\mathsf T}A)J+BC^{\mathsf T}-C^{\mathsf T}B\\
(1,3)&(C^{\mathsf T}A-AC^{\mathsf T})J+BD^{\mathsf T}-D^{\mathsf T}B\\
(2,3)&(AB^{\mathsf T}-B^{\mathsf T}A)J+CD^{\mathsf T}-D^{\mathsf T}C.
\end{array}
```

Every difference is zero because the circulants commute. The lower blocks
are their transposes. We have therefore proved

```{math}
:label: eq-construction-gs-gram
\operatorname{GS}(f)\operatorname{GS}(f)^{\mathsf T}
 =I_4\otimes\sum_{j=0}^3L_{f_j}L_{f_j}^{\mathsf T}.
```

The Kronecker product is indexed by pairs:
$(X\otimes Y)_{(i,a),(j,b)}=X_{ij}Y_{ab}$. Summing independently over the
two middle indices proves
$(X\otimes Y)(Z\otimes W)=XZ\otimes YW$ whenever the dimensions match.

The implementation uses reversal of the second matrix index, exactly as
in {eq}`eq-construction-gs`:

```{literalinclude} ../../experiments/scalar_lift/williamson_scalar_endpoints.py
:language: python
:pyobject: gs
```

(construction-design)=
## Construct four disjoint coefficient matrices

Let $x_0,x_1,x_2,x_3$ be commuting real variables and set

```{math}
M(x)=\begin{pmatrix}
 x_0&-x_1&-x_2&-x_3\\
 x_1&x_0&x_3&-x_2\\
 x_2&-x_3&x_0&x_1\\
 x_3&x_2&-x_1&x_0
\end{pmatrix}.
```

Every column has squared length $\sum_jx_j^2$. The six mixed column
products are

```{math}
\begin{aligned}
&-x_0x_1+x_1x_0-x_2x_3+x_3x_2,&
&-x_0x_2+x_1x_3+x_2x_0-x_3x_1,\\
&-x_0x_3-x_1x_2+x_2x_1+x_3x_0,&
&x_1x_2+x_0x_3-x_3x_0-x_2x_1,\\
&x_1x_3-x_0x_2-x_3x_1+x_2x_0,&
&x_2x_3-x_3x_2+x_0x_1-x_1x_0.
\end{aligned}
```

They vanish term by term. Thus
$M(x)^{\mathsf T}M(x)=(\sum_jx_j^2)I_4$.
For the input function define $f^x(u)=M(x)b(u)$. Equation
{eq}`eq-construction-t` gives

```{math}
\sum_u f^x(u)^{\mathsf T}f^x(u+s)
 =\sum_u b(u)^{\mathsf T}M(x)^{\mathsf T}M(x)b(u+s)
 =p\Bigl(\sum_jx_j^2\Bigr)\delta_{s,0}.
```

By {eq}`eq-construction-gs-gram`,

```{math}
D(x):=\operatorname{GS}(f^x)=\sum_{j=0}^3x_jX_j,
\qquad D(x)D(x)^{\mathsf T}
 =p\Bigl(\sum_jx_j^2\Bigr)I_{4p},
```

where $X_j=\operatorname{GS}(f^{e_j})$ are integer matrices of order $4p$.
Substituting $x=e_j$, then $x=e_i+e_j$, proves

```{math}
:label: eq-construction-x
X_jX_j^{\mathsf T}=pI_{4p},\qquad
X_iX_j^{\mathsf T}+X_jX_i^{\mathsf T}=0\quad(i\ne j).
```

There is also an entry assertion. Since $b(u)$ is a signed basis vector,
$f^x(u)$ is a signed column of $M(x)$. Each of its entries is exactly one
signed variable. The block array only copies, negates and permutes these
entries. Hence at every matrix position exactly one $X_j$ is nonzero,
and that entry is $1$ or $-1$. The Gram identity alone would not establish
this partition of supports.

`coefficient_sequences` encodes the displayed matrix $M(x)$ through the
literal variable/sign table `L`. `verify_orthogonal_design` checks all
$10p$ diagonal and mixed polynomial correlation coefficients, then checks
the partition of supports of the four constructed matrices:

```{literalinclude} ../../experiments/scalar_lift/williamson_scalar_endpoints.py
:language: python
:pyobject: coefficient_sequences
```

```{literalinclude} ../../experiments/scalar_lift/williamson_scalar_endpoints.py
:language: python
:start-at:     checks = 0
:end-at:                 checks += 1
```

Checking these coefficients proves the complete matrix polynomial identity
through the already proved circulant and block-array calculation. The
certificate records this proof route separately from dense multiplication
of all coefficient matrices.

(construction-sixth-root)=
## Obtain the sixth-root frame and order 7572

Define the fixed integer matrix

```{math}
\Delta=\begin{pmatrix}
0&1&1&1\\-1&0&1&-1\\-1&-1&0&1\\-1&1&-1&0
\end{pmatrix},\qquad a(u)=\Delta b(u),\qquad d=4p.
```

Its transpose is $-\Delta$, its columns have squared length $3$, and their
mixed dot products are zero. Thus $\Delta^{\mathsf T}\Delta=3I_4$.
For every shift,

```{math}
\begin{aligned}
\sum_u a(u)^{\mathsf T}a(u+s)&=3p\delta_{s,0},\\
\sum_u\bigl(a(u)^{\mathsf T}b(u+s)+b(u)^{\mathsf T}a(u+s)\bigr)
 &=\sum_u b(u)^{\mathsf T}(\Delta^{\mathsf T}+\Delta)b(u+s)=0.
\end{aligned}
```

Put $U=\operatorname{GS}(a)$ and $V=\operatorname{GS}(b)$. Applying
{eq}`eq-construction-gs-gram` to $a$, $b$ and $a+b$ gives

```{math}
UU^{\mathsf T}=3pI_d,\qquad VV^{\mathsf T}=pI_d,\qquad
UV^{\mathsf T}+VU^{\mathsf T}=0.
```

The occupied coordinate of $b(u)$ is precisely the zero coordinate of
$a(u)$: it selects the same column in $I_4$ and $\Delta$. The other three
entries of that column of $\Delta$ are signs. Thus the corresponding
entries of $U,V$ have complementary supports and signed nonzero entries.
It follows that

```{math}
P=U-V,\qquad Q=U+V
```

are both sign matrices. Expanding their products proves

```{math}
:label: eq-construction-pq
PP^{\mathsf T}=QQ^{\mathsf T}=dI_d,\qquad
PQ^{\mathsf T}+QP^{\mathsf T}=dI_d.
```

These three identities produce a rational sixth-root action. Let
$\xi=\sqrt{-3}$, $\zeta=(1+\xi)/2$ and $K=\mathbb Q(\zeta)$. Since
$(2r-1)^2=-3$ has no rational solution, $t^2-t+1$ is irreducible over
$\mathbb Q$. Every element of $K$ therefore has a unique form $r+s\zeta$.
Set

```{math}
R=QP^{\mathsf T}/d.
```

From $P^{-1}=P^{\mathsf T}/d$ we also have $P^{\mathsf T}P=dI_d$.
Using {eq}`eq-construction-pq`,

```{math}
RR^{\mathsf T}=I_d,\qquad R+R^{\mathsf T}=I_d,\qquad RP=Q.
```

The first two imply $R^2-R+I_d=0$. Hence

```{math}
\rho(r+s\zeta)=rI_d+sR
```

preserves multiplication: both $\zeta$ and $R$ satisfy the same quadratic
relation. It is injective, since $rI_d+sR=0$ with $s\ne0$ would make a
rational scalar root of that irreducible polynomial. Finally,
$\bar\zeta=1-\zeta$ and $R^{\mathsf T}=I_d-R$ give
$\rho(\bar z)=\rho(z)^{\mathsf T}$.

The field degree is $[K:\mathbb Q]=2$. The action has rational dimension
$d=4p$ and module rank $d/2$ over $K$; these quantities have different
meanings. Its trace form in the actual rational coordinates is described
in {ref}`restriction-representation-dimension`.

Use the complete Butson source

```{math}
N=\begin{pmatrix}-1&\zeta&\zeta\\\zeta&-1&\zeta\\\zeta&\zeta&-1\end{pmatrix}.
```

Each row has squared norm $3$, and the product of distinct rows is
$1-\zeta-\bar\zeta=0$. Thus $NN^*=3I_3$. Its actual scalar entries have
the sign block images $\rho(-1)P=-P$ and $\rho(\zeta)P=Q$. Consequently

```{math}
:label: eq-construction-h7572
H=\begin{pmatrix}-P&Q&Q\\Q&-P&Q\\Q&Q&-P\end{pmatrix},
\qquad HH^{\mathsf T}=3dI_{3d}.
```

A diagonal Gram block is $PP^{\mathsf T}+2QQ^{\mathsf T}=3dI_d$.
An off-diagonal block is
$-PQ^{\mathsf T}-QP^{\mathsf T}+QQ^{\mathsf T}=0$.
This also verifies the final identity directly. Every entry of $H$ is a
sign, so its order and Gram scale both equal $3d=12p$. For the retained
$p=631$ witness this is $7572$.

The certified palette extends to $\{1,\zeta,-1,-\zeta\}$, with images
$P,Q,-P,-Q$. It does not include every sixth root: since
$\zeta^2=\zeta-1$,

```{math}
\rho(\zeta^2)P=Q-P=2V.
```

Each row of $V$ has $p$ signed nonzero entries among $4p$ positions.
Thus this block has zeros and entries of absolute value two. The actual
palette is part of the construction, as explained in
{ref}`restriction-flatness`.

`hadamard_frame_search.verify_endpoint` checks the source in exact
$\mathbb Z[\zeta]/(\zeta^2-\zeta+1)$ arithmetic, every coefficient shift,
the full $P,Q$ and cross identities, and every row pair of the assembled
endpoint. Its binary row words are exact encodings of signs: a row product
is $d-2\operatorname{Ham}(r,s)$, where $\operatorname{Ham}$ counts differing
bits.

```{literalinclude} ../../experiments/scalar_lift/hadamard_frame_search.py
:language: python
:start-at:     P,Q = gs(a-b),gs(a+b)
:end-at:     gram(rows,3*d)
```

The independently replayable order-7572 bundle is
{repo-file}`artifacts/hadamard-7572 <artifacts/hadamard-7572/README.md>`.

(construction-group-algebra)=
## Reuse the coefficient matrices with a group-algebra source

Let $A_0,A_1,A_2,A_3$ be sign matrices of order $w$ with

```{math}
\sum_jA_jA_j^{\mathsf T}=4wI_w,\qquad
A_iA_j^{\mathsf T}=A_jA_i^{\mathsf T}\quad(i\ne j).
```

The second condition is called pairwise amicability here. Set

```{math}
H=\sum_{j=0}^3X_j\otimes A_j.
```

At a position $((i,a),(k,b))$ exactly one $X_j(i,k)$ is nonzero. The entry
of $H$ is its product with the sign $A_j(a,b)$, proving the entry alphabet.
For $i<j$, the two mixed Gram terms are

```{math}
\begin{aligned}
&X_iX_j^{\mathsf T}\otimes A_iA_j^{\mathsf T}
 +X_jX_i^{\mathsf T}\otimes A_jA_i^{\mathsf T}\\
&\qquad=(X_iX_j^{\mathsf T}+X_jX_i^{\mathsf T})
                \otimes A_iA_j^{\mathsf T}=0.
\end{aligned}
```

The diagonal terms sum to
$pI_{4p}\otimes\sum_j A_jA_j^{\mathsf T}=4pwI_{4pw}$.
Thus $H$ is a real Hadamard matrix of order $4pw$.

For the retained circulant sources this is a literal scalar restriction.
Let $G$ now be a finite abelian group of order $w$. Define its rational
algebra with basis $c_g$, product $c_gc_h=c_{g+h}$, involution
$c_g^*=c_{-g}$, and trace functional
$\tau(\sum_g a_gc_g)=a_0$. In the inverse-group basis $\beta_g=c_{-g}$,
multiplication by $a=\sum_g a_gc_g$ has matrix

```{math}
:label: eq-construction-group-action
\rho(a)_{u,v}=a_{v-u}.
```

Indeed, $a\beta_v=\sum_g a_gc_{g-v}$ has coefficient $a_{v-u}$ at
$\beta_u$. Composition of multiplication maps gives
$\rho(ab)=\rho(a)\rho(b)$; inversion of the group index gives
$\rho(a^*)=\rho(a)^{\mathsf T}$. The action is faithful because applying
multiplication by $a$ to $c_0$ recovers $a$. Its matrix trace is $wa_0$, so
$\tau$ is the normalized regular trace. Also

```{math}
\tau(\beta_u^*\beta_v)=\delta_{u,v}.
```

The trace metric is therefore exactly $I_w$ in this declared basis.
Take $f_j=\sum_g a_j(g)c_g$ with all $a_j(g)$ signs and with

```{math}
f_if_j^*=f_jf_i^*,\qquad \sum_jf_jf_j^*=4w c_0.
```

Then $S=\sum_jX_jf_j$ has order $4p$ over $\mathbb Q[G]$. The same paired
term calculation proves $SS^*=4pwI_{4p}$. Its entries are signed copies
of the four $f_j$. Restriction sends them to signed copies of
$A_j=\rho(f_j)$ and gives precisely the displayed matrix $H$.

The source order is $4p$, its Gram scale is $4pw$, and the scalar algebra
has rational dimension $w$. Restriction keeps that Gram scale, which
equals the expanded order $4pw$. The algebra can be a product of fields;
its factor degrees and multiplicities are explained in
{ref}`restriction-products`.

`verify_williamson` checks the retained sign rows, self-adjointness,
amicability, the summed integer Gram, and an independent cyclic
convolution identity. `ScalarEndpoint.entry` is the literal final recipe:

```{literalinclude} ../../experiments/scalar_lift/williamson_scalar_endpoints.py
:language: python
:pyobject: ScalarEndpoint.entry
```

`ScalarEndpoint.row` reconstructs a full row without storing the complete
matrix. `ImplicitScalarEndpoint` uses the same entry proof when a source
quadruple itself has an implicit row recipe. These are representations of
the proved matrix, not approximations to it.

(construction-verification)=
## Replay the argument in implementation order

The construction has four finite gates:

1. `verify_t_rows` proves the alphabet, support partition and all input
   correlations.
2. `verify_orthogonal_design` proves the full coefficient identities and
   disjoint supports used in substitution.
3. `verified_williamson` or the corresponding retained source verifier
   proves its complete exact norm and paired identities.
4. The entry recipe and the Gram calculation above prove the whole final
   matrix. Additional dense audits are recorded separately when performed.

{repo-file}`The endpoint tests <tests/test_williamson_scalar_endpoints.py>`
compare all dense coefficient products with the polynomial proof for a
small length-seven control. They then construct every row and check every
row product for source widths $5,9,17$. They also reject corrupted input
and source data and check the field-factor dimensions.

```sh
PYTHONPATH=src:. .venv/bin/python -m pytest -q \
  tests/test_williamson_scalar_endpoints.py
```

The finite input is the only missing construction ingredient for a target
such as length $223$. Conditional on such an input, the same proof gives
orders $2676$, $4460$ and $11596$ from source widths $3$, $5$ and $13$.
These conditional formulas do not assert that those inputs have been found.
The retained manuscript and its Lean development contain the completed
existence results; the search records keep unfinished targets separate.
