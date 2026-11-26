# From a finite T-matrix witness to a Butson-tower endpoint

This note explains the construction proved in the manuscript. The search
finds its finite input. Once that input is verified, the subsequent matrix
identities are exact construction proofs. The submitted manuscript and
its frozen packages are preserved separately.

## 1. The required endpoint

A real Hadamard matrix of order $N$ is a sign matrix $H$ satisfying
$HH^{\mathsf T}=NI_N$. Two facts must be established: every entry is a
sign, and the Gram scale is exactly its order.

A Butson matrix of order $n$ has roots of unity as entries and satisfies
$BB^*=nI_n$. Restriction of scalars transports its norm identity with the
appropriate trace metric. The entry alphabet, trace metric, representation
dimension and Gram scale must all be retained. Flatness of the final
entries is a further condition.

## 2. The finite input: periodic T-matrices

Let $p$ be an odd prime, and let $e_0,e_1,e_2,e_3$ be the standard basis of
$\mathbb Z^4$. We seek
\[
b:\mathbb Z/p\mathbb Z\longrightarrow\{\pm e_0,\pm e_1,\pm e_2,\pm e_3\}
\]
with
\[
C(s):=\sum_x b(x)^{\mathsf T}b(x+s)=p\delta_{s,0}.
\]
The coordinate functions $b_j$ take values in $\{0,1,-1\}$, and at each
position exactly one is nonzero. Define the circulant matrices
$T_j(u,v)=b_j(v-u)$, with all indices reduced modulo $p$. Then
\[
\left(\sum_jT_jT_j^{\mathsf T}\right)_{u,v}=C(v-u),
\qquad \sum_jT_jT_j^{\mathsf T}=pI_p.
\]
Their supports are disjoint and cover every matrix position. These are the
periodic T-matrices used in the paper.

For the successful $p=631$ search, multiplication by $458$ has order $21$.
We impose $b(458x)=\Pi b(x)$, where $\Pi$ fixes $e_0$ and cycles
$e_1,e_2,e_3$. One signed basis vector per multiplication orbit determines
$b$; there are 31 such orbits including the residue zero. The search uses
15 nonzero shift representatives and verifies all 631 shifts at the end.

## 3. Finding the input by exact modular lifting

First choose the occupied coordinate at each position or orbit. Modulo two,
its sign is irrelevant. The unsigned correlation parities are quadratic in
the binary encodings of these four coordinate labels.

For the successful symmetric chart at 631, an explicitly verified invertible
binary coordinate change puts all parity conditions into the single equation
\[
AD+BC=1,\qquad A,B,C,D\in\mathbb F_{2^{15}}.
\]
Every coefficient of this quadratic identity is checked. We sample this
locus directly: for $A\ne0$, choose $A,B,C$ and set
$D=(1+BC)/A$; for $A=0$, choose $B\ne0$, set $C=B^{-1}$, and choose $D$.
Thus no time is spent sign-lifting coordinate patterns that already fail
parity. The symmetric 223 chart uses $\mathbb F_{2^{37}}$. The unrestricted
167 and 179 charts use $\mathbb F_{2^{83}}$ and $\mathbb F_{2^{89}}$.
Their maps have their own complete coefficient and expanded-shift checks.

For a fixed coordinate pattern, write each sign as $1-2x_i$, with
$x_i\in\{0,1\}$. Each correlation can be written
\[
C_\ell=C_\ell^{\rm abs}-2\sum_{i<j}w_{ij,\ell}(x_i\mathbin{\oplus}x_j).
\]
The nonnegative integer $w_{ij,\ell}$ counts the relevant equal-coordinate
pairs. Modulo four, the equations are linear over $\mathbb F_2$. Binary
elimination gives $x=x^{(0)}+By$. Modulo eight, the remaining conditions
are quadratic Boolean equations in the free bits $y$.

We also impose the exact coordinate sums. If $r=\sum_xb(x)$, then summing
all correlations gives $\|r\|^2=p$. For a fixed coordinate support of size
$n_j$, the number of negative signs is exactly $(n_j-r_j)/2$. All four
counts are imposed in the unrestricted search.

The original 631 chart typically left 13 free bits, so complete native
enumeration was cheap. The current 223 chart leaves roughly 32--35 bits.
The unrestricted 167 and 179 models typically leave 80 and 86 bits;
SAT and CP-SAT represent the exact constraints without requiring an initial
enumeration of all $2^{80}$ or $2^{86}$ assignments. Their search can still
be expensive. An accepted witness is reconstructed independently and checked
at every integer shift. A timeout is an unfinished fixed-pattern search.

## 4. From the finite input to a disjoint orthogonal design

For real variables $x_0,x_1,x_2,x_3$, put
\[
M(x)=\begin{pmatrix}
x_0&-x_1&-x_2&-x_3\\
x_1&x_0&x_3&-x_2\\
x_2&-x_3&x_0&x_1\\
x_3&x_2&-x_1&x_0
\end{pmatrix}.
\]
Direct multiplication gives $M(x)^{\mathsf T}M(x)=(\sum_jx_j^2)I_4$.
Consequently the four functions forming $f^x(u)=M(x)b(u)$ satisfy
\[
\sum_u f^x(u)^{\mathsf T}f^x(u+s)
=p\left(\sum_jx_j^2\right)\delta_{s,0}.
\]

For the circulants $C_0,C_1,C_2,C_3$ of four functions and the reversal
$J_{u,v}=\delta_{u+v,p-1}$, the manuscript uses
\[
\operatorname{GS}(C_0,C_1,C_2,C_3)=\begin{pmatrix}
C_0&C_1J&C_2J&C_3J\\
-C_1J&C_0&-C_3^{\mathsf T}J&C_2^{\mathsf T}J\\
-C_2J&C_3^{\mathsf T}J&C_0&-C_1^{\mathsf T}J\\
-C_3J&-C_2^{\mathsf T}J&C_1^{\mathsf T}J&C_0
\end{pmatrix}.
\]
Circulants and their transposes commute, and $JC_jJ=C_j^{\mathsf T}$.
These identities cancel every off-diagonal Gram block; each diagonal
block is $\sum_j C_jC_j^{\mathsf T}$.
Applying this array to $f^x$ gives
\[
D(x)=\sum_{j=0}^3x_jX_j,\qquad
D(x)D(x)^{\mathsf T}=p\left(\sum_jx_j^2\right)I_{4p}.
\]
In particular,
\[
X_jX_j^{\mathsf T}=pI_{4p},\qquad
X_iX_j^{\mathsf T}+X_jX_i^{\mathsf T}=0\quad(i\ne j).
\]
Since $b(u)$ is a signed basis vector, $M(x)b(u)$ is a signed column of
$M(x)$. Every entry is exactly one signed variable. This remains true
through the block array. Hence exactly one $X_j$ is nonzero at each
position, and that entry is $1$ or $-1$.

## 5. The original order-7572 construction

Define
\[
\Delta=\begin{pmatrix}0&1&1&1\\-1&0&1&-1\\-1&-1&0&1\\-1&1&-1&0\end{pmatrix},
\qquad a=\Delta b.
\]
Here $\Delta^{\mathsf T}=-\Delta$ and $\Delta^{\mathsf T}\Delta=3I$.
Because $b(u)$ has one signed nonzero, both $a(u)-b(u)$ and $a(u)+b(u)$
have sign entries. Moreover
\[
(\Delta\pm I)^{\mathsf T}(\Delta\pm I)=4I,\qquad
(\Delta-I)^{\mathsf T}(\Delta+I)+(\Delta+I)^{\mathsf T}(\Delta-I)=4I.
\]
These identities, the correlation identity for $b$, and polarization of the
linear block array give, with $d=4p$,
\[
P=\operatorname{GS}(a-b),\quad Q=\operatorname{GS}(a+b),\quad
PP^{\mathsf T}=QQ^{\mathsf T}=dI,\quad
PQ^{\mathsf T}+QP^{\mathsf T}=dI.
\]
Set $R=QP^{\mathsf T}/d$. Then
\[
RR^{\mathsf T}=I,\qquad R+R^{\mathsf T}=I,\qquad R^2-R+I=0,
\qquad RP=Q.
\]
Thus $R$ supplies a rational action of $\mathbb Q(\zeta_6)$, with complex
conjugation represented by transpose. The field degree is two; $d$ is the
representation dimension.

Take the verified Butson source
\[
B=\begin{pmatrix}-1&\zeta_6&\zeta_6\\\zeta_6&-1&\zeta_6\\\zeta_6&\zeta_6&-1\end{pmatrix},
\qquad BB^*=3I_3.
\]
Its off-diagonal Gram entry is $1-\zeta_6-\overline{\zeta_6}=0$.
Replacing its scalars by the rational action gives blocks $-I_d,R$ and
preserves Gram scale 3. Multiplication by $I_3\otimes P$ then gives
\[
H=\begin{pmatrix}-P&Q&Q\\Q&-P&Q\\Q&Q&-P\end{pmatrix},
\qquad HH^{\mathsf T}=3dI_{3d}.
\]
All entries are signs because the actual used blocks are $-P,Q$.
The frame contributes the factor $d$ in the Gram scale. For $p=631$,
$d=2524$ and $3d=7572$. For $p=223$, the identical proof would give 2676.
There is no assertion that every sixth-root palette passes this entry test.

## 6. Reuse with other verified sources: the group-algebra restriction

Let $G$ be a finite abelian group of order $w$. In $\mathbb Q[G]$, write
$c_gc_h=c_{g+h}$, $c_g^*=c_{-g}$ and $\tau(\sum_g a_gc_g)=a_0$.
In the inverse-group basis $\beta_g=c_{-g}$, multiplication has matrix
\[
\rho(a)_{u,v}=a_{v-u},\quad
\rho(ab)=\rho(a)\rho(b),\quad
\rho(a^*)=\rho(a)^{\mathsf T},\quad
\tau(\beta_u^*\beta_v)=\delta_{u,v}.
\]
Thus the declared trace metric is exactly $I_w$.

Our known sources supply elements $f_0,\ldots,f_3$ with sign coefficients,
\[
f_if_j^*=f_jf_i^*,\qquad \sum_jf_jf_j^*=4w\,c_0.
\]
The matrix $S=\sum_jX_jf_j$ has order $4p$ over this algebra and satisfies
$SS^*=4pwI_{4p}$. Cross terms cancel by the displayed paired identities;
the diagonal terms sum to the stated norm. Its actual entries are signed
copies of the four $f_j$.

Restricting scalars replaces each $f_j$ by $A_j=\rho(f_j)$:
\[
H=\sum_jX_j\otimes A_j,\qquad HH^{\mathsf T}=4pwI_{4pw}.
\]
The $A_j$ have sign entries. Disjointness of the $X_j$ ensures exactly one
summand contributes to each final entry. This proves flatness separately
from the transported Gram identity.

Here the source order is $4p$, its Gram scale is $4pw$, and the scalar
algebra has dimension $w$. Restriction preserves that Gram scale, which
equals the expanded order $4pw$. For $G=C_w$ the algebra decomposes as
$\prod_{\ell\mid w}\mathbb Q(\zeta_\ell)$. All factors, traces and
reconstruction maps are retained; individual factor coordinates need not
be flat.

These steps explain the target orders:

| T-matrix length $p$ | Verified source width $w$ | Resulting order $4pw$ |
|---:|---:|---:|
| 223 | 3 | 2676 |
| 223 | 5 | 4460 |
| 167 | 7 | 4676 |
| 179 | 11 | 7876 |
| 167 | 13 | 8684 |
| 179 | 13 | 9308 |
| 223 | 13 | 11596 |

Each row is conditional on its missing T-matrix input. The small source
quadruples are already verified.

## 7. Towers, families and certification

A Butson tower is our name for successive scalar restrictions. Norms are
carried through each step with their forms, involutions and scalar actions.
The final entry test applies after the complete reconstruction, so
intermediate matrices and individual field factors may be nonflat.

The manuscript also uses known source recurrences. Reusing the same T631
coefficient matrices with source widths $w_r$ gives order $2524w_r$.
For example, $w_r=7\cdot13^r$ gives $17668\cdot13^r$. The other Spence
families and the sixth-root tripling family have their stated exact
recurrences and all-parameter proofs in the manuscript and Lean development.

For a new finite input we require: a literal signed-coordinate witness;
all integer correlations; the coefficient-matrix identities; the exact
source identities and used alphabet; the resulting order and Gram scale;
and the final entry test. Only then is there a completed real matrix.
Literature novelty and a new Lean literal-witness proof are subsequent
checks. A fixed-pattern UNSAT result or a timeout is not a nonexistence
result for a target order.
