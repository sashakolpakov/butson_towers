(known-orders)=
# Recovered known orders and matrix comparisons

The finite inventory contains **34 new orders and 31 recovered known orders**,
for **65 distinct completed orders**. The nine infinite families are counted
separately; this finite count does not enumerate their members.

A known order means that an earlier construction provides a Hadamard matrix
of that size. Whether our construction gives an equivalent matrix is a
further question. Two real Hadamard matrices are equivalent if independent
row and column permutations and sign changes transform one into the other.
Changing a recipe alone does not prove inequivalence.

## Positive prior constructions at all 31 orders

In the table, Paley I has order $q+1$ for $q\equiv3\pmod4$, and Paley II
has order $2(q+1)$ for $q\equiv1\pmod4$. All displayed $q$ are primes.
A Paley-block product has order $(q+1)t$, with $q\equiv3\pmod8$ and
T-matrices of length $t$. A Cooper–Wallis product has order $4tw$, with
T-matrices of length $t$ and a Williamson-type quadruple of order $w$.
The Turyn–Whiteman rows obtain the latter quadruple with $w=(q+1)/2$.
The table abbreviates Cooper–Wallis as CW and Turyn–Whiteman as TW.

The prior sources are [Paley](https://doi.org/10.1002/sapm1933121311),
[Seberry–Yamada, Theorem 4](https://documents.uow.edu.au/~jennie/WEBPDF/139_1990.pdf)
for the Paley-block product,
[Cooper–Wallis](https://doi.org/10.1017/S0004972700045019),
[Turyn](https://doi.org/10.1016/0097-3165(72)90095-7) and
[Whiteman](https://doi.org/10.1016/0097-3165(73)90010-1), and
[Đoković–Golubitsky–Kotsireas, §3.3](https://arxiv.org/abs/1301.3671v2)
for the cyclic supplementary difference sets.

| Order | Prior construction | Parameters |
|---:|---|---|
| 1524 | Paley I | $q=1523$ |
| 2524 | Cyclic supplementary difference sets | $v=631$ |
| 3252 | Paley I | $q=3251$ |
| 6564 | Paley I | $q=6563$ |
| 8636 | CW | $t=17, w=127$ |
| 12620 | Paley I | $q=12619$ |
| 27764 | Paley I | $q=27763$ |
| 37860 | Paley-block product | $q=12619, t=3$ |
| 47956 | Paley II | $q=23977$ |
| 53004 | Paley I | $q=53003$ |
| 63100 | Paley-block product | $q=12619, t=5$ |
| 68148 | Paley I | $q=68147$ |
| 83292 | Paley-block product | $q=27763, t=3$ |
| 88340 | Paley I | $q=88339$ |
| 103484 | Paley I | $q=103483$ |
| 113580 | Paley-block product | $q=12619, t=9$ |
| 123676 | Paley II | $q=61837$ |
| 138820 | Paley-block product | $q=27763, t=5$ |
| 143868 | Paley II | $q=71933$ |
| 159012 | Paley-block product | $q=53003, t=3$ |
| 164060 | Paley-block product | $q=12619, t=13$ |
| 239780 | Paley I | $q=239779$ |
| 335692 | CW/TW | $t=7, w=11989, q=23977$ |
| 456844 | Paley II | $q=228421$ |
| 885924 | Paley I | $q=885923$ |
| 3611844 | CW/TW | $t=9, w=100329, q=200657$ |
| 6564924 | CW/TW | $t=9, w=182359, q=364717$ |
| 11996572 | CW/TW | $t=7, w=428449, q=856897$ |
| 16257084 | CW/TW | $t=113, w=35967, q=71933$ |
| 18629644 | Paley II | $q=9314821$ |
| 23859372 | Paley I | $q=23859371$ |

The T-inputs of odd lengths through 19 in this table are classical.
For $t=113$, the published Turyn-type sequences of length 38 give base
sequences of lengths 75,75,38,38, and hence T-sequences of length 113; see
[Best et al., Corollary 5.1](https://arxiv.org/abs/1206.4107).
The order-8636 row uses previously known Williamson-type matrices of order
127, recorded by [Đoković, §3](https://arxiv.org/abs/0912.5091). That input differs from this project's new periodic T127 construction.
The order-2524 cyclic supplementary difference sets use a cyclic group of order $v=631$
and have block sizes 315,330,330,330 and
common difference multiplicity 674; Section 3.3 of the cited paper gives
the correct parameter.

The {repo-file}`provenance checker <scripts/audit_known_order_provenance.py>`
replays the integer hypotheses and binds the known orders to the positive
endpoint certificates. Its {repo-file}`record <certificates/scalar-lift/known-order-provenance.json>`
contains all reference links and selectors. These checks verify the
parameters of cited constructions; they do not reconstruct every earlier
published matrix or establish equivalence to it.

## A proved distinction at order 885924

Our sixth-root representative at order **885924** is inequivalent to the
standard Paley I matrix of that order. Their ranks over the field
$\mathbb F_3$ differ:

| Representative | Rank over $\mathbb F_3$ |
|---|---:|
| Our source of order 351 with the order-2524 frame | 246090 |
| Paley I, $q=885923$ | 442962 |

Signed permutation matrices are invertible modulo 3, so equivalent
Hadamard matrices must have the same rank over this field. The proof also
excludes equivalence to any skew Hadamard matrix of this order. It does not
classify all previously known classes at order 885924. Equivalence for the
other thirty recovered orders remains unsettled.

### Why the large matrix is unnecessary

The integral frame satisfies

```{math}
PP^{\mathsf T}=QQ^{\mathsf T}=dI,\qquad
PQ^{\mathsf T}+QP^{\mathsf T}=dI,\qquad d=2524.
```

Put $M=P+Q$. Then $MM^{\mathsf T}=3dI$. Its row space modulo 3 is
contained in its orthogonal complement, so its rank is at most $d/2$.
Taking determinants gives $v_3(\det M)=d/2$, where $v_3$ counts factors
of 3. If an integral matrix has rank $r$ modulo 3, elimination leaves
$d-r$ rows divisible by 3; therefore $d-r\leq v_3(\det M)$.
The two bounds force rank $1262$.

The sixth-root action is $R=QP^{-1}$, with $R^2-R+I=0$. Thus modulo 3,
$E=R+I$ satisfies $E^2=0$ and has rank 1262. Its image equals its kernel.
Choosing a basis of this image and preimages of its basis vectors splits
the frame space into 1262 identical two-dimensional summands. On each,
$R=-I+E$ and $E$ has matrix
$\left(\begin{smallmatrix}0&1\\0&0\end{smallmatrix}\right)$.

Let $A_3,B_3$ be the printed order-13 sign pair after three triplings.
Their order is 351. The source is
$\zeta_6^{-1}(A_3+\sqrt{-3}B_3)/2$. On the preceding two-dimensional
summand write $\zeta_6=-1+\epsilon$, $\epsilon^2=0$. Since
$\sqrt{-3}=2\zeta_6-1=2\epsilon$ modulo 3, its action has the same rank
as

```{math}
T=\begin{pmatrix}A_3&B_3\\0&A_3\end{pmatrix}
\quad\text{over }\mathbb F_3.
```

The invertible scalar factor and a basis rescaling remove the factor 2
in the upper-right block. Exact elimination gives rank 195. Consequently
the endpoint rank is $1262\cdot195=246090$: the actual endpoint is
$H^{(3)}=\rho(N_3)(I_{351}\otimes P)$, and its right factor is invertible
modulo 3.

The {repo-file}`comparison checker <experiments/scalar_lift/known_order_comparisons.py>`
reconstructs this calculation from the literal source and verifies both
rank bounds with a factorization certificate. If $U$ consists of the
recorded 195 columns of $T$, it checks $T=UV$ and $LU=I_{195}$ over
$\mathbb F_3$. The first identity bounds the rank above; the second,
together with the selected-column condition, bounds it below. The
factors are regenerated by the specified elimination algorithm; the selected
column indices and factor hashes are recorded. No matrix larger
than 702 by 702 is required.

For a skew Hadamard matrix $K$ of order $h$ divisible by 3,
$KK^{\mathsf T}=0$ over $\mathbb F_3$ gives rank at most $h/2$, while
$K+K^{\mathsf T}=2I$ gives $h\leq2\operatorname{rank}K$.
Its rank is therefore exactly $h/2$. Paley I has a skew representative,
so at $h=885924$ its rank is 442962.

This comparison has a mathematical proof and an exact executable finite
certificate. It is separate from the existing Lean existence proofs;
no Lean formalization of the finite-field rank comparison is claimed.
The full argument appears in the manuscript's section on recovered known
orders.

Replay both records from the repository root, after installing the project
dependencies:

```sh
python scripts/audit_known_order_provenance.py --check
python -m experiments.scalar_lift.known_order_comparisons \
  --check certificates/scalar-lift/known-order-paley-comparison.json
```

The first command needs only the Python standard library. The second uses
NumPy for exact arithmetic modulo 3 and checks a separate finite-field rank
calculation with FLINT when that optional package is installed. Its
{repo-file}`retained certificate <certificates/scalar-lift/known-order-paley-comparison.json>`
has the same content with or without FLINT.
