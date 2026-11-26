# Williamson scalar endpoint order records

The complete periodic T631 design combines with the 27 literal Williamson
quadruples in [the input data](../experiments/scalar_lift/williamson_inputs.json)
to give completed sign-matrix recipes of order `2524*w`. The
[exact executor](../experiments/scalar_lift/williamson_scalar_endpoints.py)
checks the full orthogonal-design coefficient identity, every Williamson
Gram and amicability identity, the scalar source norm, and the actual sign
images. The weighted source is over `Q[C_w]`, a product of cyclotomic number
fields; it is not called a Butson matrix.

Order labels follow the [author-knowledge convention](HADAMARD_HIGHER_ORDER_RECORDS.md).
The [machine-readable status data](../experiments/scalar_lift/williamson_order_status.json)
give all 27 results and the prior construction parameters for each known
order. They were checked on 24 November 2026.

The [current source inventory](SOURCE_INVENTORY.md) extends these original
27 inputs with 19 compatible quadruples, twelve of which add new orders.
The new inputs can be symmetric back-circulants; their exact scalar algebra
is audited separately rather than identified with a cyclic group algebra.

## Complete source and scalar-transfer proof

Let `b=(b0,b1,b2,b3)` be the verified periodic T631 generating rows. Apply
the same Goethals–Seidel convention as in the standalone H7572 recipe to

\[
t_0=b,\quad t_1=(-b_1,b_0,b_3,-b_2),\quad
t_2=(-b_2,-b_3,b_0,b_1),\quad t_3=(-b_3,b_2,-b_1,b_0).
\]

The resulting four coefficient matrices `X_j` have order `4p=2524`,
entries `0,±1`, and disjoint supports covering every position. The full
periodic correlation identities and the commuting circulant GS identity give

\[
X_jX_j^T=pI,\qquad X_iX_j^T+X_jX_i^T=0\quad(i\ne j).
\]

The executor checks all ten polynomial Gram coefficients at every one of
the 631 shifts. Thus `D(x)=sum_j X_j*x_j` is a complete orthogonal design
with Gram `p*sum_j x_j^2*I` for commuting self-adjoint scalar variables.

For a literal symmetric Williamson quadruple of order `w`, write
`A_j[r,s]=a_j[(s-r) mod w]` and put

\[
K_w=\mathbb Q[c]/(c^w-1),\quad \overline c=c^{-1},\quad
f_j(c)=\sum_{k=0}^{w-1}a_j[k]c^k.
\]

The complete source checks give `f_j=conjugate(f_j)` and
`sum_j f_j^2=4w` in `K_w`. Consequently the **weighted source** `D(f)` has
algebraic order `4p`, actual alphabet the recorded distinct polynomials
`±f_j`, and Gram scale `4pw`. These entries need not be roots or have equal
absolute values in the component fields.

The algebra decomposes exactly as `K_w = product_(d|w) Q(zeta_d)`.
Use the transported positive trace form

\[
\tau(a)=[c^0]a=\frac1w\sum_{d\mid w}
\operatorname{Tr}_{\mathbb Q(\zeta_d)/\mathbb Q}(a_d),\qquad
\langle u,v\rangle=\tau(\overline u v).
\]

In the actual inverse group basis `beta_k=c^(-k)`, its metric is exactly
`I_w`, and multiplication satisfies
`rho(f)[r,s]=f_coefficients[(s-r) mod w]` and
`rho(conjugate(f))=rho(f)^T`. In the separate cyclotomic power bases the
metric is the block diagonal matrix with entries
`RamanujanSum_d(j-i)/w`; the independent audit verifies the exact change
of basis. The field degrees sum to `w`; `w` is not called a field degree.

Restriction therefore gives the literal endpoint

\[
H=\sum_{j=0}^3 X_j\otimes A_j,\qquad HH^T=4pwI.
\]

Every `X` position chooses exactly one signed `A_j` block, so **every
endpoint entry is a sign**. The source Gram scale remains `4pw` under
restriction; the rational order also becomes `4pw`, and the terminal
frame is `I_w` with Gram scale one. This proves all endpoint Gram entries.
Selected actual rows supplement that proof as interface checks.

The [retained certificate](../certificates/scalar-lift/williamson-scalar-endpoints.json)
contains all 27 complete source and endpoint proofs. The
[independent trace audit](../certificates/scalar-lift/williamson-cyclotomic-trace-audit.json)
checks exact cyclotomic actions, conjugation and metric transport and full
small endpoint Grams through order 588.

Replay with Python and NumPy:

```sh
python -m experiments.scalar_lift.williamson_scalar_endpoints \
  --output /tmp/williamson-scalar-endpoints.json
```

`endpoint(631,17)` returns the order-42908 object with exact `entry(i,j)`
and `row(i)` access. The script also runs outside the repository beside
its neighboring T631, Williamson-input and order-status JSON files.

## Additional new orders

Eight orders extend the previously completed results:

| Williamson order w | Real Hadamard order 2524*w |
|---:|---:|
| 17 | 42908 |
| 23 | 58052 |
| 29 | 73196 |
| 31 | 78244 |
| 37 | 93388 |
| 43 | 108532 |
| 51 | 128724 |
| 61 | 153964 |

The inputs `w=3,7,9,13,39` give new orders 7572, 17668, 22716, 32812 and
98436 already constructed by the sixth-root route. The Williamson recipes
are additional constructions of these orders; their equivalence classes
have not been compared.

## Positive prior constructions

The following 14 results have explicit prior constructions. Paley I and
II are [Paley's constructions](https://doi.org/10.1002/sapm1933121311),
also stated in the [2026 survey, Theorem 2](https://arxiv.org/html/2610.27175v2#S2.Thm2).
Here every displayed Paley parameter is prime.

| w | Endpoint order | Prior construction |
|---:|---:|---|
| 1 | 2524 | Explicit Goethals–Seidel construction in [Đoković–Golubitsky–Kotsireas, 2013](https://arxiv.org/abs/1301.3671v2) |
| 5 | 12620 | Paley I, q=12619 |
| 11 | 27764 | Paley I, q=27763 |
| 15 | 37860 | Paley sixteen-block product, q=12619, s=3155, t=3 |
| 19 | 47956 | Paley II, q=23977 |
| 21 | 53004 | Paley I, q=53003 |
| 25 | 63100 | Paley sixteen-block product, q=12619, s=3155, t=5 |
| 27 | 68148 | Paley I, q=68147 |
| 33 | 83292 | Paley sixteen-block product, q=27763, s=6941, t=3 |
| 41 | 103484 | Paley I, q=103483 |
| 49 | 123676 | Paley II, q=61837 |
| 55 | 138820 | Paley sixteen-block product, q=27763, s=6941, t=5 |
| 57 | 143868 | Paley II, q=71933 |
| 63 | 159012 | Paley sixteen-block product, q=53003, s=13251, t=3 |

The Paley block theorem gives `H(4*s*t)` from classical T-matrices of
order `t` and a prime power `q=4*s-1` congruent to 3 modulo 8. Its
hypotheses and the requisite small T-sequences are recorded in
[Tables 3 and 4 of the 2026 survey](https://arxiv.org/html/2610.27175v2#S4.T4).
Only `t=3,5` are needed above. In particular, 83292 is a known order,
although it is not reached by a direct Paley test at that order.

The arithmetic comparison also tested the Turyn Williamson family
`WT((q+1)/2)` with `q` a prime power congruent to 1 modulo 4, combined
with [Cooper–Wallis](https://doi.org/10.1017/S0004972700045019), and
Yamada's two order `4*(q+2)` constructions. For `w=31`, the latter
parameter is the prime 19559, congruent to 7 modulo 8, so neither of
those Yamada hypotheses applies. The JSON retains exact factorizations
of `N-1`, `N/2-1` and `N/4-2` for all 27 endpoints.

## The T127 control

`T127` with the order-17 Williamson source gives H8636, a known order.
The prior construction swaps the ingredients: classical T17 and known
Williamson-type matrices of order 127 give the same order by
Cooper–Wallis. [Đoković's 2009 paper, Section 3](https://arxiv.org/html/0912.5091#S3)
explicitly records `WT(127)` and includes `2159=17*127` among its
corrected known odd parts. This prior Williamson-type ingredient is
different from the new periodic T127 design.

The same argument makes all `4*127*w` with classical T(w), including
every one of the 27 input orders here, known real Hadamard orders.

## Input provenance

The literal data come from the
[official Sage table pinned at commit 8de8dad](https://github.com/sagemath/sage/blob/8de8dad7bcb38e1baa43508cbe944bf04a2d0401/src/sage/combinat/matrices/hadamard_matrix.py#L577).
The input JSON records the original source keys, exact first rows, row
sums, source commit and SHA256 hashes. Independent integer checks confirm
all entries are signs, all four circulants are symmetric, their summed
Gram is `4*w*I`, and every pair is amicable.

The dictionary item labelled 45 contains words of length 43 in that
pinned source. It is excluded and recorded as a rejected literal entry;
the 27 verified orders are used exactly as supplied.

## First expanded recurrence members

The [literal Spence recurrence](../experiments/scalar_lift/williamson_spence_family.py)
and its [complete source proofs](../certificates/scalar-lift/williamson-spence-family.json)
provide five families with source order `(q+1)*q^r/2`, for `r>=0` and
`q=13,17,61,101,121`. The source recurrence is
[Spence's Theorem 4 (1977)](https://doi.org/10.1017/S1446788700020255);
substitution into the verified T631 orthogonal design gives endpoint
order `1262*(q+1)*q^r`.

The first expanded members have separate
[finite order-status records](../experiments/scalar_lift/williamson_family_order_status.json):

| q | r | Source order | Real endpoint order | Status |
|---:|---:|---:|---:|---|
| 13 | 1 | 91 | 229684 | New |
| 17 | 1 | 153 | 386172 | New |
| 61 | 1 | 1891 | 4772884 | New |
| 101 | 1 | 5151 | 13001124 | New |
| 121 | 1 | 7381 | 18629644 | Known: Paley II, prime q=9314821 |

For the four new entries, the direct Paley and Yamada parameters are
not prime powers. The status data also check every factorization
`N=4*s*t` with `631` not dividing `t`: none passes either the Turyn
Williamson prime-power gate or the Paley sixteen-block prime-power gate.
These labels concern the five completed finite members; the recurrence
gives an infinite construction without assigning a novelty label to
every order in it.
