# Complete sources and certified scalar endpoints

The [completed-order inventory](../certificates/scalar-lift/completed-order-inventory.json),
updated on 26 November 2026, contains **65 distinct completed real Hadamard
orders: 34 new to the best of the author's knowledge and 31 known recovered
orders**. Each entry points to a positive construction certificate and its
dated order-status record. The count concerns orders; it does not count
inequivalent matrices, incomplete searches or new ingredients as additional
real orders.

## Usable Butson seeds and their actual palette

The certified rational frame has width `d=2524`, with
`PP^T=QQ^T=dI` and `PQ^T+QP^T=dI`. Its sixth-root action sends

\[
1\mapsto P,\quad \zeta_6\mapsto Q,\quad
-1\mapsto-P,\quad-\zeta_6\mapsto-Q.
\]

Thus the admitted exponent palette modulo six is `{0,1,3,4}`. The other
two roots have images `Q-P` and `P-Q`, with nonflat entries. The chosen
complete seeds have orders **3, 7, 9 and 13**:

| Seed order | Retained source | Actual exponents after its recorded gauge |
|---:|---|---|
| 3 | `unreal3-q6` | 1,3 |
| 7 | `catalogue-p7-q6` | 0,1,3,4 |
| 9 | `aalto-9-6-000` | 0,1,3,4 |
| 13 | `unreal13-q6` | 0,1,3,4 |

These are complete selected sources with exact source Grams and literal
gauge witnesses, rather than a classification of their Butson classes.
Their certificates are [H7572](HADAMARD_7572.md) and
[frame reuse](../certificates/scalar-lift/hadamard-7572-frame-reuse.json).
The [unreal tripling recipe](TRIPLED_UNREAL_SCALAR_ENDPOINTS.md) gives the
additional completed order-39, -117 and -351 sources and an all-exponent
cubic output family.

## Larger complete Butson recipes

The compact `forward_source_registry.sources()` registry also retains 20
complete source recipes at other phases. Listing their descriptors does
not allocate their matrices; `load_source(id)` replays the appropriate
exact Gram proof before returning a verified source. The recorded phase
need not be the smallest possible phase after gauging.

| Recorded root order k | Matrix orders n in the registry |
|---:|---|
| 6 | 3, 13 |
| 8 | 1000, 2000 |
| 16 | 1000, 2000 |
| 20 | 5, 160 |
| 24 | 13, 39, 117, 351, 1053, 5329 |
| 40 | 40, 80 |
| 80 | 40 |
| 120 | 15, 35, 58081 |

In particular, there are large character-source recipes of orders 5329
and 58081, and order-1000/2000 sources over eighth or sixteenth roots.
These are genuine larger starting sources. They need a compatible scalar
action and final entry test for their actual palettes; the width-2524
sixth-root frame does not automatically admit them. The stored phase-24
cubic sources have an explicit gauge to the admitted sixth-root palette.
The separate higher-phase recipes remain available for other exact
reductions without being counted as additional real endpoints.

## Finite Williamson inputs

The original [literal input table](../experiments/scalar_lift/williamson_inputs.json)
contains **27 symmetric circulant Williamson quadruples**, at component
widths

`1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,37,39,41,43,49,51,55,57,61,63`.

For width `w`, the four sign matrices are self-adjoint commuting
multiplication images of elements of `Q[C_w]`, and their squared norms
sum to `4w`. The normalized regular trace gives metric `I_w` in the
inverse-group basis. Substitution into the complete T631 orthogonal
design gives a weighted algebra-valued source of **order 2524** and Gram
scale **2524w**. Restriction produces a real sign endpoint of order and
Gram scale `2524w`. The component width `w` and the weighted source order
2524 are distinct quantities. All 27 outputs have
[complete endpoint proofs](../certificates/scalar-lift/williamson-scalar-endpoints.json).

The [additional literal table](../experiments/scalar_lift/near_williamson_inputs.json)
supplies **19 further imported quadruples**, at widths

`35,45,47,53,59,65,67,69,73,79,91,93,95,97,103,109,127,133,181`.

Width 69 is an ordinary symmetric circulant Williamson quadruple from
the [Bright–Kotsireas–Ganesh data](https://doi.org/10.5281/zenodo.1195693).
The other 18 come from near-Williamson data: reversing the first circulant
gives four symmetric commuting sign matrices. Their scalar algebra is
`Q[A_0,A_1,A_2,A_3]` acting on `Q^w`, with identity conjugation and module
metric `I_w`. The actual multiplication matrices and source alphabet are retained in the
[endpoint certificate](../certificates/scalar-lift/near-williamson-scalar-endpoints.json).
The separate [algebra audit](../certificates/scalar-lift/near-williamson-scalar-algebras.json)
retains every factor, module multiplicity and reconstruction metric.
The same weighted-source order 2524 and Gram scale `2524w` apply. These
inputs give 12 additional new orders, six known orders and the already
counted order 229684; the
[finite status table](../experiments/scalar_lift/near_williamson_order_status.json)
records each decision.

## Explicit source recurrences and output families

There are now **11 verified Spence-source recurrences**, with prime-power
parameters
`5,13,17,53,61,97,101,109,113,121,137`.
The source method is Spence's
[1977 Theorem 4, p.254](https://doi.org/10.1017/S1446788700020255).
Its stronger starter identities, Paley matrices and sign-entry recurrence
are checked explicitly. For `q=p^e` the component source at step `r` has
order `(q+1)*q^r/2` and group
`C_((q+1)/2) x (C_p)^(e*r)`. Every output has an exact recursive entry
recipe; a large dense matrix is unnecessary.

The [original certificate](../certificates/scalar-lift/williamson-spence-family.json)
and [additional certificate](../certificates/scalar-lift/additional-spence-families.json)
keep the two implementation groups separate. The new module leaves the
original five parameter definitions and frozen certificates unchanged.

| Additional q | All-step real order, r>=0 | Completed r=1 order | Finite status of r=1 |
|---:|---|---:|---|
| 5 | `7572*5^r` | 37860 | Known |
| 53 | `68148*53^r` | 3611844 | Known |
| 97 | `123676*97^r` | 11996572 | Known |
| 109 | `138820*109^r` | 15131380 | New to the best of the author's knowledge |
| 113 | `143868*113^r` | 16257084 | Known |
| 137 | `174156*137^r` | 23859372 | Known |

The nine output families singled out as new are the original five
Spence-source families, the cubic family, and
`7572*5^r`, `174156*137^r`, **`15131380*109^r`**, for integers `r>=0`.
The last formula starts at source-recursion step one. Each contains an
individually new finite order; this does not label every order within it
new. Parameters 53,97,113 supply further fully constructive recurrences
and known initial-member controls, without an additional new-family claim.

Four higher members of the original source families also have complete
literal row and exact NTT norm certificates:

| q,r | Component order | Real order | Finite status |
|---|---:|---:|---|
| 13,2 | 1183 | 2985892 | New to the best of the author's knowledge |
| 13,3 | 15379 | 38816596 | New to the best of the author's knowledge |
| 61,2 | 115351 | 291145924 | New to the best of the author's knowledge |
| 17,2 | 2601 | 6564924 | Known |

The [additional finite status record](../experiments/scalar_lift/additional_spence_order_status.json)
separates these ten completed members from the universal construction
proofs. Endpoint controls use individual exact entries, while their
complete Grams follow from the verified polynomial orthogonal-design and
source identities.

## Distinct-order count and replay

The positive-certificate union is

\[
38+(19-1)+(6-1)+4=65.
\]

The first subtraction removes the repeated order 229684; the second
removes the repeated order 37860. The 34 new orders comprise the original
18, twelve additional finite Williamson orders and four additional finite
family members. The other 31 orders have positive earlier constructions.
The new periodic T127 ingredient is counted separately; its real endpoint
1524 is already included among the known orders.

The [construction and order inventory](../certificates/scalar-lift/completed-order-inventory.json)
lists every distinct order, its positive certificate selectors, input
hashes and status references. The executable replay for the additional
source families is:

```sh
PYTHONPATH=src:. .native-env/bin/python -m experiments.scalar_lift.additional_spence_families \
  --output /tmp/additional-spence-families.json
```

For individual entry access, use `source(q,r).entry(channel,i,j)` or
`endpoint(q,r).entry(i,j)` from that module. `finite_member(q,r)` performs
the full bounded-size NTT source replay; the recursive entry APIs also
work beyond that replay bound. Formal existence statements and their
axiom audits are recorded in [the formalization guide](FORMALIZATION.md).

## Known orders and matrix identity

The [separate prior-construction table](sphinx/known-orders.md) accounts for all
31 recovered known orders. Positive prior existence does not decide matrix
equivalence. Our order-885924 sixth-root representative has rank 246090 over
$\mathbb F_3$, whereas the standard Paley matrix has rank 442962, proving
inequivalence to that matrix (and to every skew representative). The other
thirty comparisons remain unsettled. This comparison has a separate exact
finite certificate; it does not alter the Lean existence package or the
34-new-order count.
