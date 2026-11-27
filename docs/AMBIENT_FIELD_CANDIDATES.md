# Actual scalar restrictions in new period ambient fields

The [executed report](../certificates/scalar-lift/ambient-field-candidates.json) contains **12 verified source/ambient-field recipes and 24 pure restriction paths**, reaching eight distinct rational orders: **4676, 6876, 7924, 8684, 9308, 9716, 10052, 11596**. Every declared finite recipe is finished. No right frame or normalization operation occurs in those paths.

**All product-basis endpoints have zero entries and therefore fail the ±1 condition.** This is a structural exclusion for the whole executed product-basis family, independent of the auxiliary conductor. Changing only the real period trace basis commutes with every source scalar and cannot remove those zeros. There is no further conductor search pending in this family. Nonproduct joint bases and unequal input/output frames are outside this computation; no new real Hadamard order is claimed.

The subsequent [mixed-basis audit](https://github.com/sashakolpakov/complex-hadamard/blob/976cede671f40c810d593c8d5815e26e11a03ed6/docs/MIXED_RESTRICTION_BASIS_AUDIT.md) proves a
stronger obstruction: every shared input/output basis change preserves a
determinant incompatible with the expanded Hadamard order. Mixing factors
alone cannot fix it. Unequal input and output frames lie outside that obstruction.

The retained sources are the exact order-13 sixth-root input, the separately acquired order-7 sixth- and twelfth-root inputs, and the existing literal order-9 sixth-root input. Their full Grams and source hashes are replayed. These endpoints occur in the attributed bounded frontier snapshot. None of its retained orders above 20,000 is divisible by any of the fixed source orders 3, 7, 13, 15, 35. The actual order-15 and order-35 common-base inputs generate Q(zeta120); replacing their field by the sextic field would be invalid.

| Real degree r | Prime conductor ell | h=(ell−1)/r | Norm witness (a,b,t) | Sources | Retained endpoints |
|---:|---:|---:|---|---|---|
| 167 | 294923 | 1766 | (3850,233,9) | 7 at phases6/12; 13 at phase6 | 4676; 8684 |
| 179 | 268501 | 1500 | (721,186,5) | 13 at phase6 | 9308 |
| 191 | 1431737 | 7496 | (27861,788,25) | 9 at phase6 | 6876 |
| 223 | 189997 | 852 | (5369,342,17) | 13 at phase6 | 11596 |
| 283 | 748819 | 2646 | (5384,165,7) | 7 at phases6/12 | 7924 |
| 347 | 3191707 | 9198 | (143224,3241,87) | 7 at phases6/12 | 9716 |
| 359 | 3856379 | 10742 | (17582,95,9) | 7 at phases6/12 | 10052 |

All norm witnesses satisfy `a²+r b²=ell t²`, with odd t; all conductors satisfy `2^h=1 mod ell`. The degree-167 conductor is new and differs from the earlier conductor157649 example.

## Literal restriction paths and exact replay

For E=L Q(zeta12), the two paths are

`E -> Q(zeta12) -> Q(zeta6) -> Q`, with relative degrees `(r,2,2)`,

and

`E -> L(zeta6) -> L -> Q`, with relative degrees `(2,2,r)`.

Their actual source orders are `(nr,2nr,4nr)` and `(2n,4n,4nr)`. Source Gram scale remains **n at every stage**, with identity

`R (I_n tensor G^-1) R* = n (I_n tensor G^-1)`.

The cyclotomic product basis is `{1,zeta6,zeta12,zeta12*zeta6}`, and the period basis is the eta_j below. The first final basis is ordered `(j,a,b)` and the second `(a,b,j)`. Every used scalar, the full trace metric, and the faithful auxiliary multiplication operator obey the certified perfect shuffle between these coordinates. The final nodes are kept distinct: permutation equivalence does not imply literal equality or basis-independent flatness.

The [pure restriction implementation](../experiments/scalar_lift/pure_period_restriction.py) stores only small scalar tables and period Kronecker factors. `PeriodRestrictionStage.entry(i,j)` returns exact coefficient vectors in the declared lower field. When L remains in the coefficient field, the returned cyclotomic coefficients embed as rational constants using `1=-sum eta_j`. No dense `n²` output is required.

`family(r, source_id, bundle=path)` in [ambient_field_candidates.py](../experiments/scalar_lift/ambient_field_candidates.py) saves the entire definition, including original source exponents, faithful period multiplication, trace metrics, scalar tables, and stage hashes. The [retained report](../certificates/scalar-lift/ambient-field-candidates.json) summarizes twelve executed bundles; generate a bundle at a fresh local path to replay its stored arrays and metadata with `SavedPeriodBundle(path)`. The method `.entry(i,j,path=0,stage=-1)` accesses a pure restriction stage. Its separately named `.supplementary_normalized_entry(i,j)` accesses only the optional normalized control described below.

Literal deduplication hashes the source, coefficient field, actual scalar action table, and actual metric in their stated coordinates. The report gives node counts and zero pending recipes. Only literal equalities are deduplicated.

## Faithful period action

Let g be the certified primitive root modulo ell, H=<g^r> the order-h subgroup, and

`eta_j = sum_{x in g^j H} zeta_ell^x`.

Because h is even, these periods are real. Their trace Gram is

`G = ell I_r - h J_r`.

It is positive definite: its eigenvalues are 1 on the constant vector and ell on its orthogonal complement. The r periods are therefore a rational basis for the fixed field L of H, of degree r. As r is prime and eta_0 is nonrational, eta_0 generates L. The compositum E=Q(zeta12)L has degree 4r because the two field degrees are coprime.

The actual multiplication matrix is retained, not replaced by an identity auxiliary action. In column coordinates,

`M[k,j] = #{x in H: g^k-x belongs to g^j H} - h*[j=0]`.

The subtracted term replaces the constant contribution `h*[j=0]` using `1=-sum eta_k`. Computing all entries uses only r h=ell-1 membership lookups. Exact checks verify `M^T G=G M`, trace M=-1, and commutation with its Galois conjugate. An independent small-field test checks this formula directly by full cyclotomic convolution.

The congruence `2^h=1 mod ell` means 2 splits completely in L. This is a local compatibility condition, not a flat-frame existence theorem.

## Exact trace normalization

Let Q be the skew Jacobsthal matrix of order r, so `Q^T=-Q`, `Q Q^T=rI-J`, and `QJ=0`. With the displayed norm witness, put

`S = [r a I - r b Q + (t-a+h r t) J] / (r t ell)`.

Then `S^T G S=I`. On the constant vector S is the identity; on its orthogonal complement it is `(aI-bQ)/(t ell)`, where the norm identity proves the assertion. All full integer matrix identities are checked using arbitrary-precision FLINT matrices. The denominator is odd, and the numerator reduces to I modulo two, so S preserves the period lattice locally at two.

The degree-716 full cyclotomic fields containing the sextic field have conductors 1077 and 2154; both fail the existing standard-metric dyadic test. There is no full cyclotomic field of degree 892 containing that field. The proper subfields above are distinct options, not alternate names for those rejected full fields.

## Supplementary normalization control: not a restriction stage

Use one of the verified complete sources of order n, with its exact source hash recorded. On rational quaternions let D be left multiplication by `i+j+k`, and J right multiplication by i. Then `D²=-3I`, `J²=-I`, and JD=DJ. The sextic action is `(I+D)/2`, and `J(I-D)/2` is a faithful primitive twelfth-root action. Tensoring these operators with the period action gives a faithful representation of E, with metric `I_4 tensor G`; the basis `I_4 tensor S` makes that metric standard.

The source entries lie in Q(zeta12), and the sextic inputs lie in its quadratic subfield. Thus their numerical operators are the four-dimensional local action tensored with `I_r`. **Changing the auxiliary field alone has not changed this repeated source matrix.** A rational quaternion normalization frame P with `PP^T=4r I_4` gives an actual integer matrix

`N = N0 tensor I_r`, with `N0 N0^T = 4nr I_(4n)`.

Every full entry is `N[i,j]=N0[i//r,j//r]*[i%r=j%r]`. The bundles contain N0 and the complete ambient-field generators; no huge dense output is needed. The nonflat integer alphabet is recorded explicitly. The omitted construction is a nonproduct joint frame whose images of the used source entries are all signs. No such frame is executed or asserted here. Neither trace normalization nor the displayed norm frame establishes that condition.
