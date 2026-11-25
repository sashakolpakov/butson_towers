# Constructive infinite families from the T631 design

The complete T631 orthogonal design can be combined with Spence's recursive
Williamson quadruples. The new
[`williamson_spence_family`](../experiments/scalar_lift/williamson_spence_family.py)
executor gives literal entries for every integer `r>=0` in these families:

| Parameter q | Williamson component order | Real endpoint order and Gram scale |
|---:|---:|---:|
| 13 | `7*13^r` | `17668*13^r` |
| 17 | `9*17^r` | `22716*17^r` |
| 61 | `31*61^r` | `78244*61^r` |
| 101 | `51*101^r` | `128724*101^r` |
| 121 | `61*121^r` | `153964*121^r` |

These are the five Spence-source members of the project's six new output
families. The sixth is the [cubic family](TRIPLED_UNREAL_SCALAR_ENDPOINTS.md)
`32812*3^r`. Each output family has a new base rank; the completed endpoint
ledger records individually known ranks within the families.

The primary recurrence is Spence,
[*An infinite family of Williamson matrices* (1977), Theorem 4, printed page 254](https://doi.org/10.1017/S1446788700020255).
The executor checks all four recurrence formulas against the paper and
independently verifies its required starter identities. The source method
is known; the new frame supplies the endpoint multiplier 2524 from the
component dimension. Relative to the classical Hadamard rank `4*w_r`,
the output rank `2524*w_r` has the new odd multiplier 631.

## Exact starters and induction

Put `w0=(q+1)/2`. The starter quadruple consists of symmetric commuting
sign matrices `D,E,F,G` of order `w0`, satisfying

\[
M_0=(D-E)^2+(F-G)^2=4I,\qquad
P_0=(D+E)^2+(F+G)^2=4qI.
\]

These stronger identities are checked explicitly: an arbitrary Williamson
quadruple of the same order cannot simply be substituted. For `q=13,17`
the certificate supplies literal symmetric circulant conference halves
`R,S` with `R^2+S^2=qI`. Taking `D=E=S`, `F=I+R`, `G=-I+R` gives the
two identities. For `q=61,101,121`, explicitly recorded signs and channel
permutations of the retained complete Williamson inputs give the starters.

Let `Q` be the symmetric Paley character matrix of order `q` and `J=J_q`.
The executor checks `Q^2=qI-J` and `QJ=JQ=0` exactly. For `q=121` it uses
`GF(11)[u]/(u^2-2)` and indexes its additive group as `C11 x C11`.
Apply the same operation to each paired component `(A,B)=(D,E)` or `(F,G)`:

\[
A'=\frac{A-B}{2}\otimes J+\frac{A+B}{2}\otimes(Q+I),\qquad
B'=\frac{A-B}{2}\otimes J+\frac{A+B}{2}\otimes(Q-I).
\]

For any paired sign entries exactly one half-sum or half-difference is
nonzero. The zero diagonal of `Q` and the disjoint diagonal/off-diagonal
supports of `I,Q` therefore give only sign entries. Symmetry and
commutativity are preserved. Direct expansion gives

\[
M_{r+1}=P_r\otimes I,\qquad
P_{r+1}=qP_r\otimes I+(qM_r-P_r)\otimes J.
\]

Consequently `M_r=4q^r I`, `P_r=4q^(r+1) I`, and

\[
D_r^2+E_r^2+F_r^2+G_r^2
=\tfrac12(M_r+P_r)=4w_r I,\qquad w_r=w_0q^r.
\]

This proves every source Gram in the family. The recursive entry recipe
uses `O(r)` operations and never constructs a dense component matrix.

## Scalar algebra, metric and flatness

The components are group-circulant on
`G_r=C_w0 x (C_p)^(e*r)`, where `q=p^e`. Their sign first rows represent
self-conjugate commuting elements `a_0,...,a_3` of `Q[G_r]` with
`sum a_j^2=4w_r`. For the disjoint T631 design coefficients `X_j`,

\[
X_jX_j^T=631I,\qquad
X_iX_j^T+X_jX_i^T=0\quad(i\ne j),
\]

the complete weighted algebra-valued source `S=sum X_j a_j` has order
2524 over `Q[G_r]` and Gram `2524*w_r I`. Its actual entries are signed
copies of the four literal group polynomials; it is not asserted Butson.

The algebra is commutative semisimple. More explicitly, with
`m_r=(q^r-1)/(p-1)`,

\[
\mathbb Q[G_r]\cong
\prod_{d\mid w_0}\mathbb Q(\zeta_d)
\;\times\;
\prod_{d\mid w_0}\mathbb Q(\zeta_{pd})^{m_r}.
\]

Conjugation inverts group elements and is complex conjugation on each
field factor. The normalized regular trace is the identity coefficient,
or `1/w_r` times the sum of all field traces with those multiplicities.
It gives metric `I_w_r` in the inverse-group basis; each factor power
basis has trace-metric block `RamanujanSum_d(j-i)/w_r`.

Restriction in that group basis sends each `a_j` to its literal sign
group-circulant matrix `A_j`. The result

\[
H=\sum_{j=0}^3X_j\otimes A_j
\]

has order and Gram scale `2524*w_r`. Each position in the `X_j` is occupied
by exactly one sign, so every entry of `H` is a sign. Any zero field-factor
entries are preserved during the factor restrictions; the complete group
basis reconstruction is dense and flat.

## Finite replay

The retained
[certificate](../certificates/scalar-lift/williamson-spence-family.json)
independently executes the first recursive step for all five parameters.
It checks every group-correlation coefficient of the summed Gram and both
strong identities at source orders `91,153,1891,5151,7381`. The exact
integer NTT modulus exceeds the stated integer residual bound; there is no
floating-point rounding. Source orders 91 and 153 also receive complete
dense integer Gram checks. The resulting endpoint orders are
`229684,386172,4772884,13001124,18629644`; large dense endpoints are not
allocated by this source executor.

```sh
PYTHONPATH=src:. .native-env/bin/python -m experiments.scalar_lift.williamson_spence_family \
  --output /tmp/williamson-spence-family.json
PYTHONPATH=src:. .native-env/bin/python -m experiments.scalar_lift.williamson_spence_family \
  --q 61 --r 20 --entry 2 0 100
```

`SpenceSource(q,r)` provides `.order`, `.group_dimensions`,
`.entry(channel,i,j)`, `.entries(i,j)`, `.row(channel,i)`, and
`.coefficient_rows()`. The last two materialize the requested full rows;
individual entry access stays implicit at arbitrary `r`.
