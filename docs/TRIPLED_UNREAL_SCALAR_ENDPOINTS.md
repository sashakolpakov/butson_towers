# Completed endpoints from the unreal tripling family

The known unreal sixth-root source of order 13 gives complete sources of
orders `13*3^k`. The new
[`tripled_unreal_scalar_endpoints`](../experiments/scalar_lift/tripled_unreal_scalar_endpoints.py)
executor verifies the literal sources of orders **39, 117 and 351** and reuses
the certified T631 frame to provide exact entry and row access to these
completed real endpoints:

Status labels follow the [author-knowledge convention and literature ledger](HADAMARD_HIGHER_ORDER_RECORDS.md).

| Source order and Gram scale | Frame scale | Endpoint order and Gram scale | Status |
|---:|---:|---:|---|
| 39 | 2524 | 98436 | New |
| 117 | 2524 | 295308 | New |
| 351 | 2524 | 885924 | Already known: Paley I at the prime 885923 |

These are complete literal matrix recipes with an exact proof of every
Gram entry. The executor stores bounded-index entry and row access and uses
the complete source and frame identities for the full Gram proof.
The retained [certificate](../certificates/scalar-lift/tripled-unreal-scalar-endpoints.json)
contains all exact source checks, gauges, frame bindings and interface controls.

The registered next source at order 1053 would give H2657772, which is
already known by Paley I with prime parameter 2657771. The
[current source inventory](SOURCE_INVENTORY.md) records larger compatible
sources and their positive prior-order checks.

## Source, gauge and family rule

The raw order-13 seed is the published
[M13A matrix](https://raw.githubusercontent.com/matrix-toolbox/chm/main/CHM/M13A.m).
Its [catalogue entry](https://raw.githubusercontent.com/matrix-toolbox/chm/main/catalogue/1306B.html)
credits R. Nicoară, A. LaClair, N. Geist and A. Wintenberg, May 2016.
The retained seed JSON records the literal row and column phase gauges
that make every entry an unreal sixth root.

Write an unreal source as

\[
U=(A+\sqrt{-3}B)/2,\qquad
AA^T+3BB^T=4nI,\qquad AB^T=BA^T,
\]

where `A,B` have sign entries. For

\[
C=\begin{pmatrix}0&1&-1\\-1&0&1\\1&-1&0\end{pmatrix},
\quad J=J_3,\quad
A'=B\otimes J,\quad B'=A\otimes I_3+B\otimes C,
\]

the identities `C^T=-C`, `CC^T=3I-J`, `CJ=JC=0` and `J^2=3J` give

\[
A'A'^T+3B'B'^T=3(AA^T+3BB^T)\otimes I_3=12nI,
\]

\[
A'B'^T-B'A'^T=(BA^T-AB^T)\otimes J=0.
\]

Every entry of `A'` and `B'` is a sign. Thus `(A'+sqrt(-3)B')/2` is an
unreal `BH(3n,6)`. This proves the family rule for every `k>=0`; the finite
executor separately checks all row and column Gram entries through order 351.
The known family is described in the publisher's abstract of
[Compton–Craigen–de Launey, *Unreal BH(n,6)'s and Hadamard matrices*,
Designs, Codes and Cryptography 79 (2016), 219–229](https://link.springer.com/article/10.1007/s10623-015-0045-y):
unreal orders `3^k*h*m`, with `h` a real Hadamard order and
`m` in `{1,7,10,13}`.

The existing phase-24 registry sources have literal exponent matrices
`E24=(4*S6+2) mod24`, where `S6` is the unreal sixth-root source. The
explicit global phase gauge `E24 -> E24-6` puts every exponent in
`{0,4,12,16}`. Division by four gives the actual sixth-root input
`E6=(S6-1) mod6`, with used alphabet exactly `{0,1,3,4}`. The executor
checks equality with each complete registry source and checks every row
and column Gram of the actual gauged input in the exact basis `(1,zeta_6)`.

## Scalar restriction and exact endpoint proof

The reused frame has rational dimension **2524**, field degree **2**, and
`Q(zeta_6)`-module rank **1262**. Its action is `R=Q*P^T/2524`, with
`R^T=I-R` and `R^2-R+I=0`; conjugation is transported by transpose.
The canonical trace metric is `[[2,1],[1,2]]`. In the actual coordinates,

\[
h(u,v)=u^Tv/2-(u^T(2R-I)v)\sqrt{-3}/6,
\qquad \operatorname{Tr}(h(u,v))=u^Tv.
\]

The complete `P,Q` Grams and symmetric cross identity are replayed from the
immutable H7572 recipe. Their packed hashes bind the exact frame. The
actual rational entry alphabet of `R` is retained in the new certificate.
Before applying the frame, the scalar blocks use `{0,1,-1}` and the
recorded rational `R` entries and their negatives. The actual used source
symbols then map as `0->P`, `1->Q`, `3->-P`, `4->-Q`, all literal sign blocks.

For a verified source `M` of order `n`, put `H_ij=rho(M_ij)P`. Exactly,

\[
(HH^T)_{ij}
=2524\,\rho\left(\sum_k M_{ik}\overline{M_{jk}}\right)
=2524n\,\delta_{ij}I.
\]

Hence the endpoint is dense with alphabet `{-1,1}`, order `2524n`, and
Gram scale `2524n`. This is a complete proof of all endpoint Gram entries,
using the complete checked source and frame identities. The general
family therefore has real orders `32812*3^k`, while only the three finite
entries in the table are newly executed here. None of these orders greater
than four can be a nontrivial Kronecker product of real Hadamard orders:
each is four times an odd number, whereas such a product would be divisible
by eight. This statement does not settle equivalence to another construction.

Replay the finite certificate:

```sh
PYTHONPATH=src:. .native-env/bin/python -m experiments.scalar_lift.tripled_unreal_scalar_endpoints \
  --output /tmp/tripled-unreal-scalar-endpoints.json
```

`endpoint(39)`, `endpoint(117)` and `endpoint(351)` return bounded-index
entry and row recipes. Selected full rows and their integer inner products
serve as implementation controls; the complete Gram proof is the exact
identity above. The [dated status ledger](HADAMARD_HIGHER_ORDER_RECORDS.md)
records the literature comparison for the two new orders and the prior
Paley construction at order 885924.

## Additional complete sources and the palette gate

The separate
[`additional_sixth_root_source_gates`](../experiments/scalar_lift/additional_sixth_root_source_gates.py)
executor records two concrete negative acquisition results in its
[certificate](../certificates/scalar-lift/additional-sixth-root-source-gates.json).
They concern compatibility with this fixed four-symbol frame, not source
existence or alternative scalar actions.

Szöllősi's thesis, [Theorem 1.4.41, printed page 30](https://arxiv.org/pdf/1110.5590),
gives `BH(p^2,6)` matrices `P tensor P+J tensor I+omega*I tensor J` from the
Paley character matrix, with `omega=zeta_6^2`. These have actual exponent
alphabet `{0,1,2,3}`. For every odd prime `p>=5`, their cubic-color matrix
contains four diagonal `2x2` blocks `[[1,2],[2,1]]`, with off-block zeros.
The resulting `8x8` pattern has no compatible gauge: fixing the first row
shift to zero, an independent enumeration checks all `3^7=2187` row gauges;
each leaves a column using all three cubic residues, so no column phase
can omit a class. The certificate verifies the literal embedding and full
row/column Grams for source orders 25, 49 and 121. The source of order nine
is retained as a positive gauge control. This does not exclude other
Butson matrices at those orders.

The literal `BH(19,6)` in
[Szöllősi's 2013 paper, page 33](https://ajc.maths.uq.edu.au/pdf/55/ajc_v55_p031.pdf)
also passes its complete exact row and column Grams, but its complete
sixth-root gauge test rejects the fixed representative. Other order-19
sources remain outside that finite conclusion. None of these rejected
inputs is counted as a completed scalar endpoint.
