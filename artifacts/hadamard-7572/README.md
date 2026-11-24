# Exact real Hadamard matrix of order 7572

This directory contains an explicit sign matrix, a complete sequence recipe,
and an exact verifier requiring only Python 3.10+ and NumPy. Relative to the
three primary records audited on **2026-11-24**, this settles a previously
unresolved real Hadamard order. The dated evidence and precise scope are in
[novelty.json](novelty.json); no absolute publication or method priority is
claimed.

Run from this directory:

```sh
python3 recipe.py
python3 recipe.py --write /tmp/h7572-replayed.npz
```

The verifier reconstructs the matrix, checks the source and every expanded
coefficient correlation, verifies both 2524-square frame Grams and their
symmetric cross identity, and checks all **28,663,806** distinct endpoint row
pairs using exact integer Hamming distances. It also compares the result
with the saved matrix and verifies its hash. It uses no search code or
repository imports. A fresh reconstructed output must have a new filename.
The `--write` replay also works with only `recipe.py` and `t631.json`, with
no matrix archive present. It always checks the fixed packed sign-bit hash;
when the reference archive is present, it additionally checks that archive.

`t631.json` gives the 31 orbit states and all four sequences of each of
`a,b`. State `2*axis+negative_bit` denotes a signed coordinate vector.
Starting with residue 0, then the least unused residue, multiply by
`mu=458` modulo 631; its order is 21. Each multiplication fixes coordinate
0 and cycles coordinates `1 -> 2 -> 3 -> 1`. This expands `b`; then
`a=D4*b`, with the explicit skew conference matrix in `recipe.py`.
The actual entries of `a,b` are `0,±1`, with disjoint supports and
`a_i(x)^2+b_i(x)^2=1`. Their row sums are respectively
`(21,-22,-22,-22)` and `(22,7,7,7)`.

Let `P=GS(a-b)`, `Q=GS(a+b)` using the literal Goethals–Seidel convention
in `recipe.py`, and put `d=2524`. Both matrices have entries `±1`, and the
exact checks establish

\[
PP^T=QQ^T=dI,\qquad PQ^T+QP^T=dI.
\]

The complete Butson source is the 3-square matrix with diagonal `-1` and
off-diagonal `zeta_6`, where `zeta_6^2-zeta_6+1=0`. Its actual alphabet is
`{zeta_6,-1}`, and its Gram scale is 3: off-diagonal row products are
`1-zeta_6-conjugate(zeta_6)=0`.

For restriction of scalars from `Q(zeta_6)` to `Q`, define the rational
action and frame by

\[
R=QP^T/d,\quad F=P,\quad RF=Q.
\]

The pair identities give `RR^T=I`, `R^T=I-R` (complex conjugation), and
`R^2-R+I=0`. Thus the rational dimension is 2524 and the field-module rank
is 1262. With `S=2R-I`, the explicit positive Hermitian form, linear in
its second argument, is

\[
h(u,v)=\tfrac12u^Tv-\tfrac16(u^TSv)\sqrt{-3},
\qquad \operatorname{Tr}_{\mathbb Q(\zeta_6)/\mathbb Q}h(u,v)=u^Tv.
\]

The field trace metric in the basis `(1,zeta_6)` is `[[2,1],[1,2]]`.
The rational blocks `-I,R` preserve the source Gram scale 3 before the
flat frame is applied. The actual rational entry alphabet of `R` is
recorded as exact numerators over 2524 in `certificate.json`.
The frame scale is 2524; the used source entries become `-P,Q`.
Consequently the endpoint is the literal block matrix

\[
H=\begin{pmatrix}-P&Q&Q\\Q&-P&Q\\Q&Q&-P\end{pmatrix},
\qquad HH^T=7572I,
\]

with actual alphabet `{-1,1}`. Source scale 3, frame scale 2524 and
endpoint order/scale 7572 are distinct quantities.

`h7572.npz` is a byte-for-byte copy of the originally verified artifact.
It contains `order=7572`, `bitorder="little"`, and the unsigned-byte array
`negative_entry_bits` of shape `(7572,947)`. Bit `j%8` of byte `j//8`
encodes a negative entry; unused high bits in each last byte are zero.
To obtain the dense matrix:

```python
import numpy as np
with np.load("h7572.npz", allow_pickle=False) as data:
    H = 1 - 2*np.unpackbits(data["negative_entry_bits"], axis=1,
                          count=int(data["order"]), bitorder="little").astype(np.int8)
```

NPZ SHA256: `cb0541a9019643d1e76698405c260731458627c01c8516bcb6d4d13c626d6364`.
Packed sign-bit SHA256:
`bc3cb6fca9243a6ef01c82de2a446080c8866e0f30d76c0f23ef0508ae5e2e9b`.
