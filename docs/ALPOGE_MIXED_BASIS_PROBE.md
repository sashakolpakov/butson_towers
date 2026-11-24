# Mixed two-basis restrictions using published Hadamard inputs

The [executed probe](../certificates/scalar-lift/alpoge-mixed-basis-probe.json)
uses the published order-668 width-12 matrix and order-716 width-4 matrix
as fixed input-basis frames. Their full sign Grams are independently
verified. Their construction is replayed from published seeds; no new source
search or constraint solver is used.

This follows the correction in the
[basis audit](https://github.com/sashakolpakov/complex-hadamard/blob/976cede671f40c810d593c8d5815e26e11a03ed6/docs/MIXED_RESTRICTION_BASIS_AUDIT.md): **a shared change of basis
cannot turn a pure scalar restriction into an expanded-order Hadamard**.
Similarity preserves the wrong determinant. Different input and output bases
are legitimate coordinates for the restricted linear map and avoid this
invariant.

For source entries in `{1,zeta6,-1,-zeta6}`, choose an output basis with
orthogonal scalar action `R=Q R0 Q^T` and input basis equal to that output
basis times the published matrix F. The restricted scalar blocks are

`B(1)=F`, `B(zeta6)=G=Q R0 Q^T F`.

Here `R0=(I+L_(i+j+k))/2 tensor I_r`, `d=4r`. The required conditions are

`F F^T = G G^T = d I`, `G F^T+F G^T=d I`, and `F,G` have only sign entries.

The Gram and cross identities hold for every rational orthogonal Q. The
remaining literal sign condition is tested exactly. These blocks define the
complete restriction of the verified order-7 and order-13 sextic inputs;
`restricted_entry` reads any entry without allocating the full target.
Targets are **4676, 8684 and 9308**.

## Actual ambient fields and metric

The fields are `Q(zeta12) L_167` with conductor294923 and
`Q(zeta12) L_179` with conductor268501. Their faithful period generators,
trace standardizers, and hashes are inherited from the exact period-field
certificate. The auxiliary operator is
`I4 tensor (S^-1 M_eta S)`, conjugated by the same Q as all local operators.
Thus the field has not been replaced by an unspecified multiplicity space.

The standard metric here is an explicit **weighted** trace form, not the
unweighted canonical trace. With `D=L_(i+j+k)`, `J=right_i`, use
`theta=(3-DJ)/12`; both real embeddings are positive. The basis

`1, J, (1+D-J+DJ)/2, (-1+D-J-DJ)/2`

has weighted trace Gram I4. All sixteen entries are checked exactly. Tensor
this with the certified normalized period basis. The new output basis has
weighted metric I; its input basis has weighted metric `F^T F=d I`.

## Completed reflection batch and exact pruning rule

For each input the batch evaluates the unchanged output basis and 64 dense
rational Householder changes `Q=I-2vv^T/(v^T v)`. It finishes **130 candidate
visits, 130 distinct partner matrices, and no sign partner**. Each candidate
stores its actual vector, exact partner hash, entry-class census, failed
entry, and a tensor-minor witness when its basis column is nonproduct.
Evaluation uses rank-one identities, with no dense cubic multiplication.

The result is stronger than failure of those samples. Put `A0=R0 F`, which
is integral. If Q is orthogonal and `rank(Q-I)<=k`, then

`rank_Q(G-A0)<=2k`.

If G were a sign matrix, its reduction modulo two would be the all-ones
matrix J. For the integral matrix `G-A0`, reduction cannot increase rank:

`rank_F2(A0-J)<=rank_Q(G-A0)<=2k`.

Exact binary elimination gives:

| Published F | rank over F2 of A0-J | Necessary rank(Q-I) |
|---:|---:|---:|
| 668 | 167 | at least84 |
| 716 | 178 | at least89 |

Consequently **all single-reflection changes, and every product of fewer
than84/89 Householder reflections respectively, are excluded**, regardless
of how dense their vectors are. Future basis probes must change sufficiently
many independent directions; another shallow reflection tree cannot work.
This does not exclude general two-basis restrictions or other scalar-action
classes.

The code includes a successful known-order positive control: the width-4
two-basis restriction of the verified order-13 sextic source produces a full
sign Hadamard of order52. Independent rational matrix multiplication checks
the dense reflection formula, Gram, and cross identities on a small example.

Replay:

```sh
PYTHONPATH=src:. ./.native-env/bin/python -m experiments.scalar_lift.alpoge_mixed_basis_probe --check certificates/scalar-lift/alpoge-mixed-basis-probe.json
```
