# Exact forward frame and gauge search

`experiments/scalar_lift/forward_frame_search.py` provides terminal sign frames
for genuine Butson inputs. It does not replace the input's exact Gram proof.
Its 16 frames cover phases 3, 5, 6, 10, 12, 20, 24, 40, 48, 80, 96 and 120.
The phase-120 frames have width 32, so an admitted source of order 1499 would
give a real matrix of order 47968. No such source is constructed here.

For each frame the module constructs every integer block `B_e`, including
the excluded, nonflat phases. With `F=B_0` and width `d`, it verifies

\[
FF^T=dI,\qquad C=B_1F^T,\qquad CC^T=d^2I,
\qquad CB_e=dB_{e+1}.
\]

It also verifies every cyclotomic power-basis relation among the `B_e`.
Consequently `rho(zeta)=C/d` is an orthogonal representation of the field,
`B_e=rho(zeta^e)F`, and
`B_e B_f^T=d rho(zeta^(e-f))`. A verified source Gram therefore transfers
exactly. Flatness is checked separately on every entry of every block.
The low-phase frames compose the literal twelfth/fifth transfers in
`src/hadamard/dense_lift.py` with the Gaussian-to-sign transfer.

`global_affine_fit` exhausts unit multipliers and global root phases using
the source's actual alphabet. `row_column_gauge_search` optionally solves
the larger family using one CP-SAT worker and an explicit time/entry budget.
The frame alphabets exclude one residue at each of the primes 3 and/or 5.
Their row/column constraints therefore separate by the Chinese remainder
theorem. Dyadic coordinates need not change. A Galois multiplier gives no
additional feasibility in this larger search: multiplication rescales the
gauges, and a global translation restores the omitted color. Only proven
`INFEASIBLE` excludes this gauge family; `UNKNOWN` and `NOT_RUN` retain it.

For the implicit character sources, `character_cubic_obstruction` checks the
actual eight-point source submatrices and reuses the complete finite cubic
gauge census. It does not materialize their enormous square matrices.

The certificate records a genuine BH(15,120) scrambled by nonconstant row
and column gauges. Every global affine transformation fails to fit its
alphabet, while explicit inverse row/column gauges recover an admitted
source and an independently checked real H480. The tests additionally
recover gauges through CP-SAT and check the resulting full sign Gram.
This is a construction control, not a new order.

APIs used by the forward search are `frame_options`, `frame_proof`,
`frame_block`, and `check_frame`. The block tables are read-only and the
endpoint can be accessed lazily. Failure for these frames never excludes
other frames or a restriction path with weighted, nonflat intermediate
matrices; see [the forward formulation](FORWARD_SCALAR_DESCENT.md).

Replay:

```sh
PYTHONPATH=src:. ./.native-env/bin/python -m experiments.scalar_lift.forward_frame_search \
  --report certificates/scalar-lift/forward-frame-search.json
```
