# Forward Butson descent: implementation and replay

The mathematical formulation and generator conventions are in
[FORWARD_SCALAR_DESCENT.md](FORWARD_SCALAR_DESCENT.md). This implementation
starts with a complete verified source and applies the explicit edges

\[
 BH(n,pq)\longrightarrow BH(pn,q),\qquad p\text{ prime},\quad p\mid q.
\]

This is the known
[Ó Catháin–Swartz construction](https://arxiv.org/html/1904.10771#S1.Thmtheorem1.1),
not a new construction theorem. Intermediate Butson matrices are a
sufficient forward subfamily; the general metric-retaining formulation
does not require every intermediate stage to be flat.

`experiments/scalar_lift/forward_butson_descent.py` supplies:

- `edge_parameters`, `scalar_blocks`, and `scalar_certificate`: exact
  prime-step actions, cyclotomic polynomial identities, and every scalar
  cross-Gram identity.
- `ForwardTower`: a nested literal matrix whose entries can be obtained
  without allocating its full expansion.
- `verify_source` and `proof`: source Gram verification and proof records
  tied to the precise starting matrix hash, order, and root alphabet.

The default `generator_mode="primary"` uses the generator of the
\(p\)-primary field factor. `generator_mode="primitive"` instead uses
the full primitive upper root. Both are independently certified at
phase 12. Their literal blocks differ; no arbitrary partial-alphabet
preservation is asserted. In primary mode the selected generator has
order \(p^a\), while the separately represented primitive upper root
has order \(Q\). Both cyclotomic polynomial checks are recorded.

The certificate independently verifies the complete Grams along

\[
 F_8\longrightarrow BH(16,4)\longrightarrow H_{32},
 \qquad
 F_{27}\longrightarrow BH(81,9)\longrightarrow BH(243,3).
\]

The second terminal matrix remains complex; `sign_entry` rejects it.
The high-order control replays the existing complete conference-family
\(BH(1000,8)\), including its actual order-500 skew-conference identities,
symbolic order-1000 Gram, saved hash, and primitive eighth-root cross
ratio. It gives the fully defined forward chain

\[
 BH(1000,8)\longrightarrow BH(2000,4)\longrightarrow H_{4000}.
\]

Its final Gram follows exactly from the source and scalar certificates.
The full 4000×4000 matrix need not be allocated to define or access any
of its sign entries. Order 4000 is already known; no new order is claimed.

```python
from experiments.scalar_lift.even_petrescu_conference_family import source_exponents
from experiments.scalar_lift.forward_butson_descent import ForwardTower

source = source_exponents(1000, 8)
tower = ForwardTower(source, 8)
assert tower.orders == [1000, 2000, 4000]
assert tower.phases == [8, 4, 2]
value = tower.sign_entry(3999, 3998)  # A literal integer +1 or -1.
```

`entry(row,column,stage)` gives the root exponent at a chosen stage.
`materialize(stage,max_entries=...)` enforces an explicit allocation
budget. `proof(source_proof)` requires a matching verified source hash;
an arbitrary root-exponent array alone does not establish orthogonality.

Replay:

```sh
PYTHONPATH=src:. python \
  -m experiments.scalar_lift.forward_butson_descent \
  --check certificates/scalar-lift/forward-butson-descent.json
```

Tests check the complete certificate, agreement with the existing
color-separated quadratic blocks, literal lazy versus dense entries,
the actual final sign Gram, source-hash binding, and rejection of an odd
terminal alphabet as a real matrix. No inverse search or native solver
is involved.
