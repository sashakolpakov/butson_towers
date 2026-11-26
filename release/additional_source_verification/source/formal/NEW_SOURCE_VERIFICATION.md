# Additional finite sources and Spence families

The new modules leave the earlier literal source modules and submitted
verification bundles unchanged. From the repository root, check the
literal-data bindings with:

```sh
python3 scripts/check_lean_additional_sources.py
```

From `formal`, check the complete new proof dependencies with:

```sh
lake build ButsonTowers.NearWilliamsonEndpoints ButsonTowers.AdditionalSpenceFamilies ButsonTowers.AdditionalResults
lake env lean AuditAdditionalSources.lean
```

`Data/NearWilliamson/Core.lean` proves that replacing the first circulant
by its product with the reversal `R[x,y]=1` iff `x+y=-1` preserves its
Gram contribution and produces an amicable quadruple when the other
three rows are symmetric. The width-69 source uses four ordinary
circulants. The twelve separate literal modules check signs, the needed
symmetries and every summed periodic correlation by `decide +kernel`.
`NearWilliamsonEndpoints.lean` applies the existing T631 substitution
theorem to prove all twelve endpoint orders.

The six `Data/AdditionalSpence/Q*.lean` modules check the signed literal
starter rows, both strong norm identities and each complete Paley row.
`AdditionalSpenceFamilies.lean` applies the existing generic Spence
induction for every natural-number parameter. The statements concern
existence; novelty labels are external literature claims.
`AdditionalResults.lean` proves the four further finite specializations
and the shifted family `15131380*109^r`.

The standalone audit prints the dependencies of all 48 additional named
theorems. Only the standard `propext`, `Classical.choice` and `Quot.sound`
axioms may occur. No native evaluator, unchecked axiom, or omitted proof
is used. A portable package must include these new modules, their local
import closure, both public input JSON files and the new literal checker;
the earlier frozen packages retain their original scope.
