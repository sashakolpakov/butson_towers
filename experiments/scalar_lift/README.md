# Active scalar-reduction work

Start with [the identities and scope](../../docs/SCALAR_REDUCTION.md).
Every starting design must have an exact Gram proof. The active operations
restrict its scalars while keeping the source, metric, alphabet and literal
matrix recipe. General weighted nodes are useful outputs in their own right.

| Task | Entry module |
|---|---|
| New real H7572: exact T631, scalar frame, saved matrix and full integer verification | [Standalone recipe and certificate](../../artifacts/hadamard-7572/README.md) |
| New periodic T127, full T identities and known-order H1524 scalar control | [Standalone recipe and dated scope](../../artifacts/t-matrices-127/README.md) |
| Screen the remaining frontier below 12000 using both odd coordinate actions and bounded orbit counts | `low_gap_target_scan --check`; [finite necessary-condition certificate](../../certificates/scalar-lift/low-gap-target-scan.json) |
| Active SL2 parity-axis generation, exact mod4/mod8 sign lifting, finite pool controls and completed endpoint verification | `hadamard_frame_search`; [method and replay](../../docs/HADAMARD_FRAME_SEARCH.md) |
| Search for length-223 periodic T-matrices, targeting missing order 2676 first | `large_parity_map`, `quadratic_sign_lift`, `exact_t_campaign`; [target, exact checks and replay](../../docs/LOW_ORDER_SEARCH.md) |
| Validated length-739 three-cycle chart targeting 2956 with source width one | `exact_t_campaign`; [map scope and replay](../../docs/LOW_ORDER_SEARCH.md#the-next-lower-gap-order-2956) |
| Identity-action multiplier charts at lengths 223 and 491; conditional orders 2676/4460/11596 and 5892 | `identity_multiplier`, `identity_t_campaign`; [distinct geometry, 30-bit T223 control and replay](../../docs/LOW_ORDER_SEARCH.md#identity-coordinate-action-at-lengths-223-and-491) |
| Exact unrestricted sign searches at lengths 167 and 179 | `unrestricted_parity_map`, `unrestricted_t_sign_model`, `unrestricted_t_sign_sat`, `unrestricted_t_sign_cp`, `unrestricted_t_campaign`; [solver scope and replay](../../docs/LOW_ORDER_SEARCH.md#unrestricted-sign-solvers-for-lengths-167-and-179) |
| Select compact multiplier charts; verified maps at 1951, 3571 and 4219 | `multiplier_target_scan`, `large_parity_map`, `hadamard_frame_search`; [finite scan and target selection](../../docs/HADAMARD_FRAME_TARGET_SCAN.md). Use `exact_t_campaign` and `quadratic_sign_lift` for the larger orbit masks. |
| New real orders 17668/22716/32812 from the verified width-2524 frame and further admitted Butson sources | `reuse_hadamard_7572_frame` |
| Apply all four admitted source orders to a newly completed prime frame | `reuse_hadamard_frame` |
| New real orders 98436/295308 and known order 885924 from complete unreal-source triplings | `tripled_unreal_scalar_endpoints`; [proof and scope](../../docs/TRIPLED_UNREAL_SCALAR_ENDPOINTS.md) |
| 27 exact Williamson sources restricted through cyclotomic group algebras; eight additional new real orders | `williamson_scalar_endpoints`; [inputs, completed orders and prior constructions](../../docs/WILLIAMSON_ORDER_RECORDS.md) |
| 19 further compatible inputs; twelve additional new real orders with full field/module reconstruction | `near_williamson_scalar_endpoints`, `commutative_scalar_algebra`; [source inventory](../../docs/SOURCE_INVENTORY.md) |
| Original five literal Spence families for every integer r>=0; four additional new first expanded orders and known H18629644 | `williamson_spence_family`, `williamson_scalar_endpoints`; [induction and scalar proof](../../docs/SPENCE_SCALAR_FAMILIES.md) |
| Six additional prime recurrences and higher finite members; three additional announced families | `additional_spence_families`; [current family scope](../../docs/SPENCE_SCALAR_FAMILIES.md) |
| Positive prior constructions for all 31 known orders; exact order-885924 inequivalence certificate | `scripts/audit_known_order_provenance.py`, `known_order_comparisons`; [scope and replay](../../docs/sphinx/known-orders.md) |
| Fixed published order-19 and prime-square source palette checks | `additional_sixth_root_source_gates` |
| Symbolic endpoint flatness, nonflat trace towers, real H768 examples and complex order-42/294 families | `projective_flat_descent` |
| Small exact reference for carrying the complete operator and incoming form | `exact_norm_descent` |
| Exact relative fields, actions and trace metrics | `relative_scalar_restriction` |
| Resumable weighted divisor-field trees | `weighted_scalar_tree` |
| Certified repeated-prime root-preserving reductions | `forward_butson_descent` |
| Saved literal matrices and implicit large descendants | `large_butson_descendants`, `large_butson_factory` |
| Complete source recipes and acquired published inputs | `forward_source_registry`, `catalogue_source_acquisition`, `broad_source_acquisition`, `expanded_gauge_inputs` |
| Exact uniform-alphabet obstruction | `scalar_restriction_obstructions` |
| Partial third/fifth-root endpoint reductions | `partial_fifth_descent` |
| Ambient fields and alternate scalar bases | `ambient_field_candidates`, `pure_period_restriction`, `mixed_basis_restriction` |
| Certified fixed endpoint frames and source-gauge admission | `forward_frame_search`, `forward_descent_search` |

The remaining modules supply exact source recipes and arithmetic for these
constructions. Superseded search branches are recoverable from
[Git history](../../docs/RESEARCH_HISTORY.md).
The [complete norm-descent reference](../../docs/EXACT_NORM_DESCENT_EXECUTOR.md)
also handles a nonidentity incoming Hermitian form. Its small dense controls
separate norm transport, Euclidean normalization and sign flatness.

From the repository root, replay the complete known-order descent:

```sh
python -m experiments.scalar_lift.forward_butson_descent \
  --check certificates/scalar-lift/forward-butson-descent.json
```

Run a bounded weighted tree from a verified small source, using fresh paths:

```sh
python -m experiments.scalar_lift.weighted_scalar_tree \
  --source unreal3-q6 --database /tmp/hadamard-scalar-control.sqlite3 \
  --seconds 10 --report /tmp/hadamard-scalar-control.json
```

The [projective endpoint construction](../../docs/PROJECTIVE_FLAT_DESCENT.md)
also gives complete flat families, including 36 independent unit phases at
order 294 and five pairwise inequivalent sixth-root specializations. It also
constructs five pairwise inequivalent real H768 matrices from quaternary
sources, distinguished by binary-code dimensions 17 through 21.

These inventories preserve orthogonality with their recorded metric. Claim a
Butson or real Hadamard endpoint only after the required ordinary Gram scale
and root/sign entries have been proved. The
[completed results and status convention](../../docs/HADAMARD_HIGHER_ORDER_RECORDS.md)
record **65 completed real orders in the main campaign: 34 new and 31 known
recovered orders**, plus the separate new periodic T127 ingredient. The
nine announced infinite constructions have individually audited finite members;
their entire order sets are not labelled new. Earlier source inventories
and reference chains retain their separate known-order control scope.
The identity-action maps and their affine sign controls are search reductions;
they provide no new signed T223/T491 witness or additional completed order.
