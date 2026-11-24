# Optional higher-root input preflight

The [expanded report](../certificates/scalar-lift/phase-independent-expanded-frontier.json)
adds94 explicitly gauged twelfth-root representatives to the original555
fixed input bindings. The old default report and its source catalogue are
unchanged. This is a separate optional run, not a replacement checkpoint.

The new representatives come from the
[higher-root acquisition audit](../certificates/scalar-lift/expanded-gauge-inputs.json).
Of455 published sixth-root classes,14 already admitted the original
alphabet test and80 previously excluded classes admit a fixed twelfth-root
gauge. Those80 consist of four order9 and76 order13 inputs. The other361
remain excluded in that specific higher-root gauge search. The94 inputs
are new coordinate representatives of known classes, not94 newly
constructed projective Butson classes.

All649 input bindings receive the same full156-order frontier scan and
all distinct used-root anchor tests as the original run:

| Executed quantity | Count |
|---|---:|
| Source/target divisibility checks | 101244 |
| Literal input classes modulo global phase | 629 |
| Distinct raw/anchor alphabets | 1216 |
| Explicit used-entry anchor recipes | 4292 |
| Divisible anchor/target records | 2908 |
| Necessary-test surviving anchor/target records | 1374 |
| Surviving source bindings | 135 |
| Surviving literal source classes modulo global phase | 126 |
| Surviving target orders | 22 |

The source hashes and full-Gram proof hashes bind every representative.
The expanded acquisition certificate is itself hash-bound in the report.
Exact global-phase duplicates are reported; arbitrary row/column
equivalence is not used to prune uniform-alphabet cases. The actual entry
field, anchor alphabet, sign-block width, odd integer relation and existing
standard-metric representation filters are evaluated separately.

The22 retained target orders are unchanged. The only one above20000
is still23960, through the fifth-root order5 source or the quaternary
order10 control. No target basis or sign endpoint is constructed by this
preflight; its retained records are necessary-test survivors only.

The [implementation](../experiments/scalar_lift/phase_independent_frontier.py)
keeps its previous default. Optional entry points are
`all_sources(include_expanded=True)`,
`route_records(include_expanded=True)`, and `report(include_expanded=True)`.

```sh
PYTHONPATH=src:. ./.native-env/bin/python -m experiments.scalar_lift.phase_independent_frontier --expanded-gauges --check certificates/scalar-lift/phase-independent-expanded-frontier.json
```

The original command without `--expanded-gauges` still replays the old555
input certificate exactly. No broad Walsh probe or additional source
expansion was run as part of this optional preflight.
