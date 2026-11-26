# Previous higher-order research

The 637 superseded higher-order files are preserved as ordinary Git files in
[snapshot 976cede](https://github.com/sashakolpakov/complex-hadamard/tree/976cede671f40c810d593c8d5815e26e11a03ed6).
That snapshot was pushed before the local archive was removed. Those commits belong to the original `complex-hadamard` repository.
This project retains the active verified sources and scalar-reduction engines;
its separation provenance is recorded in `MIGRATION.json`.

Read a particular historical file without restoring the retired branches:

```sh
git -C ../hadamard show 976cede:docs/JOINT_CYCLOTOMIC_INCIDENCE.md
```

To restore a specific file when needed:

```sh
git -C ../hadamard restore --source=976cede -- experiments/scalar_lift/joint_incidence.py
```

The [historical inventory](https://github.com/sashakolpakov/complex-hadamard/blob/976cede671f40c810d593c8d5815e26e11a03ed6/archive/higher-order-2026-11-23.manifest.json)
records the retired paths and original checksums. Historical chart exclusions,
known controls and incomplete searches keep their original mathematical scope.
Ordinary work starts with [scalar reduction](SCALAR_REDUCTION.md).

Subsequent unsuccessful local, Fourier, neighborhood and fixed-basis completion
searches were removed from the construction repository as well. Their source
and finite run records remain recoverable in
[the preceding snapshot, 5b71ca4](https://github.com/sashakolpakov/butson_towers/tree/5b71ca4).
That cleanup removed 228 obsolete files (about 40 MiB) from the working tree.
The active
[lower-order work](LOW_ORDER_SEARCH.md) uses exact determinant coordinates
and modular sign lifting.

The later cleanup is also reversible from the pushed
[snapshot/pre-sphinx-cleanup-20261126](https://github.com/sashakolpakov/butson_towers/tree/0d187ac3eff0e1db61cd1fcf85e0e6d647e70797).
That snapshot includes the construction sources and the first complete Sphinx
manual before the remaining obsolete experiment reports and notes were removed.
The later cleanup removed 20 additional files (about 12 MiB).
