# A basis-independent obstruction for uniform scalar restrictions

For the two fixed character sources, **no uniform scalar-restriction tower
can end in equally sized ±1 blocks**, regardless of intermediate fields,
uniform scalar bases, final frames, or representation multiplicity.
This is stronger than failure of the registered frames. It does not exclude
nonuniform coordinate maps, and no row/column gauges are part of this test.

Let the actual scalar alphabet be `U`, and suppose it has an exact relation

\[
\sum_{z\in U} c_z z=0,\qquad c_z\in\mathbb Z,
\qquad \sum_z c_z\equiv1\pmod2.
\]

Every uniform rational linear block map `B` preserves this relation. If all
`B(z)` have entries ±1, reducing the matrix relation modulo 2 gives the
all-one matrix equal to zero, a contradiction. A composition of scalar
restrictions and uniform left/right changes of scalar bases is still such
a linear map. Repeating the representation does not change the argument.

The actual order5329 phase24 and order58081 phase120 sources both contain
`1`, `zeta6`, and `zeta6^-1` in row zero. They satisfy

\[
1-\zeta_6-\zeta_6^{-1}=0.
\]

The coefficient sum is −1. The certificate records the exact source hashes,
literal row/column occurrences, and the zero vector in the integral
cyclotomic basis. Thus these are source-specific, exact exclusions, not
time limits or numerical residuals. Their weighted field restrictions can
still be computed as intermediate matrices; those computations cannot
remove the stated obstruction to a uniform sign endpoint.

The implementation additionally gives a **complete test for this particular
odd-relation obstruction** whenever the phase has no prime factors outside
2, 3, and 5. It uses the already proved
[complete minimal-sum catalogue](VANISHING_KNOWN_CATALOGUE.md), not a
modulo-2 kernel test:

1. Replace negative coefficients by occurrences of the opposite root,
   obtaining a positive vanishing multiset in `U union (-U)`.
2. Decompose it into minimal vanishing multisets. An odd total weight forces
   at least one odd-weight minimal member.
3. Square-free relative conductor puts that member, up to rotation, in
   the complete phase30 catalogue. The code exhausts its rotated embeddings
   in `mu_lcm(phase,30)` and maps every matched root back to a signed actual
   source root.
4. Independently verify the resulting relation over the integers in the
   original cyclotomic power basis.

The classification inputs are Lam–Leung, Corollary 3.2, and
Christie–Dykema–Klep, Theorem 3.2, cited precisely in the linked catalogue
document. Other phases can still receive a verified positive witness, but
an unsuccessful match is explicitly incomplete there. A finite embedding
budget is also reported as `SKIP`. No `SKIP` result asserts a sign endpoint
exists; it simply supplies no obstruction.

The fixed order1000 phase8 conference source and order35 phase120 source
have no such odd relation. They are retained as controls; their already
verified positive scalar transfers remain valid. No source was constructed
or altered by this audit.

The separate full-registry manifest
`certificates/scalar-lift/scalar-restriction-all-sources.json` audits all
20 currently registered fixed inputs: the two character sources are
excluded and the other 18 are retained by this gate. The retained list is
not a claim that a real endpoint exists. Replay it by adding
`--all-registry` and using that manifest with `--check` below.

The implementation is
`experiments/scalar_lift/scalar_restriction_obstructions.py`, with API
`alphabet_obstruction(phase, used)`. The saved report is
`certificates/scalar-lift/scalar-restriction-obstructions.json`.

```sh
PYTHONPATH=src:. ./.native-env/bin/python -m experiments.scalar_lift.scalar_restriction_obstructions \
  --check certificates/scalar-lift/scalar-restriction-obstructions.json
```
