# Scaling boundary and reduction scheduling

## What the local experiments establish

The direct Laurent-ideal strategy has a sharp wall between orders 4 and 5 in
the current implementation.

| Computation | Result in this environment |
|---|---|
| Full unsplit order-4 ideal | dimension and 9 minimal primes complete |
| Full unsplit order-5 ideal, 16 variables and 40 ordered row/column equations | saturated dimension did not finish in 120 seconds |
| Row-only order-5 ideal, 16 variables and 10 unordered equations | saturated dimension did not finish in 120 seconds |
| Order-5 translation cover | 65 row-canonical patterns, 3 exact full-symmetry orbits, all 3 certified empty |
| Order-5 square cover | 3,876 row-canonical patterns, 30 exact full-symmetry orbits, completely exhausted |
| `R2+R3` alternatives | odd-order and square-power filtering leave 43 canonical jobs covering 73 literal branches; all empty |
| Regular-pentagon alternatives | square-power filtering leaves 41 canonical jobs covering 528 compatible labelings |
| Surviving order-5 assembly | 6 point presentations, 24 after Galois completion, one Fourier monomial class |

The complete computation is run by
`scripts/run_torsion_order.py -n 5 --until-complete`.
The same generic entrypoint and stopping rule handle `N=4` and `N=6`; none is
selected by a special closure frontend.

## The bottleneck that was removed

The old scheduler did this:

1. scan a cover branch;
2. factor many generators in separate Singular processes;
3. generate a weight-specific short mod-2 cut;
4. ask Singular separately whether every candidate equation was already in
   the ideal;
5. only much later ask whether the resulting branch was empty.

The exact membership calls were only an optimization. Adding a relation that
is already implied does not alter the variety, and the normalized visited-state
key prevents a loop. On an order-5 square branch, those calls repeatedly
recomputed the same Gröbner basis.

The large-frontier scheduler now follows this order:

1. reject unit ideals and exact translated-lattice contradictions;
2. reduce every nonbinomial generator modulo the visible lattice;
3. row-reduce Laurent-unit-normalized generators and feed newly exposed
   monomials or binomials back into the lattice step;
4. align sparse generator pairs by Laurent monomials when a support-count
   test shows that exact cancellation could leave at most two terms, then
   feed every verified result back into the lattice step;
5. apply a visible toric-univariate quotient if one is exposed;
6. otherwise use the recursive exact relation catalog;
7. quotient labeled embeddings by the exact selected-row stabilizer;
8. keep one polygon labeling per exact pointwise coefficient-stabilizer
   Galois orbit and store the omitted orbit size;
9. impose disjoint direct sums of prime polygons simultaneously, so their
   combined binomial lattice is visible at once and the same Galois subgroup
   acts diagonally on the labeling product;
10. split a fully binomial leaf into all primitive translated-lattice
   components;
11. reserve saturation, Gröbner bases, and minimal-prime computation for
   terminal assembly rather than putting them in every search node.

The production order-six frontier also made a bounded Singular
emptiness query only after structural reduction and, by default, only at
active rank at most 21. A positive unit-ideal answer is an exact certificate.
A negative answer or timeout never prunes a branch, and the unchanged ideal
is not queried again. In recorded samples from the order-six frontier,
seventeen of twenty rank-18 leaves and all ten sampled rank-20 leaves were
certified empty in about a tenth of a second each. Sampled rank-21 leaves
instead timed out and were retained. A separate twenty-leaf hard sample timed
out at both 0.2 and 0.5 seconds. The production gate is rank 21 with a
0.2-second bound: fast certificates are retained while inconclusive probes
return promptly to exact structural reduction. A measured low-rank tier uses
one second at active rank at most 13.  On the sealed hard-tail frontier,
30/30 sampled rank-11--13 leaves certified empty in 0.126--0.781 seconds;
the old 0.2-second production bound missed examples in every one of those
ranks.  Sampled rank-14--18 hard leaves remained inconclusive even at 1.5
seconds, so the longer default is not extended to them.  These timings choose
when to ask Singular, not what counts as proof: every timeout retains the
unchanged branch.

The first full 3,000-node round with that tier visited 78,000 nodes across 26
open branches.  It reduced the sealed frontier from 755 to 634 leaves, added
1,058 terminal survivors, and peaked at 2.49 GiB process-tree RSS.  The
preceding round at the same node budget reduced the frontier by 44 leaves;
the measured tier reduced it by 121.  This comparison is evidence for the
scheduling choice, not by itself a classification. At that stage the frontier
still had to reach zero and the final component digest had to be reproduced;
the completed primary and replica campaigns described below later satisfied
both conditions.

Process scheduling is also part of the fix. Each unresolved exact relation
branch is one process-pool task, and the pool persists across rounds. This
writes progress as soon as a branch finishes, preserves exact-expansion
caches, and prevents a hard branch from holding either a static multi-pattern
shard or its sibling branches hostage. Completed branches are validated by
the collector and are not rescheduled. Repeated exact relation realizations
are cached within each worker across successive tasks.

## Order-5 diagnostic

For the zero-sign square branch, every coefficient ratio has odd order.
Consequently the exact length-5 relation has only the regular-pentagon
alternative: every `R2+R3` construction would force an odd-order term ratio
to equal `-1`. Signed representatives whose coefficients already contain
`-1` may still have `R2+R3` branches; the invariant test decides this per
relation.

Over `Q(zeta_5)`, nine resulting equations split completely in the
coefficient field. The implementation now recognizes such a split without a
root-order search: it evaluates the polynomial on the finite root-of-unity
group already contained in the field and proceeds only if the number of
distinct exact roots equals the degree.

After normalizing one nontrivial row to

```text
(1, zeta_5, zeta_5^2, zeta_5^3, zeta_5^4),
```

Singular finds a zero-dimensional saturated ideal with exactly six minimal
points. They are the `3!` permutations of the remaining Fourier rows.

## Why this closes order 5

For the nauty-chosen sign representatives, the termwise square companion
removes pentagon labelings and pair-plus-triangle embeddings that cannot also
satisfy the square-cover equation.  Joint stabilizer/Galois canonicalization
then represents the 528 compatible pentagon labelings by 41 exact jobs.
Odd-order and square-power filtering leave 73 literal pair-plus-triangle
branches, represented by 43 jobs. With the three translation jobs this is
87 exact ideals covering 604 symmetry-weighted branches.

The three translation jobs and 83 of 84 square jobs are empty. The remaining
square job assembles to six exact points. Completing their Galois orbits
gives 24 presentations; the dephased exponent cross-difference canonical
form puts all of them in the Fourier monomial class. Consequently `BH(5,q)`
has one equivalence class if `5 | q`, and no class otherwise.

## Order 6: complete all-exponent Butson classification

The square-cover frontend has now been exhausted combinatorially at order 6.
There are

```text
2^25 = 33,554,432                  dephased sign patterns,
C(36,5) = 376,992                 patterns after sorting the five nonfirst rows,
242                               row/column/transpose orbits.
```

The last number is an exact nauty canonical-label computation, not an
invariant bucket.  A sign pattern is encoded by its mod-2 row/column
cross-ratio tensor as a colored incidence graph.  The same encoder gives 10
orbits at order 4 and 30 at order 5, agreeing with the independent exhaustive
orbit traversal at both regression orders.

Before using coefficient parity, the general weight-six catalog has the
exhaustive 40-way split:

- 15 instances of an `(R5:R3)` relation;
- 10 partitions into two triples;
- 15 partitions into three pairs.

The former 4,989-job/9,680-embedding ledger expanded all 40 shapes on every
sign orbit. The relation engine first tests relative coefficient orders.
Whenever they are all odd, every shape containing `R2` is impossible. It then
checks whether the ideal visibly contains the termwise square
`t_1^2+...+t_6^2` of the same source. Every numerical realization of
`(R5:R3)` has nonzero square sum, while each `R3` component has zero square
sum. Thus the zero sign pattern has one job of orbit size 10, covering only
`R3+R3`; its 15 glued embeddings are certified impossible before algebra.
Signed patterns need not display the same termwise square and may retain the
other types.

Across all 242 representatives the current frontend has 1,454 exact algebra
jobs of total symmetry weight 2,690. These totals are recorded by the
manifest-hashed expected-branch index and regenerated from the sign manifest
and relation engine; the historical unfiltered numbers are not a completion
invariant.

A measured v8 continuation advanced 13 hard branches by 1,000 nodes each,
closed none, and grew the one branch with a recorded baseline from 39 to 89
frontier leaves. The later scheduler leaves the established
lattice/univariate order unchanged but considers a weighted power-companion
source first among relations of equal weight. On a clean zero-sign pilot its
first eight 1,000-node batches had frontier sizes

```text
169, 144, 128, 105, 129, 143, 116, 86.
```

The pilot filtered 80 incompatible relation types and visited the same scale
of tree several times faster than the old ordering. These figures are
historical scheduling evidence, not completion counts.

The next order-independent reduction checks every visible translated
character `x^v=c`. A torsion branch requires `c` to be a root of unity, and
an odd-order branch requires its order to be odd. On the sealed v11 hard-tail
checkpoint this exact test rejected 746 of 1,730 active leaves immediately
(43.1 percent), before any new relation search. Version 12 migrates the
remaining prefix and applies the same check during assembly to terminal
survivors inherited from older versions. This is a reduction statistic, not
a completeness claim.

The resulting v12 square campaign is now exhaustive. All 1,454 expected jobs
are present, accounting for 2,690 symmetry-weighted relation branches. The
frontiers visited 11,348,216 nodes and retained 112,561 terminal survivors.
Exact branch assembly closes all 1,454 jobs with zero unresolved leaves and
emits 104,380 branch-deduplicated component presentations: 96,186 points,
7,909 curves, and 285 surfaces. These are pre-global presentation counts;
they are not numbers of inequivalent or irredundant components.

The translation cover is smaller. Its singleton test leaves 4,370
row-canonical patterns, which collapse to 24 exact sign orbits.  The first
local Singular pass resolves all 24 jobs: 20 are empty, while four have
dimensions 2, 1, 0, and 1 and respectively 6, 4, 8, and 15 minimal primes.
The general frontier path now terminates all 24 translation jobs (20 empty,
4 nonempty). Residual binomial primes are split geometrically after their
minimal-prime calculation, rather than accepted merely because they are prime
over the working coefficient field. The corrected exact assembly closes all
24 jobs with zero unresolved leaves and emits 134 serialized primitive
presentations. Exact primitive-coset normalization joins these presentations
to the square cover.  The completed primary run takes 104,514 presentations
as input and removes 103,085 exact primitive-key duplicates before the
expensive global steps.  Galois and transpose completion generate 15 and 366
presentations, respectively.  Exact row/column monomial equivalence gives 10
classes before containment; six are proper subcosets of another class,
leaving four irredundant classes of dimensions 2, 2, 1, and 0.  Their
component-union digest is
`ad72fefc0b1c51184d255fb6674bc4415c6ba5e7bda6e2316712eecf7d0d4835`.
Every search and mathematical-completeness obligation in the primary
certificate is true; only the deliberately separate independent-reproduction
flag is false in that first record.  A distinct fresh checkpoint tree closes
the same 24 translation and 1,454 square jobs and produces the same normalized
component-union digest.
Its raw traversal is different---3,887,731 square-frontier nodes and 99,211
terminal presentations---which is expected because those are reducer-path
counts rather than invariants.  The fresh certificate sets mathematical
completeness, independent reproduction, and overall classification
completeness to true.  Its Smith-coordinate fixed-`q` replay matches all 16
determinate order-six cells in the reference table.  Consequently neither
134, 104,380, nor 99,211 is itself a classification count; the invariant
result is the four-class component union above.

The finite sign frontend has now also been replayed independently of the
stored representative lists. All 4,370 translation signs and 376,992 square
signs were recanonicalized with nauty, recovering exactly 24 and 242 orbits
and the same orbit minima. The resulting combined audit digest is
`760dbed18d0ec9f96f6ff1970fd7dd46388c5ec18bdcdcce6a80f559178ad93b`.
Separately, the indexed raw artifacts are sealed by distinct 2,963-entry
roots, `68de5ab1...d826` for the primary and `288aa63c...effe` for the
replica. The latter are byte-integrity commitments, not a replay of historical
reduction transitions.

Exact named-family postprocessing identifies those four irredundant maximal
torsion-coset classes as

```text
dimension 2: F23 = F_{2,3}^{r,s}
dimension 2: F32 = F_{3,2}^{r,s}
dimension 1: H6  = H_6^t
dimension 0: T6  = T_6
```

The two Fourier--Diţă surfaces are transpose partners, while the Haagerup
curve and Tao point are transpose-stable. The four maximal classes are an
irredundant union, not a disjoint partition: proper torsion subcosets can lie
in family intersections.

This also proves automatic regularity at order six. The relation engine did
not assume that every six-term row quotient was a sum of `2+2+2` or `3+3`
cycles: the exhaustive cover included the nonregular `(R5:R3)` type and
resolved every surviving signed instance. Since the final union contains only
the three regular parametrized families and the regular Tao point, no
nonregular order-six Butson matrix remains.
The conclusion concerns root-of-unity entries; the unrestricted complex
Hadamard classification in order six remains outside its scope.

[Banica--Bichon--Schlenker](https://doi.org/10.1016/j.jfa.2009.04.013) had
already published the same three parametrized families and Tao point and said
that an omitted computer calculation removes the regularity hypothesis. The
present contribution is an exact, independently reproduced all-exponent proof
of that stated result, not a new classification statement. The compact records
replay the finite sign cover and fixed-`q` counts and recheck the final ideals
and witnesses; a full reduction-and-assembly audit still requires the retained
multi-gigabyte checkpoints or a rerun. The naming code is generic
postprocessing of primitive torsion-coset certificates. It neither changes the
component union nor adds an order-six case to the frontend; the same cover and
reducer used at orders four and five produced the order-six certificate.

### Why Gram positivity was not the missing order-six filter

For a fully assigned unit-modulus matrix, positivity adds nothing to the
Hadamard equations: its Gram matrix is already required to be `6 I`. During
the all-exponent calculation the variables live in a complex algebraic torus,
where entrywise inversion agrees with conjugation only at torsion points. A
positive-semidefinite inequality on an arbitrary complex point of an
intermediate variety would therefore not be a sound algebraic pruning rule.

Positivity is useful in a different search regime. If a fixed-`q` enumerator
has assigned some columns `A`, the unassigned columns would have to realize
the residual Gram matrix `R = 6 I - A A*`; hence `R` must be positive
semidefinite and its rank cannot exceed the number of remaining columns. That
is a strong finite-search filter, but it does not replace the uniform
root-order torsion-closure argument used here. The filters that materially
closed order six were exact sign-orbit reduction, odd-order parity, visible
termwise-power companions, translated-character rejection, and lattice
closure.

This pilot also found a backend bottleneck.  Saturating a reduced
20-variable ideal by the product of all variables timed out, while successive
saturation by the individual variables finished.  These operations are
equal:

```text
I : (x1 ... xd)^infinity
  = (...((I : x1^infinity) : x2^infinity)...) : xd^infinity.
```

The Singular bridge now uses the successive form.

## Running order 6 locally or on a CPU host

The computation is CPU- and memory-bound.  Nauty canonical labeling,
Singular Gröbner bases and primary decomposition, PARI number-field
arithmetic, and integer-lattice reduction do not use the GPU.  More CPU cores
and RAM may help, but no remote instance is assumed by this project or by the
validation record.

Install `nauty` as well as the native algebra stack. The checked-in exact
manifests can be reused, and the same parameterized runner can be assigned the
number of workers actually available:

```bash
micromamba create -f environment.yml
micromamba activate hadamard
python -m pip install -e '.[dev,native,search]'

.native-env/bin/python scripts/run_torsion_order.py \
  -n 6 --output-directory runs/order6-v3 \
  --run-label local-order6-primary --workers 20 --until-complete \
  --translation-manifest \
    runs/order6-v3/translation-manifest.json \
  --square-manifest \
    runs/order6-v3/square-manifest.json \
  --max-nodes 3000 --max-depth 80 --timeout 180 \
  --memory-cap-gb 30 --host-memory-reserve-gb 2 \
  --worker-process-memory-gb 6.5 \
  --singular .native-env/bin/Singular \
  --frontier-empty-probe-timeout 0.2 \
  --frontier-empty-probe-active-rank 21 \
  --frontier-empty-probe-low-rank-cutoff 13 \
  --frontier-empty-probe-low-rank-timeout 1.0
```

The launcher never treats a requested cap larger than the host as usable
memory. It reserves 4 GiB for the operating system and unrelated processes
by default, then uses the smaller of `--memory-cap-gb` and installed physical
memory minus that reserve. The displayed continuation explicitly lowers the
reserve to 2 GiB and therefore enforces a 22 GiB process-tree ceiling on a
24 GiB host. The requested and effective values,
detected physical memory, reserve, and worker count are written to
`resource-policy.json`; raise the reserve with
`--host-memory-reserve-gb` on a shared host.

The 3,000-node slice is measured rather than guessed. From the same sealed
hard-tail frontier, one complete 20-worker round visited 93,000 nodes, added
1,109 terminal survivors, and peaked at 4.90 GiB process-tree RSS. Relative
to the preceding 1,000-node round it produced more terminal certificates per
node, approximately halved net frontier growth per node, and had slightly
higher node throughput. These are scheduling measurements, not mathematical
evidence; the completeness gate is unchanged.

Every completed frontier and assembly is an atomic checkpoint. Rerunning the
same command validates and reuses it. After the primary certificate closes,
run a second fresh output tree with a distinct label and pass
`--reproduction runs/order6-v3/order-classification.json`.
Mathematical cover closure and independent reproduction are separate fields;
only both together set `classification_complete`.

## Order 7: the theorem gap is gone, the scale gap is not

Order 7 has `2^36 = 68,719,476,736` dephased sign patterns and
`C(69,6) = 119,877,472` row-canonical square-cover patterns. The relation
stage itself is now exhaustive. At weight seven the recursive catalog gives

```text
R7
(R5:R3,R3)
R2 + R5
R2 + R2 + R3.
```

For labeled term positions these give 1, 105, 21, and 105 symbolic
embeddings. The closure runner expands the top-level `R7` into all `6! = 720`
cyclotomic labelings, for 951 covered branches per sign pattern before
parity filtering. On the zero sign pattern all relative term orders are odd,
so the two `R2`-containing types are impossible. The termwise-square companion
also rejects every `(R5:R3,R3)` realization. Exact column/Galois symmetry
therefore leaves one `R7` job with stored orbit size 720.

Thus an order-seven frontend no longer stops because somebody omitted a
weight-seven case. The obstacle is multiplying a large sign-orbit computation
by hard multivariate algebra. The next responsible step is a bounded,
checkpointed nauty shard benchmark, followed by algebra pilots on selected
orbits; a full launch should be based on those memory and wall-time
measurements. The GPU is not useful for the current kernels, though a GPU
machine with many CPU cores and adequate RAM can still be a suitable host.

The exact ordinary-relation catalog has a proved boundary at weight twelve.
That is a theorem boundary shared by all orders, not one specialized frontend
per order. A branch above that weight can still close through another
reduction, but the catalog alone cannot mark it exhaustive.
