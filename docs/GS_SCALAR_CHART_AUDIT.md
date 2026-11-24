# Goethals–Seidel as a joint scalar-restriction chart

This note treats one **sixth-root source family**, not a necessary phase
for Hadamard matrices, Goethals–Seidel arrays, or scalar restriction.
The selected order-7 and order-13 inputs have alphabet
`{+1,-1,+zeta6,-zeta6}`; that choice explains the cubic-pair equations below.
The same candidate order4676=7*668 is degree-compatible with both the
verified `BH(7,6)` input and the genuinely twelfth-root `BH(7,12)` input.
The latter has different scalar and alphabet constraints and must not be
forced into this sextic-pair chart. Degree compatibility alone does not
certify a flat restriction.

For the genuine twelfth-root order-7 input, an explicit global shift by
`zeta12^-1` gives the palette `mu4 union zeta12*mu4`. Write
`J=rho(i)`, `R=rho(zeta12)`. The corresponding orthogonal scalar actions
satisfy `J^2=-I`, `JR=RJ`, and `R^2-JR-I=0`. The required sign blocks are
`F, JF, RF, JRF` and their negatives, with `FF^T=dI`. These are different
constraints from the sextic pair below. For general twelfth-root inputs,
the exact used palette, including its chosen gauge, must be tested; this
four-block palette must not be assumed universally.

There is a concrete compatible use of Goethals–Seidel already implemented
in [joint_gs_223.py](../experiments/scalar_lift/joint_gs_223.py). Its
`goethals_seidel_array(sequences)` assembles four circulant blocks.
`verify_sequences(states,prime,multiplier)` checks the complete integer
coefficient correlations and then both full sign-block Grams. This need
not become a new source-construction operation in the restriction tree:
the array can instead parameterize the unknown scalar-coordinate images
inside the joint sign/Gram intersection problem.

The useful unknown object at width `d=4r` is a real orthogonal design
`OD(4r;3r,r)`, also written with its two weights in the opposite order.
Write it as `D(x,y)=xU+yV`, with complementary ternary supports, and require

`UU^T=3rI`, `VV^T=rI`, `UV^T+VU^T=0`.

Then `F=U-V` and `G=U+V` are sign matrices with

`FF^T=GG^T=4rI`, and `GF^T+FG^T=4rI`.

The exact rational matrix `R=GF^T/(4r)` is orthogonal and satisfies
`R^2-R+I=0`. Therefore `B(1)=F`, `B(zeta6)=G` defines the required two-basis
sextic scalar map; the other blocks are `G-F,-F,-G,F-G`. Substituting the
already verified sixth-root sources is scalar restriction in those bases.
No new Butson source is assumed or constructed by that substitution.

At width 668 the missing condition is specifically
`OD(668;501,167)`; the order of these weights is immaterial if the variables
are relabelled consistently. This would supply the pair for source orders
7 and 13, targeting 4676 and 8684. Width 716 requires
`OD(716;537,179)` for the order-13 route to 9308. Width 892 requires
`OD(892;669,223)` for the order-3 route to 2676 and the order-13 route
to 11596. These are conditional targets; no such flat pair at these
widths is certified here.

Alpöge's published H668 supplies one F, not the second G or this formal
design decomposition. If F is held fixed, the chart additionally imposes
`U-V=F`. Every border row, border column and corner must satisfy the same
three coefficient identities. A valid signed Goethals–Seidel matrix by
itself does not establish them. The earlier
[M-structure product audit](https://github.com/sashakolpakov/complex-hadamard/blob/976cede671f40c810d593c8d5815e26e11a03ed6/experiments/scalar_lift/mstructure_668_product.py)
therefore checked full target-array row Grams instead of assuming that its
matrix-valued blocks commute.

The recovered R supplies a sextic field representation. Identifying it
with a particular previously chosen full ambient
`L_r Q(zeta12)` representation, with its prescribed compatible metric and
auxiliary operators, is a further condition. The formal design alone must
not be claimed to certify that full ambient-field embedding.

## Preserving the norms before imposing flatness

Both frames can vary while keeping all three coefficient Gram identities.
Let `a,b` be four real cyclic sequences of length `r`, with unnormalised
Fourier transforms `a_hat,b_hat`. At every character impose

```
sum_i |a_hat_i|^2 = 3r,
sum_i |b_hat_i|^2 = r,
Re(sum_i a_hat_i conjugate(b_hat_i)) = 0.
```

These are exactly the Fourier transforms of the required summed periodic
correlations. Thus their Goethals–Seidel arrays `U,V` satisfy the three
coefficient identities above even when their entries are nonflat. No
prescribed Hadamard frame is needed. At a nonreal character the normalised
pair is a real orthonormal two-frame in `C^4`, with inner product
`Re(x* y)`; complex orthogonality would impose an extra condition. The
zero character is real, and opposite characters must be conjugates.

The remaining terminal equations are just, at every sequence coordinate,

```
a_i b_i = 0,     a_i^2 + b_i^2 = 1.
```

Indeed `(a_i +/- b_i)^2 - 1 = (a_i^2+b_i^2-1) +/- 2a_i b_i`.
The signed entries of `F,G` are selected from these sums and differences,
so this identity factors their flatness equations through the scalar chart.
Homogenising replaces `1` by `h^2`. The norm parameterization does not
establish that this real algebraic intersection has a rational point.

For any point of the norm parameterization, set `d=4r` and
`R=G F^T/d`. Then `R R^T=I`, `R+R^T=I`, `R^2-R+I=0` and `R F=G`.
At an integer terminal point this is an exact rational scalar action,
with conjugation represented by transpose, frame scale `d`, and the
verified source scale retained separately. In particular the shifted
order-3 source, with diagonal exponent 3 and off-diagonal exponent 1,
gives blocks `-F` and `G` and Gram `3d I`.

At a rational point, the trace metric can be specified without assuming a
canonical basis is Euclidean. Put `S=2R-I` and `xi=sqrt(-3)`, acting by
`S`. The Hermitian form
`h(u,v)=(u^T v)/2-(u^T S v)xi/6` on this
`Q(xi)`-module has `Tr(h(u,v))=u^T v`. Its module rank is `d/2` and its
rational dimension is `d`. This retains the incoming Hermitian form;
it does not identify it with the identity form over `Q(xi)`.
Floating-point projections onto the Fourier constraints are search tools;
only exact coefficient, frame, source and entry checks certify an endpoint.

## Real variables and complex torsion are different checks

The formal identity uses **commuting variables fixed by the involution**.
The ordinary real Goethals–Seidel layout does not automatically become a
complex Hadamard array when arbitrary roots of unity are substituted.
For example, take four length-one blocks `(1,i,1,1)`. Every entry in the
ordinary array is a fourth root, and its Gram diagonal is 4, but the Gram
entry `(0,1)` is `2i`. Extra identities or an appropriate complex array
would be necessary.

Likewise, substituting a nonreal phase for x in the real formal design
does not preserve its transpose-based polynomial Gram automatically. The
safe route above derives the rational sextic action from the verified
pair F,G and checks its conjugate/adjoint relations explicitly.

## Exact control rerun

[gs_scalar_chart_audit.py](../experiments/scalar_lift/gs_scalar_chart_audit.py)
replays the preexisting literal `OD(12;9,3)` control. It verifies all three
coefficient identities, complementary supports, the recovered sextic
operator and complete sign Grams after applying the existing source orders
3, 7 and 13. The resulting control orders are 36, 84 and 156, all known;
no new source-search edge or new-order claim is made. The same audit checks
the complex negative example exactly.

[Certificate](../certificates/scalar-lift/gs-scalar-chart-audit.json).
The nine existing GS/T-sequence tests also pass. No large solver or random
search was started in this audit.

```sh
PYTHONPATH=src:. ./.native-env/bin/python -m experiments.scalar_lift.gs_scalar_chart_audit --check certificates/scalar-lift/gs-scalar-chart-audit.json
```
