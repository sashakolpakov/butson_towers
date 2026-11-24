# Exact norm descent: the invariant and the remaining search

The user's correction is right: an actual Butson matrix already supplies
the complete norm identity. Scalar restriction should transport that
identity, rather than search for it again. The previous palette/polar
search and the search imposing a skew matrix's quadratic norm equation
were searching for orthogonality; they were not an executor of this
invariant.

`experiments/scalar_lift/exact_norm_descent.py` is a new exact reference
executor. It accepts an actual complete Butson source, verifies its norm,
and carries the operator **and its Hermitian form** through every edge.
It never computes an objective measuring a Gram error.

## What restriction preserves

Let the current field be (E), with conjugation, and let the actual
current state satisfy

\[
 M^*GM=nG.
\]

Initially (M=H\in BH(n,q)) and (G=I_n). For a restriction (E/K)
with basis (e_1,\ldots,e_d), the new operator is literal block
multiplication (M'=\rho(M)). The new form is

\[
 G'_{(a,u),(b,v)}=
 \operatorname{Tr}_{E/K}(\bar e_uG_{ab}e_v).
\]

Consequently (M'^*G'M'=nG'), exactly. If (G_e) is the relative trace
Gram of the basis, an implementable formula is

\[
 G'=(I\otimes G_e)\rho(G).
\]

**The incoming form cannot be reset to an identity or discarded.** This
formula handles its nonconstant field coefficients. Trace transitivity
proves that a tower and a direct restriction agree after their exact
basis change. No new orthogonality equations occur along the tower.

The existing `relative_scalar_restriction.py` already provides the exact
relative multiplication tables, basis maps, trace Gram, and adjoint
certificate. The new executor composes those data, including nonidentity
incoming forms.

## Normalization and signs are separate

At a rational endpoint choose two nonsingular rational coordinate bases

\[
 B=T^{-1}MS,\qquad
 G_{\rm out}=T^TGT,\quad G_{\rm in}=S^TGS.
\]

Every choice gives the automatic identity

\[
 B^TG_{\rm out}B=nG_{\rm in}.
\]

To obtain the ordinary target norm (ndI), sufficient endpoint conditions
are (G_{\rm out}=cI) and (G_{\rm in}=dcI), for one positive rational
(c). These are explicit lattice/frame data, not a consequence of field
restriction alone. Having obtained them, search **only** for the entries
of (B) to be signs. In particular, a common-basis similarity alone
keeps the determinant's norm scale (n); it cannot provide the required
factor (d>1).

There is also a useful block realization that keeps the ordinary norm
exact throughout the terminal search. Start with a verified rational
star representation

\[
 \rho(\bar z)=\rho(z)^T
\]

of the source field (or the required scalar subalgebra), a rational norm
frame (WW^T=dI), and rational orthogonal matrices (U,V). Set

\[
 B_z=U\rho(z)VW.
\]

Then (B_zB_w^T=dU\rho(z\bar w)U^T). Substitution into the **fixed,
actual** source gives (\mathcal B\mathcal B^T=ndI) at every iterate.
No scalar-norm or matrix-Gram objective is needed. An arbitrary trace
form need not admit this rational normalization: its representation and
frame must be supplied and checked, rather than presumed.

An economical restricted branch anchors one used source root (z_0) at
an actual sign Hadamard frame (F), (FF^T=dI), and varies a rational
orthogonal (Q):

\[
 B_z(Q)=Q\rho(zz_0^{-1})Q^TF.
\]

It has the same exact complete norm at every iterate and (B_{z_0}=F).
Only the roots **actually occurring** in the gauged source require sign
blocks. Fixing one particular (F) deliberately restricts the search;
it is not a completeness claim over all frames. A rational matrix with
fixed row norm (d) is a sign matrix if all its entries are odd integers.
Thus, for a block written (A/q), the terminal condition is the exact
congruence (A\equiv q\pmod {2q}), entrywise. This gives an alphabet
objective independent of a Gram objective.

## Concrete reuse, and what this does not solve

`ambient_field_candidates.source_coordinate_view` already constructs
an actual **integer** source-derived endpoint with complete Gram (ndI)
from verified sixth/twelfth-root sources. Its rational norm frame uses
the four-square identity for (r), with (d=4r); its compact output is
the exact operator (N_0\otimes I_r). It has zeros, hence is not a
Hadamard. This is a valid norm-exact starting object for basis/lattice
search, whereas a random sign array with a large Gram residual is not.
Its auxiliary field action must be retained and transported if a child
is claimed to arise from that full field, rather than just from the
source's scalar subfield.

`inverse_blocks.known_order_twelve_pair()` supplies a further exact
small normalization control. The existing known frames at width668 can
anchor larger source-derived endpoints without optimizing their Gram.
The open problem is to find allowed rational lattice changes making
all used blocks odd integral. No efficient complete search for these
changes has been established. Generic rational rotations are not an
efficiency argument.

Previously proved family obstructions still apply: at width (4r), a
sextic flat pair needs genuine nonintegrality at2 and denominators at
the relevant odd prime (r). Dyadic changes, 2-integral changes, the
eight complete Hadamard-diagonal conjugator families, and the old
two-term Paley chart cannot simply be restarted under a new name.
Norm-preserving rational isometries outside those families remain
legitimate candidates, but their alphabet search needs a new arithmetic
parameterization or a justified bounded enumeration. This executor
corrects the invariant; it does not supply such an efficiency theorem.

## Persistent node and tested control

A genuine branch stores or references: the actual source and its exact
norm certificate; the current field inclusion and exact basis maps;
the restricted operator; the full incoming/transferred Hermitian form;
the involution; input/output lattices and their forms; the used scalar
palette; the parent hash; and terminal normalization/flatness results.
Nonflat intermediate matrices are retained as weighted norm objects,
not described as new Butson matrices. Large nodes should be lazy block
DAGs referring to source and scalar tables, not copies of the full
expanded matrix. `ExactNormState` is deliberately a small dense reference
implementation for validating that production representation.

The positive control starts from the actual Fourier (BH(4,8)), follows
(8\to4\to2), proves exact agreement with direct restriction after its
basis permutation, and uses unequal terminal bases to recover known
(H_{16}). A separate (30\to6\to2) test transports a genuinely
nonconstant incoming metric. A rational isometry control preserves the
complete ordinary norm while destroying sign flatness, explicitly
testing that the two certificates are separate. Four tests pass.

```
PYTHONPATH=src:. OPENBLAS_NUM_THREADS=1 .native-env/bin/python -B \
  -m experiments.scalar_lift.exact_norm_descent \
  --check certificates/scalar-lift/exact-norm-descent-control.json
```

The certificate and controls assert no previously unknown order.
