# Symmetry-first Williamson sources

All ten recent Williamson runs are paused at the user's request for strategy
review. The direct searches took too long without a verified new source;
this is neither a nonexistence result nor a reason to exclude these classes
from a later retry. Complete frozen inputs, checkpoints, random-generator
states and logs are preserved in a verified private archive under
`.local/search-archives`. The earlier eleven one-hot T-input workers also
remain paused. The two unstarted choices are withdrawn from the queue.
No search starts automatically.

The [next selection strategy](https://sashakolpakov.github.io/butson_towers/echelons.html) uses three order
tiers, below 10000, 10000–50000 and 50000–100000, and chooses reusable
ingredients and compatible symmetries before allocating computation.
The proof and the dated below-10000 selection remain available below.

Start with the dated list of 81 unresolved orders, factor each target as
$N=4pt$, and retain a prime $p$ only when a complete T-input of length $t$
is already available. Select compatible multiplier actions and orbit counts
before building an engine. The missing object in this route is a Williamson
source of length $p$; the T-input used in its final substitution is known.

## The actual eight-state alphabet

Seek four even sign functions $w_j:\mathbb Z_p\to\{\pm1\}$, indexed by
$j=0,1,2,3$, with $w_j(0)=1$. Let $A_j=\sum_xw_j(x)g^x$ in
$\mathcal A=\mathbb Q[C_p]$, where $g^p=1$ and $(g^x)^*=g^{-x}$.
Their circulants $W_j$ are symmetric and commute. The source norm is
$\sum_j A_j^2=4p\,1$, equivalently $\sum_jW_jW_j^{\mathsf T}=4pI_p$.

The pointwise product condition follows from this norm. Write
$J_p=\sum_xg^x$, $A_j=J_p-2X_j$, where $X_j$ records the negative positions,
and $K=\sum_j|\operatorname{supp}X_j|$. Then
\[
\sum_jX_j^2=p\,1+(K-p)J_p.
\]
Modulo two, squaring permutes the group basis because $p$ is odd and fixes
both $1$ and $J_p$. Hence
$\sum_jX_j=1+(K+1)J_p$. Its coefficient at zero is zero, so $K$ is even.
Thus $\sum_jX_j=1+J_p$: the number of negative signs is odd at every $x\ne0$.

Put $F(x)=(w_0(x),\ldots,w_3(x))^{\mathsf T}$ and $D=J_4-2I_4$,
where $J_4$ is the all-ones matrix and $e_i$ denotes a coordinate vector.
The eight permitted off-zero vectors are exactly the columns of $D$ and
their negatives. Therefore
\[
F(0)=\mathbf1,\qquad F(x)=\sigma(x)D e_{a(x)}\quad(x\ne0),
\qquad a(x)\in\{0,1,2,3\},\quad\sigma(x)\in\{\pm1\}.
\]
Choose a residue multiplier $\mu$ of even order $m$ and impose
$F(\mu x)=\Pi F(x)$. Here $\Pi$ is the identity, or fixes coordinate zero
and cycles coordinates $1,2,3$. In the latter case require $6\mid m$.
Then $\Pi^{m/2}=I$ and $\mu^{m/2}=-1$, giving the required even functions.
There are $r=(p-1)/m$ nonzero orbits.

## The signs are forced by the axes

Since $D^{\mathsf T}D=4I$, set $b=D^{\mathsf T}F/4$. Then
$b(0)=\mathbf1/2$ and $b(x)=\sigma(x)e_{a(x)}$ off zero. The earlier
inversion obstruction assumed a one-hot value at zero; its hypothesis
fails here. These four entries equal to $1/2$ must be retained in every correlation.

For $s\ne0$, let $E_a(s)$ count residues $x\ne0,-s$ for which
$a(x)=a(x+s)$. Pairing $x$ with $-x-s$ leaves the single fixed point
$-s/2$, so $E_a(s)$ is odd. If
$C_b(s)=\sum_x b(x)^{\mathsf T}b(x+s)$, the same pairing gives
\[
C_b(s)\equiv E_a(s)+\sigma(s)\pmod4,
\qquad \sigma(s)=(-1)^{(E_a(s)+1)/2}.
\]
Thus every axis word determines exactly one sign word compatible with
$C_F(s)=4C_b(s)\equiv0\pmod{16}$. The formula respects multiplier orbits.
After computing these signs, reconstruct all four rows and check every
integer correlation, symmetry, the entry alphabet and the source norm.

The selected charts below pass the row norm $\sum_j(\sum_xw_j(x))^2=4p$
and exact joint orbit-count/product tests. They remain construction targets.

| Target orders | $p$ | $m$ | $r$ | Coordinate action | Status |
|---|---:|---:|---:|---|---|
| 9260 | 463 | 22 | 21 | Identity | Paused for review |
| 2676, 4460 | 223 | 6 | 37 | Identity or three-cycle | Paused for review |
| 2956 | 739 | 18 | 41 | Three-cycle | Paused for review |
| 7924 | 283 | 6 | 47 | Identity or three-cycle | Paused for review |
| 5892 | 491 | 10 | 49 | Identity | Paused for review |
| 3004 | 751 | 10 | 75 | Identity | Unstarted; selection under review |
| 5692 | 1423 | 18 | 79 | Three-cycle | Unstarted; selection under review |

## From a completed source to the endpoint

A verified T-input gives coefficient matrices $X_j$ of an orthogonal design
of order $4t$, with Gram $t\sum_jz_j^2I_{4t}$. Substitute the completed
$A_j$ to obtain the weighted algebra-valued source $Y=\sum_jX_jA_j$.
Its order is $4t$ and its Gram scale is $4tp$. Restrict through the full
regular action $\rho:\mathcal A\to M_p(\mathbb Q)$; conjugation becomes
transpose. The normalized regular trace is $\tau(a)=[g^0]a$ and
$\tau((g^u)^*g^v)=\delta_{uv}$, so its basis metric is $I_p$.
Both factors $\mathbb Q\times\mathbb Q(\zeta_p)$ are retained; in these
components $\tau(u,v)=(u+\operatorname{Tr}(v))/p$.
Every entry of the design is one signed variable. Thus every block of
$\rho(Y)$ is $\pm W_j$, with sign entries, and
$\rho(Y)\rho(Y)^{\mathsf T}=4tpI_{4tp}$. Its order and Gram scale both equal
$N=4tp$. Flatness follows from these actual blocks, not from scalar restriction alone.

The known-source control $p=19,m=6,\mu=8$, with the stated three-cycle,
has 64 axis words and 512 signed words before sign elimination. Exactly
18 literal sources pass all full integer conditions in this class.
Independent enumeration, the native census and an ungauged CP census agree.
This is a control at a known order, not a classification of all Williamson19
matrices or a new Hadamard-order claim.

## Code and exact replay

The retained dense Williamson engines provide a CP-SAT model and an
independent native axis sampler and annealer. Both reconstruct the literal
source and verify every integer correlation before accepting a witness.
They remain usable for a later retry.

The core geometry, forced-sign formula and complete source verifier are in
[`williamson_symmetry`](../experiments/scalar_lift/williamson_symmetry.py).
The [selection certificate](../certificates/scalar-lift/williamson-symmetry-selection.json)
binds the dated 81-order frontier, the seven prime inputs and their nine
strongest coordinate-action charts. Its eight target orders are conditional.
The [positive controls](../certificates/scalar-lift/williamson-symmetry-controls.json)
contain a literal known source, the exact finite census and the small
complete T-inputs used by the endpoint recipes.

Generate a fresh coefficient file, without starting a search:

```sh
python -m experiments.scalar_lift.williamson_symmetry \
  --p 463 --m 22 --action identity --output /tmp/williamson-463.json
python -m experiments.scalar_lift.williamson_symmetry --check-selection
./scripts/run_tests.sh tests/test_williamson_symmetry.py
```

The coefficient generator checks every nonzero shift table, the even
pair coefficients and the odd fixed-point contribution. No large coefficient
tables or unsuccessful-run logs are needed in the repository. A successful
new source still needs its endpoint certificate and literal Lean bridge
before the corresponding new existence theorem is reported as formalized.
The manuscript's completed 34-order and nine-family scope is unchanged.

For a later retry, the bounded example commands are:

```sh
python -m experiments.scalar_lift.williamson_axis_cp \
  --p 463 --m 22 --action identity --target=23,-21,-21,-21 \
  --seconds 60 --seed 20261127 --output /tmp/williamson-463-cp.json
python -m experiments.scalar_lift.williamson_axis_native \
  --coefficients /tmp/williamson-463.json --output /tmp/williamson-463-native \
  --mode anneal --seconds 60 --seed 20261127
```

The CP model imposes exact integer correlations and joint sign-state counts.
The native sampler and annealer compute the forced signs and test the same
integer identities. Their independent acceptance guards reconstruct every
literal source entry. Saved jobs retain their stated budgets, seeds and
source hashes.
