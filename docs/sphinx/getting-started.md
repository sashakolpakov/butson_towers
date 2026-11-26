(getting-started)=
# Reconstruct a complete small example

A literal length-seven input makes the whole construction small enough to
inspect and check densely. It gives an order-84 Hadamard matrix, used here as
a known-order control. The same functions implement the polynomial identities
used at the manuscript's larger orders.

## Install the arithmetic and build tools

Run commands from the repository root. The base package needs Python 3.11 or
later, NumPy and SymPy. The following creates an isolated environment and
installs the tested documentation frontend as well:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -e '.[docs,dev]' -c docs/requirements.txt
```

The pinned frontend consists of Sphinx, MyST and Furo. Exact native searches
add a C++17 compiler; CP-SAT adds the `search` extra; SAT additionally needs
the `cryptominisat5` executable. None is needed to replay the finite example
below. See {ref}`sign-solvers` before starting a solver run.

## Decode the finite input

The seven integers are

```{doctest}
>>> states = [0, 6, 0, 4, 4, 5, 2]
>>> [(s // 2, 1 - 2 * (s % 2)) for s in states]
[(0, 1), (3, 1), (0, 1), (2, 1), (2, 1), (2, -1), (1, 1)]
```

A state $s$ encodes coordinate $s\mathbin{//}2$ and sign
$1-2(s\bmod2)$. At each residue the four-vector has exactly one signed
nonzero. This is the entrywise condition that will later make the endpoint
flat. Its four coordinate sums are $(2,1,1,1)$, whose squared norm is seven.
That necessary condition alone does not prove the correlation identity.

The independent verifier reconstructs the four length-seven coordinate
sequences and computes all seven periodic correlations:

```{doctest}
>>> from experiments.scalar_lift.unrestricted_parity_map import periodic_chart
>>> from experiments.scalar_lift.unrestricted_t_sign_model import verify_states
>>> data = periodic_chart(7)
>>> proof = verify_states(states, data, (2, 1, 1, 1))
>>> proof['sequence_row_sums']
[2, 1, 1, 1]
>>> proof['all_integer_correlations']
[7, 0, 0, 0, 0, 0, 0]
```

The explicit list establishes the summed circulant Gram identity
$\sum_jT_jT_j^{\mathsf T}=7I_7$. It checks integer equalities; no numerical
tolerance or search objective is involved.

## Build the disjoint design and supply a verified source

The block construction in {ref}`construction-design` produces four
coefficient matrices $X_j$ of order $28$. Their polynomial identity is

```{math}
\left(\sum_jx_jX_j\right)\left(\sum_jx_jX_j\right)^{\mathsf T}
=7\sum_jx_j^2 I_{28}.
```

The code checks the identity's diagonal and paired cross coefficients.
Separately, it checks that the four supports partition every entry.

```{doctest}
>>> from experiments.scalar_lift.hadamard_frame_search import expand
>>> from experiments.scalar_lift.williamson_scalar_endpoints import (
...     ScalarEndpoint, verify_orthogonal_design, verified_williamson)
>>> X, design_proof = verify_orthogonal_design(expand(states, data))
>>> X.shape
(4, 28, 28)
>>> design_proof['full_polynomial_Gram_verified']
True
>>> A, source_proof = verified_williamson(3)
>>> A.shape
(4, 3, 3)
```

The four order-three source matrices are symmetric commuting signs with
$\sum_jA_j^2=12I_3$. Their literal rows and complete source checks come
from the retained input table. Substitution gives
$H=\sum_jX_j\otimes A_j$, of order $28\cdot3=84$. The paired design
identities cancel the off-diagonal terms, so its Gram scale is also 84.
Disjointness ensures that each output entry is one sign times one source
sign.

## Check the actual dense endpoint independently

```{doctest}
>>> import numpy as np
>>> H = ScalarEndpoint(X, A)
>>> H.order
84
>>> dense = np.stack([H.row(i) for i in range(H.order)]).astype(np.int64)
>>> np.unique(dense).tolist()
[-1, 1]
>>> np.array_equal(dense @ dense.T, 84 * np.eye(84, dtype=np.int64))
True
>>> H.entry(12, 17) == int(dense[12, 17])
True
```

The explicit `int64` conversion keeps the Gram calculation in a sufficiently
wide integer type. A sign array stored in `int8` must not be multiplied at
large order without widening: its products can overflow even though every
individual entry is valid.

`ScalarEndpoint` can compute an entry or row without constructing the dense
matrix. This matters for the high-order outputs and families; the exact
source and polynomial proofs establish their full Grams without storing all
entries. Dense verification here is an additional independent control.

## Replay the standalone script

The complete example is an executable file, with all imports and the literal
input visible:

```{literalinclude} ../examples/verify_small_endpoint.py
:language: python
:linenos:
```

```sh
python docs/examples/verify_small_endpoint.py
```

```{doctest}
>>> import runpy
>>> example = runpy.run_path('docs/examples/verify_small_endpoint.py')
>>> result = example['verify_control']()
>>> result['order'], result['all_entries_are_signs'], result['full_integer_Gram_verified']
(84, True, True)
```

The analogous length-631 input gives coefficient order 2524. Using the same
width-three source gives order 7572. Its large saved representative and its
independent dense audit are described in {ref}`verification`.
