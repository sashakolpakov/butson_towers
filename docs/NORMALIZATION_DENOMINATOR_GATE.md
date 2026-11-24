# Necessary denominator primes for normalized sign images

Let `F,G` be sign Hadamard matrices of order `d=4r`, where `r>1` is odd.
The unique rational operator mapping `F` to `G` is

`R = G F^T / d`.

**If no odd prime dividing r occurs in a reduced entry denominator of R,
then R is a signed permutation.** In particular, a dyadic orthogonal
operator cannot supply a nontrivial cyclotomic sign image at these widths.
The claim concerns the current uniform orthogonal scalar-action model,
where each individual used block has Gram `dI`; it does not exclude
arbitrary weighted restrictions or nonuniform changes across source
coordinates.

## Proof from row-product parity

Within a Hadamard matrix of order `4r`, every row has the same sign
product: two rows have Hamming distance `d/2=2r`, which is even.
Therefore all cross inner products between rows of G and rows of F are
congruent modulo four.

If the two common row-product signs agree, every cross product is
divisible by four, so every entry of R lies in `(1/r) Z`. Removing all
denominator primes of r makes R integral. An integral orthogonal matrix
is a signed permutation.

If the row-product signs differ, every cross product is two modulo four.
Every entry of R is consequently `k/(2r)` with k odd. If no prime of r
survives in its reduced denominator, every entry is an odd half-integer.
Each squared row norm is at least `d/4=r>1`, contradicting orthogonality.
This also proves that all nonmonomial maps require a denominator prime
of r. Denominators involving unrelated primes such as three or five do
not help at r=167,179 or191.

## A cyclotomic field action cannot be the surviving permutation

An actual field representation of a primitive m-th root satisfies
`Phi_m(R)=0`; every eigenvalue is a primitive m-th root. A signed cycle
has polynomial `x^ell-epsilon`. If epsilon is positive it has eigenvalue
one. If epsilon is negative, write `ell=2^a b` with b odd; then
`x^(2^a)+1` divides the cycle polynomial, providing a primitive
`2^(a+1)`-th root eigenvalue. Every signed permutation therefore has an
eigenvalue whose order is a power of two.

It follows that a cyclotomic field action of a root with nontrivial odd
order part cannot be this signed permutation. A faithful cyclic **group**
action is insufficient for this argument: a three-cycle is faithful but
also has eigenvalue one and does not satisfy `Phi_3`.

The used alphabet need not contain one. If two used roots u,v have a
ratio with nontrivial odd order part, their normalized sign images
`rho(u)F` and `rho(v)F` are both Hadamard matrices. Their ratio operator
is `rho(u/v)`, so the theorem applies to it. Only an actually used ratio
is needed; an unused odd part in the nominal root label does not suffice.

## Stronger sixth-root consequence

For a rational orthogonal R satisfying `R²-R+I=0`, orthogonality gives
`R+R^T=I`, hence each diagonal entry is one half. Its cross diagonal is
`d/2=2r`, which is two modulo four. This forces the opposite-parity case
above for every cross entry. **Every entry of R is odd/(2r), so R is
fully dense. If r is prime, its exact common reduced denominator is 2r.**

This gives an inexpensive precondition before applying a candidate
operator to any large frame: `dR` must be integral, and in this sixth-root
case every entry of `dR` must be two modulo four.

For example, conjugating a dyadic starting action by a dyadic Walsh
orthogonal basis keeps the action dyadic and is excluded universally.
Useful rational basis candidates must introduce an odd prime of r.
A generic basis with denominator r may leave an action denominator
`2r²`; the required cancellation can be checked before any frame product.
This is a necessary condition, not an existence result.

## Exact positive controls and API

The [replayable certificate](../certificates/scalar-lift/normalization-denominator-gate.json)
checks three nonmonomial sixth-root sign-pair controls:

* Width four, denominator two: the condition `r>1` is necessary.
* Width eight, denominator two: the condition `d=4 mod8` is necessary.
* Width twelve, denominators two and six: an odd denominator prime of r
  genuinely permits a nontrivial pair. This is the saved exact order-twelve
  pair from `inverse_blocks.known_order_twelve_pair`, with entries of R
  equal to `±1/2,±1/6`.

All complete sign Grams, sixth-root polynomials, orthogonality identities
and exact denominators are checked. The general gate is
`denominator_gate(width, localization_denominator, primitive_ratio_order)`
in [the implementation](../experiments/scalar_lift/normalization_denominator_gate.py).
Its localization denominator must include every prime inverted in the
actual operator, including those of the inverse basis if needed.

```sh
PYTHONPATH=src:. ./.native-env/bin/python -m experiments.scalar_lift.normalization_denominator_gate --check certificates/scalar-lift/normalization-denominator-gate.json
```
