# Necessary gates for signed monomial multipliers of periodic T-designs

These statements concern one-hot periodic T-designs with a fixed channel
at the origin and an x-independent signed channel-permutation action.
They do not exclude unrestricted T-designs or Hadamard matrices.
The [finite control certificate](../certificates/scalar-lift/monomial-multiplier-necessary-gates.json)
records all cases and DC tests below. No active search engine was changed.

## Inversion lemma

Let `p>1` be odd and `b:Z_p->{+-e_0,...,+-e_3}` satisfy

\[
b(0)=e_0,\qquad
\sum_x b(x)^T b(x+s)=p\delta_{s,0}.
\]

Suppose a multiplier group acts by `b(gx)=M_g b(x)`, where `M_g` is a
signed permutation representation. If the group contains `-1`, the
necessary condition is

\[
p=1+2s^2\quad\text{for an integer }s.
\]

Proof. Reduce supports modulo two in
`F2[z]/(z^p-1)`. Put `B_i=sum_{x:b_i(x)!=0} z^x`, `J=sum_x z^x`, and
write `*` for `z->z^-1`. The one-hot and correlation identities give

\[
\sum_i B_i=J,\qquad \sum_i B_iB_i^*=1.
\]

The color permutation of `M_-1` fixes channel 0 and has order at most two.
If it is the identity, then

\[
1=\sum_i B_i^2=(\sum_i B_i)^2=J^2=J\ne1,
\]

a contradiction. The other possibility is a transposition of two tail
channels. Relabel it as `1<->2`, fixing 0 and 3. Its paired terms cancel
in characteristic two, leaving `(B_0+B_3)^2=1`. Squaring is an
automorphism because `f(z)^2=f(z^2)` and multiplication by 2 permutes
`Z_p`. Hence `B_0+B_3=1`. Disjoint supports and the origin state imply
`B_0=1`, `B_3=0`: both fixed channels vanish off the origin.

The signed involution exchanges channels 1 and 2 with equal signs on its
two nonzero matrix entries, so their integer row sums have equal absolute
value. The row sums are therefore `(1,s,+-s,0)`. Summing the integer
correlation identity over all shifts gives `sum_i r_i^2=p`, proving the
necessary condition. This is not a sufficiency statement.

For prime `p`, every even-order subgroup of `F_p^*` contains `-1`.
The values `(p-1)/2` at 167, 179, 191, 223, 283, 631 and 919 are,
respectively, `83,89,95,111,141,315,459`, all nonsquares. Thus no such
even-order monomial multiplier is possible. At 19 the value is 9, so this
gate is inconclusive; no existence assertion is made at length 19.

## Odd multipliers at 167, 179 and 191

The signed channel-permutation stabilizer of `e_0` has order `2^3*3!=48`.
For an odd multiplier group order `m` not divisible by 3, its image is
trivial. Each sequence is constant on the nonzero multiplier orbits, so

\[
r_0\equiv1\pmod m,\quad r_1,r_2,r_3\equiv0\pmod m,
\quad r_0^2+r_1^2+r_2^2+r_3^2=p.
\]

At 167 and 179 the only nontrivial odd group orders are 83 and 89; their
tail row sums must be zero by `|r_i|<=sqrt(p)`, giving an impossible square
value for the prime `p`. At 191 the odd orders are 5, 19 and 95. The last
two fail in the same way. For `m=5`, the first row sum is in
`{-9,-4,1,6,11}`; reducing the norm equation modulo 25 forces `r_0=-4`.
The remaining scaled squares must sum to 7, but each is in `{0,1,4}`, and
three such squares cannot sum to 7.

Combining the odd and even gates excludes every nontrivial residue
multiplier group with the stated monomial action at 167, 179 and 191.
This obstruction is broader than the three-cycle channel chart, while
remaining a restriction on symmetry rather than on general T-designs.

## Complete small controls

Every signed inversion action fixing `e_0` was checked for the identity
color action and for one tail transposition representative, up to tail
relabeling. The enumeration freely assigns any of eight signed channels
to one position in each pair `{x,-x}` and propagates the other; fixed
channels were not excluded in advance. Every nonzero integer correlation
was tested until rejection or acceptance.

| Length | Identity action-word cases | Transposition action-word cases | Distinct accepted transposition words |
|---:|---:|---:|---:|
| 3 | 64 | 32 | 8 |
| 5 | 512 | 256 | 0 |
| 7 | 4096 | 2048 | 0 |
| 9 | 32768 | 16384 | 0 |

There are eight signed identity actions and four signed actions for the
chosen transposition. Unused-channel signs cause duplicate accepted words;
the certificate distinguishes action-word pairs from distinct words.
Length 3 supplies the positive witness `(e_0,e_1,e_2)`. Length 9 passes
the square gate but has no word in this inversion-equivariant control,
illustrating that the necessary condition is not sufficient.

The enumeration can be replayed directly with the following code:

```python
from itertools import product

for n in (3, 5, 7, 9):
    for name, perm in (("identity", (0, 1, 2, 3)),
                       ("transposition", (0, 2, 1, 3))):
        signs = ([(1, *s) for s in product((1, -1), repeat=3)]
                 if name == "identity" else
                 [(1, e, e, f) for e in (1, -1) for f in (1, -1)])
        cases, accepted = 0, set()
        for eps in signs:
            for half in product(range(8), repeat=(n - 1) // 2):
                w = [0] * n  # state=2*channel+negative_bit
                for x, state in enumerate(half, 1):
                    axis = state // 2
                    s = (1 - 2 * (state % 2)) * eps[axis]
                    w[x], w[-x] = state, 2 * perm[axis] + int(s < 0)
                cases += 1
                def corr(shift):
                    return sum(1 if ((w[x] ^ w[(x + shift) % n]) & 1) == 0 else -1
                               for x in range(n)
                               if w[x] // 2 == w[(x + shift) % n] // 2)
                if all(corr(s) == 0 for s in range(1, n)):
                    accepted.add(tuple(w))
        print(n, name, cases, len(accepted))
```

## Related gate for a local four-row clock

For the fixed skew conference matrix

\[
D_4=\begin{pmatrix}0&1&1&1\\-1&0&1&-1\\-1&-1&0&1\\-1&1&-1&0\end{pmatrix},
\qquad R_6=(I-D_4)/2,
\]

the orthogonal sixth-root action sends a sign vector `x` to a sign vector
exactly when `product_i x_i=+1`. Indeed, every row of `I-D_4` has product
`-1`, so its dot product with such an `x` is `+-2`. The other eight cube
vectors map to vectors with one entry `+-2` and three zeros. All 16 cases
are checked in the certificate; the ordinary H4 supplies a positive frame.

Applying block-diagonal copies of this local action to a real Hadamard
matrix would require a closed quadruple of rows: their entrywise product
must be constant. No closed quadruple exists at order `N=4p`, with `p>1`
odd. Normalize one of the four rows to ones. The others are `u,v,+-uv`,
and orthogonality makes their four column classes all have size `p`.
A fifth row orthogonal to those four must have sum zero within every
class, by the invertible four-by-four Walsh system. A sum of `p` odd
signs cannot be zero, a contradiction.

Thus known H668, H716 or H892 cannot automatically supply a flat pair by
this local four-row action. This excludes that shortcut, including local
signed row/column variants; other rational sixth-root actions and globally
mixed compatible frames remain outside the conclusion.
