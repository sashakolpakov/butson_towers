# An obstruction to even signed-permutation multipliers at prime length

This note proves an obstruction used in the [multiplier search selection](HADAMARD_FRAME_TARGET_SCAN.md).

**Theorem.** Let p>3 be an odd prime, and let

\[
f:\mathbb F_p\longrightarrow\{\pm e_0,\pm e_1,\pm e_2,\pm e_3\},
\qquad f(0)=e_0,
\]

where e0,…,e3 are the standard basis vectors of R^4. Suppose its periodic correlations satisfy

\[
C(t)=\sum_{x\in\mathbb F_p} f(x)^{\mathsf T}f(x+t)
=p\,\delta_{t,0}.
\tag{1}
\]

There is no signed permutation involution J such that

\[
f(-x)=Jf(x)\quad\text{for every }x\in\mathbb F_p.
\tag{2}
\]

Consequently, there is no even-order subgroup H of F_p^* with a signed permutation representation M:H→O(4) and covariance f(hx)=M_h f(x). Indeed, −1 belongs to H and J=M_{−1} would satisfy (2) and J²=I. Equivalently, this covers the usual model f(μx)=Cf(x), where μ has even order m, C is a signed permutation, and C^m=I. A generator fixing f(0)=e0 fixes that vector with positive sign.

The theorem concerns this exact signed-permutation covariance at odd **prime** length. It gives no assertion about composite lengths, general linear coordinate actions, or periodic T-functions without this symmetry.

## 1. Pairing a nonzero correlation

Fix t≠0. The map x↦−x−t is an involution of F_p with the unique fixed index x=−t/2. The two corresponding correlation summands are equal, because J is orthogonal:

\[
\begin{aligned}
f(-x-t)^{\mathsf T}f(-x)
&=(Jf(x+t))^{\mathsf T}(Jf(x))\\
&=f(x+t)^{\mathsf T}f(x).
\end{aligned}
\]

All summands are integers. Thus, modulo2, only the fixed index remains:

\[
0=C(t)\equiv f(-t/2)^{\mathsf T}f(t/2)
=f(t/2)^{\mathsf T}Jf(t/2)\pmod2.
\tag{3}
\]

For a signed basis vector v=±e_i, the integer v^T Jv is ±1 if the underlying coordinate permutation of J fixes i, and0 otherwise. Equation(3), as t/2 runs over every nonzero residue, therefore says that every coordinate used away from0 must be moved by that permutation.

Equation(2) at0 gives Je0=e0. Since J²=I, its underlying permutation fixes coordinate0 and is an involution on the remaining three coordinates. It cannot be the identity, because then no coordinate could occur at a nonzero residue. After relabelling coordinates1,2,3, its only possible moved coordinates are1 and2. The involution condition gives a common sign ε∈{±1} with

\[
Je_1=\varepsilon e_2,\qquad Je_2=\varepsilon e_1.
\]

Coordinate3 is fixed, possibly with negative sign, and cannot occur away from0. Hence f0 is the indicator of0, f3 is identically0, and the remaining two coordinate functions u=f1 and v=f2 obey

\[
u(0)=v(0)=0,\quad v(x)=\varepsilon u(-x),\quad
u(x)^2+u(-x)^2=1\quad(x\ne0).
\tag{4}
\]

The last identity states that the supports of u and its reversal are disjoint and partition F_p^*.

## 2. A circulant weighing matrix

Let W be the circulant matrix W_{xy}=u(y−x), with rows and columns indexed by F_p. The four circulant matrices associated with f are

\[
I,\quad W,\quad\varepsilon W^{\mathsf T},\quad0.
\]

The correlation identities(1) are exactly the sum-of-Grams identity for those four matrices. Since circulant matrices commute, it follows that

\[
I+WW^{\mathsf T}+W^{\mathsf T}W=pI,
\qquad WW^{\mathsf T}=\frac{p-1}{2}I.
\tag{5}
\]

Furthermore W has entries in {0,±1}, zero diagonal, and W∘W^T=0, where ∘ denotes entrywise multiplication. Thus its support has the additional disjoint-reversal property(4).

Every row and column of W has the same integer sum r. Applying(5) to the all-ones vector gives r²=(p−1)/2. Write s=|r|≥1. Necessarily

\[
p=2s^2+1,\qquad WW^{\mathsf T}=s^2I.
\tag{6}
\]

The following obstruction needs only the zero coefficient u(0)=0, not the full support condition(4).

## 3. A self-conjugacy lemma

**Lemma.** Suppose a nonzero circulant matrix W with entries in {0,±1} of prime order p has zero diagonal and WW^T=s²I. If a prime q≠p divides s and has even multiplicative order modulo p, then W cannot exist.

**Proof.** Let

\[
a(z)=\sum_{x=0}^{p-1}u(x)z^x,\qquad
\Phi_p(z)=1+z+\cdots+z^{p-1}.
\]

The Gram identity implies a(z)a(z^{-1})=s² modulo z^p−1. Reduce coefficients modulo q and then pass to

\[
R=\mathbb F_q[z]/(\Phi_p(z)).
\]

This ring is reduced: q≠p makes z^p−1 squarefree, so its factor Φ_p is squarefree as well. Equivalently, R is a product of finite fields.

Even multiplicative order gives an integer j≥1 with q^j≡−1 modulo p. In R, Frobenius therefore acts by reversal:

\[
a(z)^{q^j}=a(z^{q^j})=a(z^{-1}).
\]

Since q divides s, the Gram identity becomes a(z)^{1+q^j}=0 in R. A reduced ring has no nonzero nilpotent elements, so a(z)=0 in R. Thus Φ_p divides a in F_q[z]. The degree bound deg(a)≤p−1 implies a=cΦ_p for a scalar c∈F_q. Its constant coefficient is u(0)=0, giving c=0. This contradicts the presence of a nonzero coefficient±1. ∎

In(6), the hypothesis q≠p is automatic because p>s.

## 4. Every case p>3 has a forbidden prime divisor

It remains to show that(6) always supplies a prime q as in the lemma.

If s is even, suppose every prime divisor q of s had odd multiplicative order modulo p. Products of commuting elements of odd order also have odd order. Since2 divides s, both2 and s would consequently have odd order in F_p^*. Then2s² would have odd order. Equation(6) gives2s²=−1 modulo p, a contradiction because −1 has order2. Thus some prime divisor of s has even order, and the lemma applies.

If s is odd, primality of p=2s²+1>3 forces3 to divide s: otherwise p would be divisible by3. Now s²≡9 modulo12, so p≡7 modulo12. The standard supplementary form of quadratic reciprocity gives

\[
\left(\frac3p\right)=-1
\qquad(p\equiv7\pmod{12}).
\]

By Euler's criterion, 3^{(p−1)/2}=−1 modulo p. In particular,3 has even multiplicative order modulo p. Since3 divides s, the lemma applies again. This proves the theorem.

## 5. The prime3 exception

At p=3, take f(0)=e0, f(1)=e1 and f(2)=e2. Let J fix e0,e3 and exchange e1,e2 with positive signs. Then f(−x)=Jf(x), and every summand in either nonzero correlation pairs two different coordinate vectors, so C(1)=C(2)=0. Here s=1 and W is a permutation circulant of weight1. There is no prime divisor of s to which the lemma could apply.

Thus p=3 is a genuine exception, not an omitted boundary case.
