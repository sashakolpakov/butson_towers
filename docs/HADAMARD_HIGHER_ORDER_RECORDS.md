# Higher-order records and repeat-target audit

Completed-order ledger checked on 24 November 2026; target-status update
on 26 November 2026. The completed new real Hadamard orders are
**7572, 17668, 22716, 32812, 42908, 58052, 73196, 78244, 93388, 98436,
108532, 128724, 153964, 229684, 295308, 386172, 4772884 and 13001124**.
The periodic T-matrices of order **127** are a new ingredient.

The main higher-order campaign has **38 completed real Hadamard orders**:
**18 new** as listed above and **20 known recovered orders**. The counts
are unions of distinct orders from the exact endpoint certificates;
repeated constructions of the same order are counted once. The new
periodic T-matrix order 127 is separate; its known real endpoint H1524
is already included in the real-order count.

**Status convention:** “New” means new to the best of the author's knowledge
after the recorded literature checks. “Known” means that a prior construction
has been located. These labels apply to completed constructions; an
incomplete search remains a target.

The dated machine-readable record is
[higher-order-records.json](../certificates/scalar-lift/higher-order-records.json).
The dated evidence includes explicit unresolved entries, positive table
entries and applicable construction theorems with their hypotheses checked.
Software `Unknown` or `NULL` records an implementation gap. The public
records located here do not provide a comprehensive order-status table
extending to 100,000.

## Completed real endpoints

| Real endpoint | Complete source | Gauged source exponents modulo 6 | Status |
|---|---|---|---|
| 7572 = 2524 × 3 | `unreal3-q6`, BH(3,6) | 1,3 | New |
| 17668 = 2524 × 7 | `catalogue-p7-q6`, BH(7,6) | 0,1,3,4 | New |
| 22716 = 2524 × 9 | `aalto-9-6-000`, BH(9,6) | 0,1,3,4 | New |
| 32812 = 2524 × 13 | `unreal13-q6`, BH(13,6) | 0,1,3,4 | New |
| 98436 = 2524 × 39 | Tripled unreal order-13 source, BH(39,6) | 0,1,3,4 | New |
| 295308 = 2524 × 117 | Twice-tripled unreal order-13 source, BH(117,6) | 0,1,3,4 | New |
| 885924 = 2524 × 351 | Three-times-tripled unreal order-13 source, BH(351,6) | 0,1,3,4 | Known: Paley I, prime 885923 |

The [reuse certificate](../certificates/scalar-lift/hadamard-7572-frame-reuse.json)
records literal gauges, source hashes, complete exact source Grams, and the
already certified frame identities. All use `Q(zeta_6)/Q`, its degree-two
action `R=QP^T/2524`, transpose conjugation, and trace metric
`[[2,1],[1,2]]` in the canonical basis. The actual rational representation has
dimension 2524 and field-module rank 1262. Its frame `P` has Gram scale 2524;
the four used scalar images are `P,Q,-P,-Q`. Multiplying the source Gram
scale by 2524 gives each endpoint's order and Gram scale after the sign-entry
test. The [tripling executor and proof](TRIPLED_UNREAL_SCALAR_ENDPOINTS.md)
complete the sources of orders 39,117,351 and their endpoint recipes.

Direct Paley I, Paley II and Yamada parameter checks do not settle these three
orders. Their candidate parameters factor as follows:

| Endpoint N | N−1 | N/2−1 | N/4−2 |
|---|---|---|---|
| 17668 | 3²·13·151 | 11²·73 | 5·883 |
| 22716 | 5·7·11·59 | 41·277 | 7·811 |
| 32812 | 3·10937 | 5·17·193 | 59·139 |

Each entry in this table has at least two distinct prime factors. Thus the
prime-power hypotheses of these particular direct families fail. This is
not an exclusion of other constructions. The hypotheses used are stated in
[Egan–Kotsireas–London–Ó Catháin, v2, Theorem 2 and Table 4](https://arxiv.org/html/2610.27175v2#S4.T4)
and [Đoković 2011, Theorems 2.2–2.3, p.3063](https://www.m-hikari.com/imf-2011/61-64-2011/djokovicIMF61-64-2011.pdf#page=3).

## Williamson scalar endpoints and constructive families

The complete T631 orthogonal design has weighted scalar sources over
`Q[C_w]`, explicitly a product of cyclotomic fields. The normalized regular
trace gives metric `I_w` in the inverse group basis. The complete source
norm is `2524w I`; all actual multiplication images are sign circulants,
giving endpoints of order and Gram scale `2524w`. The
[27-source certificate](../certificates/scalar-lift/williamson-scalar-endpoints.json)
and [input and prior-construction record](WILLIAMSON_ORDER_RECORDS.md)
add these eight new orders:

| Williamson source order | New real order |
|---:|---:|
| 17 | 42908 |
| 23 | 58052 |
| 29 | 73196 |
| 31 | 78244 |
| 37 | 93388 |
| 43 | 108532 |
| 51 | 128724 |
| 61 | 153964 |

The other 19 source cases reproduce five previously completed new orders
and 14 known orders. Each known order has a positive prior construction
with exact parameters in the
[Williamson order-status data](../experiments/scalar_lift/williamson_order_status.json).
The [T127/Williamson-17 endpoint H8636](../certificates/scalar-lift/williamson-scalar-known-control.json)
is also known: prior Williamson-type matrices of order 127 and classical
T17 already give it by Cooper–Wallis. Williamson-type 127 and periodic
T127 are different ingredients.

[Spence's 1977 recurrence, Theorem 4](https://doi.org/10.1017/S1446788700020255),
with five exactly verified starters, supplies the following
[constructive infinite scalar families](SPENCE_SCALAR_FAMILIES.md).
Existence and literal entries hold for every integer `r>=0`.

| Parameter q | Endpoint order for r>=0 | First expanded member r=1 | Individual status |
|---:|---|---:|---|
| 13 | `17668*13^r` | 229684 | New |
| 17 | `22716*17^r` | 386172 | New |
| 61 | `78244*61^r` | 4772884 | New |
| 101 | `128724*101^r` | 13001124 | New |
| 121 | `153964*121^r` | 18629644 | Known: Paley II, prime 9314821 |

The [source induction and exact finite source audits](../certificates/scalar-lift/williamson-spence-family.json)
and [five completed endpoint proofs](../certificates/scalar-lift/williamson-family-endpoints.json)
retain group labels, field factors, conjugation, trace metric, actual
polynomial alphabets, source and endpoint scales, and literal entry access.
The [family order-status data](../experiments/scalar_lift/williamson_family_order_status.json)
classify these five finite members. The infinite existence statements do
not label every member's order new.

Together with the earlier known orders 1524, 3252, 6564 and 885924, the
14 known T631/Williamson endpoints, H8636 and H18629644 give the 20 known
recovered orders. Small implementation controls and older scalar inventory
outputs remain outside this main-campaign count.

## What the public catalogs actually cover

| Primary record | Version or date | Coverage relevant here | Consequence |
|---|---|---|---|
| [Egan–Kotsireas–London–Ó Catháin survey](https://arxiv.org/html/2610.27175v2) | v2, 3 November 2026 | Unresolved inventory through order 10,000; general construction theorems extend further | Its inventory does not classify the three endpoints |
| [Cati–Pasechnik](https://arxiv.org/html/2411.18897v2) | v2, 30 August 2025 | Table 4 covers odd part n≤2999, hence orders 4n≤11996; m=2 entries are omitted by an explicit convention | Odd parts 4417,5679,8203 lie outside its table |
| [Kumar's author-maintained frontier](https://github.com/vatskishan03/Hadamard/blob/d885022144ba6bf4efa8315e14871a3736947527/FRONTIER.md) | commit d885022, 6 November 2026 | Updates Cati–Pasechnik's same odd-part range | Does not extend the table to the three orders |
| [Handbook, second edition, Table 1.53, as described by Đoković](https://arxiv.org/html/0912.5091) | Book published November 2006, copyright 2007; table pp.278–279 | All odd n<10000, so 4n<40000 | Its range includes the targets, but their literal entries have not been obtained in this audit; see the lookup below |
| [Seberry–Yamada 1992 appendix](https://documents.uow.edu.au/~jennie/WEBPDF/175_1992.pdf#page=105) | 1992, printed pp.535,547–554 | Says an electronic table covered odd q<40000, but prints only q<3000 | The larger electronic table is not contained in this PDF |
| [Seberry–Yamada M-structures author PDF](https://documents.uow.edu.au/~jennie/WEBPDF/134_1989.pdf) | Author preprint; published 1990 | Announces an index for q<40000, but this 28-page file ends its displayed Hadamard index at q=1999 | The advertised bound does not supply the missing entries |
| [Seberry–Yamada book](https://doi.org/10.1002/9781119520252) | 2020, Appendix A.1 scope statement and A.3 | Updated ordinary Hadamard table covers odd q≤3000 | Does not supply the three higher odd parts |
| [SageMath Hadamard module](https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/matrices/hadamard_matrix.html) | Reference Manual 10.10; source develop HEAD 8de8dad7, 14 November 2026 | Arbitrary order argument; complete implementation claim only through 1200, partial above it | A positive construction query would suffice; `Unknown` would not classify novelty |
| [Magma Hadamard database](https://docs.magma-maths.org/Combinatorics/HadamardMatrices/hadamard-database.html) | Current handbook | Matrices of degrees through 256; many classes at small degrees | Not a table of higher admissible orders |
| [HadamardR](https://cran.r-project.org/web/packages/HadamardR/refman/HadamardR.html) | 1.0.0, 7 April 2020 | Reports package coverage below order 5000 and provides construction predicates | No exhaustive larger-order status claim |
| [hoa64 author repository](https://github.com/bbeartheancient/hoa64) | master; last pushed 16 August 2026 | Reports 810 verified orders, largest 7408; complete CSV collection is not in the public repository | “810 orders” does not mean coverage through 100,000 |

The Sage source was inspected, including its construction dispatch and
existence semantics. Two public SageCell requests failed (HTTP 520 and a
20-second timeout), so this record attributes no target-specific `True` or
`Unknown` response to Sage. Private snapshots, SHA-256 digests and exact URLs
are listed in the JSON; full third-party PDFs are not included in the clean
certificate directory.

The [Handbook editor's update page](https://site.uvm.edu/jdinitz/?page_id=404)
was also checked. It was last edited in February 2018 and supplies an update
at order 764, not a current replacement table for the three orders. Historical
references to 100,000 matrices found by MathCheck refer to numbers of objects
at smaller orders, rather than an order-status catalog through 100,000.

## Historical Table 1.53 and its published corrections

The Handbook's larger historical range is explicit in [Đoković's December
2009 paper, introduction and §3](https://arxiv.org/html/0912.5091). It describes
5,000 entries for odd `n<10000`, with a dot meaning that `H(4n)` was known,
and an exponent at least 3 in 1,006 other cases. The paper gives 138
corrections and four further constructions, reducing that count to 864.
[The August 2010 follow-up, introduction and §2](https://arxiv.org/html/1008.2043)
gives 42 disjoint further changes, leaving 822. These are historical reported
counts, rather than a current unresolved-order inventory.

The explicit lists in both papers were read and their target membership
checked exactly:

| Real order 4n | Odd part n | Literal 2007 table entry read | In the 142 changes of 2009 | In the 42 changes of 2010 |
|---|---|---|---|---|
| 17668 | 4417 | No | No | No |
| 22716 | 5679 | No | No | No |
| 32812 | 8203 | No | No | No |
| 25732 | 6433 | No | No | No |
| 14292 | 3573 | No | No | No |

These are historical lookup results. The 26 November positive construction
below settles H14292 independently of the unread Handbook entry.

The original entries were not obtained. Đoković also reports in the 2010
introduction that some positive Handbook entries lacked sufficient
construction references for independent verification. The JSON preserves
the exact correction-list membership checks and the missing baseline fields.

The [official publisher preview](https://api.pageplace.de/preview/DT0400.9781420010541_A25103714/preview-9781420010541_A25103714.pdf)
has 102 PDF pages and omits printed pp.278–279. The
[institutional author-preparation record inria-00429215](https://www-sop.inria.fr/members/David.Coudert/Biblio/Year/2006.complete.html)
is Bermond–Coudert's separate chapter “Grooming,” pp.494–496, despite its
generic Handbook cover. It does not contain Table 1.53. Kotsireas's 2025
conference slides were indexed but could not be fetched with valid TLS; an
[accessible university-hosted 2024 presentation](https://fmi.uni-sofia.bg/sites/default/files/wmi2024/presentations/a260900_presentation.pdf)
does not reproduce the table.

The larger-order entries in Kumar's repository do not amount to a complete
transcription of Table 1.53. At commit `d885022`, its
[novelty policy](https://github.com/vatskishan03/Hadamard/blob/d885022144ba6bf4efa8315e14871a3736947527/NOVELTY.md)
says that the Handbook was consulted, and its
[larger-order list](https://github.com/vatskishan03/Hadamard/blob/d885022144ba6bf4efa8315e14871a3736947527/hadamard/README.md)
records individual completed constructions with scoped labels. The shared
ordinary-table data file is the 195-row Cati–Pasechnik transcription through
odd part 2999. The complete public tree and the relevant data and prose were
inspected; no full Handbook table was obtained.

The 2020 Seberry–Yamada book's appendix introduction restricts its updated
ordinary table to `q≤3000`; it therefore does not close the three lookup gaps.
Its §3.6, pp.87–90, gives a conditional T-property construction and computed
examples through order 79. This scope was checked against the original book
text in a [publicly indexed access copy](https://dokumen.pub/hadamard-matrices-constructions-using-number-theory-and-algebra-9781119520245-9781119520276-9781119520139.html),
with the [publisher DOI](https://doi.org/10.1002/9781119520252) identifying the
underlying primary work. It is not an unconditional construction for all
prime powers of the form `a²+3b²`.

## A decisive pruning rule for repeats

Turyn's family gives Williamson matrices of order `(q+1)/2` for prime powers
`q ≡ 1 (mod 4)`. In particular, if `2p−1` is such a prime power, Williamson
matrices of order `p` are already known. The original theorem is stated in
[Whiteman 1973, abstract](https://www.sciencedirect.com/science/article/pii/0097316573900101),
which gives a new proof of Turyn's theorem; it is also stated explicitly in
[the 2026 survey, §2.1 and Table 3](https://arxiv.org/html/2610.27175v2#S2.SS1).

Cooper–Wallis then combines these with known T-matrices of orders
`w=3,7,9,13` to give **every order `4pw`** in those four source columns.
See [Cati–Pasechnik, Theorem 5](https://arxiv.org/html/2411.18897v2#S4.SS4)
and the [original Cooper–Wallis paper](https://doi.org/10.1017/S0004972700045019).
The small T inputs are classical; the original paper explicitly covers the
odd parameters through 19.

| Frame prime p | 2p−1 | Consequence for w=3,7,9,13 |
|---|---|---|
| 271 | 541 prime | All four endpoint orders already known |
| 307 | 613 prime | All four already known, including 15964 |
| 331 | 661 prime | All four already known, including 17212 |
| 547 | 1093 prime | All four already known |

This rules out H17212 as a new-order target even though its direct Paley
parameters fail. It also explains why a successful independent-prime repeat
can be a useful method control while adding no new Hadamard order.

**Positive p397 exclusion, 26 November 2026.** Paley's construction gives
a symmetric conference matrix of order `198=197+1`, since 197 is prime
and congruent to 1 modulo 4. Seberry–Yamada,
[*On the products of Hadamard matrices, Williamson matrices and other orthogonal matrices using M-structures*](https://documents.uow.edu.au/~jennie/WEBPDF/139_1990.pdf#page=24),
J. Combin. Math. Combin. Comput. 7 (1990), 97–137, Corollary 29(ii),
printed p.120, therefore gives symmetric Williamson-type matrices of
order `397=2*198+1`; that page also explicitly lists 397. This uses part
(ii), separately from the missing-hypothesis issue in part (i).
Cooper–Wallis, [Theorem 3, p.274](https://doi.org/10.1017/S0004972700045019),
combines them with classical T-matrices of order 9 to give **H14292**.
Their Corollary 4 and the table on pp.274–275 supply T-matrices for every
odd order through 19. Thus every `4*397*w` with odd `w<=19` is already
known; the p397 parity map is retained only as a method control.

## Remaining targets and exact scope

Updated on 26 November 2026.
For `p=919`, the source-order-3 and -9 orders H11028 and H33084 are already
Paley I orders, with prime parameters 11027 and 33083. H47788, from source
order 13, is Paley II with prime parameter 23893. **H25732** and **H40436**,
from Williamson widths 7 and 11, remain targets for a completed length-919 input. Here
`2p−1=1837=11·167`, so the preceding Turyn family does not settle them.

The p1453 chart targets **H63932/H75556** (Williamson widths 11/13), and p1951 targets
**H23412/H54628/H101452** (widths 3/7/13). No prior construction was
located for these seven targets in the checks recorded here. The earlier
JSON records the positive Paley constructions in the other columns.
Width 11 requires the complete T-input's four-variable orthogonal design;
it is not a new admitted sixth-root source for the P/Q frame alone.
H17436 is already known: Paley gives skew H2180 because 2179 is
prime; Yamada with prime `q=4357 ≡5 (mod8)` then gives `4(q+2)=17436`.

The [p1951 parity map](../certificates/scalar-lift/hadamard-frame-p1951.json)
is now verified, with 50 nonorigin orbits and a generic normalized affine
sign-lift estimate of 23 free bits; each axis word requires exact elimination.
Its lowest prospective output is H23412; no sign witness
or endpoint has been completed. The
[expanded multiplier scan](../certificates/scalar-lift/multiplier-target-scan-20261126.json)
also queues `p3571,m51` and `p4219,m57` as arithmetic candidates, with no
verified parity maps yet. The [target scan](HADAMARD_FRAME_TARGET_SCAN.md)
separates these stages explicitly.

For targets with explicit current unresolved records, the stronger literature
position is still a length-223 frame: H2676 is listed in the
[2026 survey's Table 5](https://arxiv.org/html/2610.27175v2#S4.T5), and H11596
is an actual unresolved row of the pinned author frontier. A length-283 frame
would similarly reach H7924, explicitly listed in that survey and frontier.
These require larger sign charts than the compact prime repeats. They are
conditional endpoint targets, not completed constructions.

## A completed periodic ingredient at order 127

The repeat search produced **new periodic circulant T-matrices of order 127**. This
is exactly the object defined immediately before
[Table 3 of Egan–Kotsireas–London–Ó Catháin v2](https://arxiv.org/html/2610.27175v2#S4.T3),
dated 3 November 2026, which lists 127 as the smallest open T-matrix order.
The four generating rows partition all 127 positions, have row sums
`(10,3,3,3)` and row weights `(46,27,27,27)`, and satisfy every periodic
correlation identity. The
[witness certificate](../certificates/scalar-lift/hadamard-frame-p127.json)
and the separate independent verification give the existence proof. Its
Butson scalar reduction gives H1524, an already known order (Paley I with
prime parameter 1523).

The ingredient is a periodic T-matrix quadruple; aperiodic T-sequences and
Williamson-type matrices are separate classes. The dated
[ingredient literature audit](../certificates/scalar-lift/periodic-t127-literature-audit.json)
pins the subsequent public repository checks: Kumar's complete `main` tree
at `d885022` (6 November) and Jones's at `4213419` (12 October) contain no
order-127 T-matrix witness. Math74's public October 2026 index is an ordinary
Hadamard-order catalog, not a T-matrix registry.

The openly available
[2017 author manuscript](https://documents.uow.edu.au/~jennie/WEB/WEB17/Xtb32_6mp1.pdf)
provides conditional searches and examples through 79, and the
[2024 composite-order paper](https://toc.ui.ac.ir/article_27062_c8e585f002d7d36b07ebc0cab4ddbe16.pdf)
combines existing special T-matrix inputs; neither supplies a T127 or T631
construction in the cases checked.
