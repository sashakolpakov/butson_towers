(references)=
# References and retained source material

The mathematical starting references and the executable source inputs have
different roles. Morphism and scalar-restriction results supply general
principles; the literal T-inputs, source quadruples and certificates establish
the specific outputs described in this manual.

## Earlier morphisms and the starting reference

- Ronan Egan and Pádraig Ó Catháin,
  [*On complete morphisms of Butson Hadamard matrices*](https://arxiv.org/abs/2610.28394v1),
  20 October 2026; [*Journal of Combinatorial Designs*, Early View,
  2 November 2026](https://doi.org/10.1002/jcd.70047).
  This was the immediate starting reference for the project.
  Its complete-morphism framework treats a fixed map on an entire stated
  root alphabet. A frame passing only an admitted palette has a narrower
  domain.
- Ronan Egan and Pádraig Ó Catháin,
  [*Morphisms of Butson classes*](https://arxiv.org/abs/1707.08815v3),
  *Linear Algebra and its Applications*577 (2019), 78–93. This develops
  morphisms that change the root alphabet while expanding the matrix size.
- Pádraig Ó Catháin and Eric Swartz,
  [*Homomorphisms of matrix algebras and constructions of Butson-Hadamard matrices*](https://arxiv.org/abs/1904.10771v3),
  *Discrete Mathematics*342 (2019),111606. This supplies the repeated-prime
  exponent-reducing construction used by the retained exact executor.

The project name “Butson towers” describes successive restrictions; it makes
no priority claim for restriction of scalars or Butson morphisms. The
{repo-file}`foundations note <docs/FOUNDATIONS.md>` records the attribution
and precise theorem references used in the manuscript.

## Source quadruples and recurrences

E. Spence, [*An infinite family of Williamson matrices*](https://doi.org/10.1017/S1446788700020255),
*Journal of the Australian Mathematical Society*, Series A 24 (1977), 252–256,
Theorem 4, supplies the classical recurrence used for eight announced source families.
The literal starters, finite-field conventions and recurrence proofs appear
in {ref}`families` and the manuscript.

The small finite Williamson rows are retained in
{repo-file}`williamson_inputs.json <experiments/scalar_lift/williamson_inputs.json>`.
Their acquisition record cites the pinned
[SageMath source table](https://github.com/sagemath/sage/blob/8de8dad7bcb38e1baa43508cbe944bf04a2d0401/src/sage/combinat/matrices/hadamard_matrix.py).
The executable verifier establishes the exact sign, symmetry, commutation
and norm identities needed here; a bibliographic entry alone is not a source
certificate.

The {repo-file}`additional literal table <experiments/scalar_lift/near_williamson_inputs.json>`
retains 18 published near-Williamson inputs and the ordinary width-69
Williamson input. The source URLs, literal row choices and reversal flags
are recorded with the data. Their exact signs, norms, amicability and
operator-algebra reconstruction are checked in the separate certificates
listed in the {repo-file}`source inventory <docs/SOURCE_INVENTORY.md>`.

The prior constructions at the 31 recovered known orders have their own
table and citations in {ref}`known-orders`. This distinguishes a previously
known order from an equivalence claim about a particular representative.

## Paper, literal data and formal proofs

| Material | Purpose |
|---|---|
| {repo-file}`Manuscript source <manuscript/hadamard_7572_announcement.tex>` | Definitions, constructions, complete finite inputs and mathematical proofs |
| [Published first version](https://hexagonmath.org/2611.01070v1) | Existing HexagonMath article, with version 2 submitted as a revision |
| {repo-file}`Original submitted PDF <release/butson_towers_v3/announcement.pdf>` | Frozen first-version manuscript |
| {repo-file}`Current revision package <release/hexagon_2611_01070_v2/README.md>` | Submitted version 2 manuscript and replay records |
| {repo-file}`Version 2 receipt <release/hexagon_2611_01070_v2/submission-status.json>` | Revision of the existing article; service status under review |
| {repo-file}`Original T631 recipe <artifacts/hadamard-7572/recipe.py>` | Standalone reconstruction and checks of the first endpoint |
| {repo-file}`T127 recipe <artifacts/t-matrices-127/recipe.py>` | Independent finite T127 input and known-order control |
| {repo-file}`Lean formalization guide <docs/FORMALIZATION.md>` | All 34 endpoint theorems, nine families and the precise formal scope |
| {repo-file}`Current portable Lean package <release/additional_source_verification/README.md>` | Complete current proof closure and 143 selected axiom reports |
| {repo-file}`Paper-to-Lean audit <docs/MANUSCRIPT_LEAN_AUDIT.md>` | Correspondence, conventions and the scope of checked assumptions |

## Further exact implementation notes

The manual gives a connected reading path. The repository also retains
specialized derivations needed by the other exact engines:

- {repo-file}`Exact relative restriction <docs/RELATIVE_SCALAR_RESTRICTION.md>`
  and {repo-file}`norm-state executor <docs/EXACT_NORM_DESCENT_EXECUTOR.md>`.
- {repo-file}`Projective endpoint flatness <docs/PROJECTIVE_FLAT_DESCENT.md>`.
- {repo-file}`Original search mechanics <docs/HADAMARD_FRAME_SEARCH.md>` and
  {repo-file}`lower-order solver guide <docs/LOW_ORDER_SEARCH.md>`.
- {repo-file}`Spence scalar families <docs/SPENCE_SCALAR_FAMILIES.md>` and
  {repo-file}`sixth-root source tripling <docs/TRIPLED_UNREAL_SCALAR_ENDPOINTS.md>`.
- {repo-file}`Research history <docs/RESEARCH_HISTORY.md>` for the locations
  of retired branches. Those branches are not part of the active reading path.
