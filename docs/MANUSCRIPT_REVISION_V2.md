# Manuscript revision 2

The current paper is **Butson towers and high-rank Hadamard matrices**.
Its canonical source is
[`manuscript/hadamard_7572_announcement.tex`](../manuscript/hadamard_7572_announcement.tex),
with the current PDF, API metadata and source archive under
[`release/butson_towers_v2`](../release/butson_towers_v2/README.md).

The abstract foregrounds eighteen new real Hadamard ranks and **six new
output families**. Earlier summaries counted the five Spence-source
families while describing the cubic family separately. The correction
adds no unproved mathematical claim: `unreal_cubic_family` was already
included in the completed 85-target Lean audit at commit `686acf9`.

The paper defines Butson towers as scalar-reduction chains retaining the
action, conjugation and trace metric. Its Spence families are the five
new T631-derived real families whose source recurrence is Spence's
classical construction. The sixth real family uses classical unreal-source
tripling. The source recurrences and output-family novelty are explicitly
distinguished. Known individual ranks 885924 and 18629644 remain marked known.

The revision expands the main theorem to all eighteen ranks and all six
families, records the completed Lean existence proofs, and cites the 2026
starting paper together with the earlier scalar-morphism constructions.
Rank and order coincide for these nonsingular square matrices.

The original announcement and construction bundles retain their bytes and
hashes. The separate Lean reproducibility packet is also unchanged. This
is a local manuscript revision; no remote Hexagon draft or submission is
created by its build and preparation scripts.
