# Completed matrices in a scalar-restriction tree

`large_butson_factory.py` starts with four **existing verified matrices**:
\(BH(5329,24)\), \(BH(58081,120)\), \(BH(1000,8)\), and
\(BH(1000,16)\). It does not search new source constructions or invoke SAT.
Only certified repeated-prime scalar restrictions create children:
\(BH(n,q)\to BH(pn,q/p)\) for \(p^2\mid q\).

The queue processes matrices breadth first. It tries both already implemented
generator conventions, retains the actual entry alphabet, and verifies each
completed matrix through its original source Gram and exact scalar identities.
Each distinct node receives a compressed literal-entry bundle and a JSON
proof record. An edge reaching an identical source-plus-composed-block-table
definition is recorded as a duplicate and does not expand another copy.
This is exact, limited duplicate detection; it does not classify arbitrary
Hadamard equivalence.

The four sources give the order/phase chains

\[
(5329,24)\to(10658,12)\to(21316,6),
\]
\[
(58081,120)\to(116162,60)\to(232324,30),
\]
\[
(1000,8)\to(2000,4)\to(4000,2),
\]
\[
(1000,16)\to(2000,8)\to(4000,4)\to(8000,2).
\]

Distinct generator conventions can give more than one matrix at a listed
order/phase. The saved run report records actual node counts and duplicate
edges. The sources use established constructions; these chains do not claim
new Butson orders or new real Hadamard orders.

The character sources are saved by their complete \(q\times q\) coefficient
arrays and group-development formula, rather than allocating \(q^2\times q^2\)
dense matrices. Their full Grams follow from checked character group-ring
identities. Descendants retain explicit scalar blocks, exact proofs, and
literal entry access. A completed node is therefore a matrix construction,
not an unsolved constraint problem.

Stopping at phase 6 or 30 means only that the implemented repeated-prime
steps have ended. Other cyclotomic restrictions, and intermediate matrices
with a retained metric that are not Butson matrices, remain open operations.
They are not excluded by this finite tree. The character sources' raw
alphabets also have a separately proved obstruction to a uniform final sign
map; it does not invalidate any of their completed Butson descendants.

The CLI accepts `--database`, `--seconds`, `--report`, `--max-db-mib` and
an optional `--nodes` limit. Its SQLite cursor resumes unprocessed matrix
recipes. Files live in the sibling `<database>.matrices` directory and count
toward the storage budget. Tests check a completed order-4000 chain,
duplicate pruning, breadth-first order, locking and resume without repeated
work. The separate descendant-adapter tests replay the literal bundles.
