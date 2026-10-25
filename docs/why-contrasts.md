# Asking why: contrasts over a round

`lupin why` (and the MCP tool `annotate_why`) answers the questions people
ask of an annotation, which are contrasts: why X rather than Y, given Z. A
cell type is asked against a rival ("why NK and not gdT?"), a gene against
another in the list it is shown in, a GO term against another term. Each
answer comes from the round's own tables, by fixed rules; nothing is
written by a language model.

## The three slots

| Slot | Meaning | Left out |
| --- | --- | --- |
| fact | what the round shows or decided: a label, a gene, a term | a label question takes the cluster's label |
| foil | the rival the asker has in mind | the strongest rival (below) |
| given | the cluster, and for a gene the list it is shown in | a gene's list is the cluster's specific genes |

An answer holds the rule that decides the question, the values on both
sides (each naming the round table it comes from), and which side they
favour: the fact, the foil, neither (the evidence does not separate them),
split (its readings disagree) or cannot tell (nothing to read it from).

## Genes

A gene is asked about within a list:

| List | Holds | Order |
| --- | --- | --- |
| `specific` | genes above the cluster's mean expression | log2 fold change over the other clusters |
| `markers:<type>` | the type's panel markers (a coarse label's: its types') | log2 fold change in the cluster, unmeasured last |
| `term:<id or name>` | the term's measured genes | log2 fold change in the cluster |

The answer gives each gene's place (rank of how many), or why it is not
shown: not measured, at or below the cluster's mean (with both values), not a
marker of the type, not in the term. It also names the types a gene marks and
the cluster's top terms that hold it. Without a foil, the foil is the gene
listed next.

## Terms

Terms rank by effect in the cluster: their genes' mean log expression
against the other annotated genes, less its median over the clusters. The
answer gives both terms' effect, z, p, q and rank, the genes that drive each
(the term's genes by fold change) and the clusters where the term's effect is
at least half this cluster's. When the order by effect and the test disagree,
a note says why: a programme many clusters share is widened in the test
(see the GO section of the README). Without a foil, the foil is the term
ranked next.

A round keeps its terms' genes in `{out}.term_genes.tsv` (term id, name,
genes) beside the term tables; a round written before this has none, and
term questions ask for it to be annotated again.

## Labels

A label's answer is a tree. Its children are the lines of evidence; each line
holds its readings and raises the gene or term question that explains it.

**The panel line.** Everything that rests on the marker panel, counted as
one line however many readings it holds:

- each method's order (enrichment, SuSiE, score, units): which of the two it
  puts first; shares within 0.01 are a tie;
- the distinguishing markers: the markers only one of the two types has,
  each set tested in the cluster against 200 random gene sets matched to its
  genes' expression strata, as GO terms are; shared markers are listed and
  left out;
- the ontology walk's placement, when the round has one: whether it put the
  cluster on either type's Cell Ontology term (the walk scores the types'
  enrichment, so it is a panel reading, not an outside one).

The line favours a side when its readings that lean all lean that way, and
is split when they lean both ways. Its child asks why the fact's top
distinct marker rather than the foil's, in the cluster's specific genes.

**The terms line.** The check from outside the panel. A term backs a side
when it is among the cluster's top 50 terms by effect, significant
(q < 0.05), and holds more of that side's distinct markers that are up in
the cluster (log2 FC ≥ 1) than of the other's. A backing term is then tested
again with every panel gene taken out of it: it counts only if what remains
still stands out (p < 0.05). A term that holds the markers and little else
repeats the panel rather than confirming it. Its child asks why the fact's
best term rather than the foil's.

**The verdict** counts the lines: it favours the fact when some line does
and none favours the foil, the foil likewise, split when lines favour both,
and neither otherwise, in which case a label over both (their group) says no
more than the evidence holds.

**The foil, when left out:** the first type the round's method puts forward
that is neither the fact nor one of its types. For a coarse label with no
rival outside it, the question it leaves open is between its own types,
which are asked instead. When the method puts forward no other type, the
other type of the smallest q.

## An example

PBMC, K0 (NK cells; labelled `natural_killer_cell`, a group of NK and gdT),
`lupin why -f pbmc.r1 -c 0`:

- the panel line favours NK: the enrichment ties (shares 0.500 against
  0.500), but NK's 8 distinct markers stand out (z 6.5) more than gdT's 2
  (z 3.5); CD160, GNLY, KLRD1 and NKG7 are shared and say nothing;
- the terms line favours NK: NK-cell mediated immunity (q 0.005) holds KLRC1
  and GZMB, and its 32 genes outside the panel still stand out (z 2.5); no
  term holds gdT's TRDC or TRGC1;
- the verdict: 2 lines for NK, none for gdT.

On the same round, K9, labelled monocyte, comes back favouring Naive CD4 T
on the panel (no monocyte marker up), and K11, labelled dendritic cell,
favouring plasma cells on both lines: labels the evidence does not back.

## Limits

- The lines are only as independent as their data: the terms line uses the
  same expression profile as the panel line, but different genes.
- A term's genes are its GO annotations: a term with an unexpected name can
  back a type through its genes (S100A9 puts "astrocyte development" among a
  monocyte cluster's terms).
- A type with few distinct markers (gdT's 2) has a weak marker test by size
  alone; the answer shows the counts.
- The test reads the clusters' pooled profiles, as the GO terms' does: no
  sample permutation.
