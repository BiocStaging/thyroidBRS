# Re-deriving the 71-gene signature from TCGA-THCA

The gene list in `brs_genes` was read off a heatmap figure. The strongest
available check is to re-run the procedure that produced it and see whether
the same genes come back. `derive_signature.R` does that; this file records
what came out.

## Setup

Following Extended Experimental Procedures 14.1: equal-size sub-samples of
the BRAF-V600E and RAS mutant tumors, limma-voom at each iteration, genes at
q < 0.01, top 100 by significance, genes that rank there consistently.

| | |
| --- | --- |
| Expression | GDC STAR counts, GENCODE v36, 505 primary tumors |
| Reference groups | 234 BRAF-V600E, 52 RAS (driver mutation) |
| Gene universe | 23,339 after `filterByExpr` |
| Design | 52 vs 52 per iteration, 200 iterations |
| Published genes in the universe | 70 of 71 (`FLJ23867` does not exist in current annotation) |

The 2014 run used the TCGA RSEM gene model, roughly 20.5k mostly
protein-coding genes. Current STAR counts add ~16k lncRNAs and pseudogenes
that compete for the same 100 slots, so the comparison is run on both the
full universe and a protein-coding-restricted one.

## The published genes dominate the consistency ranking

Across 200 iterations only **262 of the 23,339 genes** ever enter a top 100,
and **68 of the 70 usable published genes are among them** — all within the
top 200 of the consistency ranking, at a median selection frequency of 92.5%.
On the coding universe the same holds with 257 genes and a median frequency
of 97%.

Taking the 71 most consistent genes:

| Universe | Published in the top 71 | p (whole universe) | p (conditioned on genes that ever appeared) |
| --- | --- | --- | --- |
| All 23,339 | 50 | 1.1e-119 | 5.6e-22 |
| Coding 16,044 | 59 | 3.3e-143 | 4.5e-36 |

The last column is the conservative test: even restricted to the ~260 genes
differentially expressed enough to reach a top 100 at all, the published list
is recovered far beyond chance.

Three published genes never reach a top 100: `ANXA2P2`, `ASAP2` and
`FLJ23867` (the last is absent from the matrix entirely).

## The strict reading is precise but narrow

Read literally, "consistently ranked at each iteration among the most
significant 100" means every single iteration:

| Design | Universe | Genes kept | Published | Precision |
| --- | --- | --- | --- | --- |
| 52 vs 52, all RAS used | all | 25 | 22 | 88% |
| 40 vs 40, both sub-sampled | all | 14 | 13 | 93% |
| 52 vs 52, all RAS used | coding | 27 | 26 | 96% |
| 40 vs 40, both sub-sampled | coding | 16 | 15 | 94% |

The genes it keeps are almost all published ones, but it keeps far fewer than
71. The three non-published genes at the strict cut are `AC012668.3`,
`AC023424.2` (lncRNAs) and `INAFM2` (annotated as `LINC00984` in 2014). None
of them existed in the gene model the 2014 analysis ran on, so measured
against what was derivable then the precision is 22/22 and 13/13.

## Most of the shortfall is the modern annotation

Re-running on a protein-coding universe (16,044 genes, approximating the
2014 RSEM gene model; `ANXA2P2` kept so the signature is not handicapped)
recovers much of the gap. The published genes are not harder to find — they
were competing for the top 100 against ~7,300 lncRNAs and pseudogenes that
did not exist as annotated genes in 2014.

| | Full universe (23,339) | Coding universe (16,044) |
| --- | --- | --- |
| Strict cut, 200 iterations | 25 genes, 22 published (88%) | 27 genes, 26 published (**96%**) |
| Top 71 by consistency | 50 of 71 (70%) | 59 of 71 (**83%**) |
| At 5 iterations | 48 genes, 37 published | 54 genes, 46 published |
| 40 vs 40 strict cut | 14 genes, 13 published | 16 genes, 15 published |
| Top 71, 40 vs 40 | — | 57 of 71 |
| Non-published at the strict cut | 3 | 1 |

The single non-published gene left at the strict cut is `INAFM2`, which was
`LINC00984` in 2014 — so against what the original analysis could have
selected, the strict cut is 26 for 26.

## How many iterations?

The paper does not say, and the count drives the result directly: surviving
"every iteration" is a much weaker filter over 5 draws than over 200.

| Iterations | Genes kept | Published | Precision |
| --- | --- | --- | --- |
| 5 | 48 | 37 | 77% |
| 10 | 43 | 34 | 79% |
| 20 | 37 | 30 | 81% |
| 50 | 32 | 26 | 81% |
| 100 | 26 | 23 | 88% |
| 200 | 25 | 22 | 88% |

Fewer iterations move the count towards 71 as expected, but not all the way
on the full universe: even at 5 iterations only 48 genes survive. On the
coding universe 5 iterations give 54 genes, 46 of them published. The two
effects — iteration count and gene universe — together account for most of
the difference; the rest is the quantification pipeline (RSEM on hg19 then,
STAR on GENCODE v36 now), which cannot be reconstructed.

## Conclusion

The signature reproduces. The published genes are not a set that a modern
re-run struggles to find; they are the top of the ranking, by a margin with
no plausible chance explanation, and on a 2014-like gene universe the strict
criterion selects essentially nothing else (26 of 27, and the odd one out was
not an annotated gene in 2014). What does not reproduce exactly is the
*count*: getting precisely 71 depends on the iteration count, the gene
universe and the quantification pipeline, none of which the paper pins down.

This is a validation of the transcribed list, not a replacement for it. The
list in `brs_genes` remains the published one.

## Reproducing this

```sh
# derive_input.rds holds list(counts, keep, coding, grp); see the script header
THCA_COUNTS=derive_input.rds BRS_UNIVERSE=coding BRS_ITER=200 BRS_CORES=10 \
    Rscript data-raw/derive_signature.R
```

Runtime is about four minutes per design on ten cores.
