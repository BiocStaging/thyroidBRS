# Provenance of the 71-gene BRS signature

## Source

Agrawal N, Akbani R, Aksoy BA, et al. *Integrated Genomic Characterization of
Papillary Thyroid Carcinoma.* Cell. 2014;159(3):676-690.
doi:10.1016/j.cell.2014.09.050 — **Figure S7A**, "mRNA expression of the 71
genes used to derive the BRAF-V600E-RAS score (BRS) across 391 PTC samples".

The gene list was never released in machine-readable form: it exists only as
the row labels of that heatmap. The symbols in `brs_genes` were read off the
figure and then re-checked against it.

The figure is panel A of the supplementary figure file distributed with the
PubMed Central deposit of the paper (`NIHMS633017-supplement-7.pdf`). It is
copyrighted and is therefore **not** redistributed with this package; retrieve
it from the publisher or from PMC if you want to re-check the transcription.

## Checks that were run

1. **Row-by-row comparison against Figure S7A.** All 71 labels match, in
   order, in both row-clusters.
2. **Block structure.** The figure has two row-clusters, each printed in
   alphabetical order: 13 genes (`ANKRD46` … `SORBS2`) and 58 genes (`ABTB2`
   … `TMEM43`). Both runs are strictly sorted with no gene out of place, which
   is what makes an omission or a misread detectable. `13 + 58 = 71` matches
   the count stated in the figure legend. The two blocks move in opposite
   directions between the BRAF-mutant and RAS-mutant groups; the `block`
   column of `brs_genes` records the assignment.
3. **Symbol resolution** against `org.Hs.eg.db` — see `check_symbols.R` in
   this directory. Four symbols are stale (`ARNTL`→`BMAL1`,
   `FAM176A`→`EVA1A`, `PVRL4`→`NECTIN4`, `TM7SF4`→`DCSTAMP`); `FLJ23867`
   resolves to nothing and is excluded from the signature.

4. **Re-derivation from the data.** The procedure of Extended Experimental
   Procedures 14.1 was re-run on TCGA-THCA — see `derive_signature.R` and
   `signature_rederivation.md`. Across 200 iterations only ~260 of 16,000-23,000
   genes ever reach a top 100, and 68 of the 70 usable published genes are
   among them at a median selection frequency of 92-97%. On a 2014-like
   protein-coding universe the strict criterion keeps 27 genes of which 26 are
   published, and the one exception was not an annotated gene in 2014.

## Still open

The re-derivation does not reproduce the count of exactly 71. That number
depends on the iteration count, the gene universe and the quantification
pipeline, none of which the paper specifies, so it is not recoverable — see
`signature_rederivation.md` for how far each of those three moves it. The
list shipped in `brs_genes` is the published one, not a re-derived one.
