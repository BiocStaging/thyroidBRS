# thyroidBRS 0.99.0

## Breaking changes

* `predict()` now returns the score defined in Agrawal et al., Extended
  Experimental Procedures 14.1 — the difference of **normalized Euclidean
  distances** to the two centroids. Previous versions returned the difference
  of *squared* distances, which has the same sign (so classes are unchanged)
  but a different magnitude and a different rank order.
* The output gains a `brs_scaled` column: the score rescaled to `[-1, 1]`
  across the scored samples, which is the axis the published TCGA-THCA values
  use.
* `brs_fit()` now records `log2_transform` in the fitted object and
  `predict()` reuses it. Previously the two could silently disagree, which
  produced meaningless scores with no error or warning.
* `brs_fit()` accepts driver-mutation labels (`"BRAF_V600E"`, `"RAS"`) as well
  as the previous `"Braf-like"`/`"Ras-like"` values, and the documentation now
  steers towards mutation status. The reference groups in the paper are
  defined by genotype; the published `BRAF_RAF_class` column is the
  classifier's own output and fitting on it is circular.

## Fixes

* `brs_fit()` errors instead of silently producing `NaN` centroids when a
  reference group has fewer than two samples.
* `brs_fit()` warns when signature genes are absent from `expr`, instead of
  dropping them silently.
* Signature genes are now resolved against whichever symbol spelling the
  expression matrix uses. `ARNTL` was renamed `BMAL1` after the paper, so
  GENCODE v36 matrices and current ones spell it differently; both now give
  the full 70 genes.
* `predict()` warns when samples cannot be scored because of missing values,
  and when duplicated row names force it to pick one row per gene.
* A score of exactly zero is documented as `"Ras-like"` rather than left to
  the reader.
* `validate_brs()` gains a `fit` argument, reports `n_resubstituted`, and
  warns when every compared sample also defined the centroids.

## Bioconductor

* Expression can now be a `SummarizedExperiment` or an `ExpressionSet` as
  well as a matrix, with `assay=` to pick the assay and `labels` optionally
  naming a `colData()`/`pData()` column. See `?containers`.
* `predict()` and `brs_score()` gain `standardize`, choosing how `newdata` is
  put on the centroids' scale: `"reference"` (the default, unchanged
  behaviour), `"cohort"` or `"rank"`. The default keeps a sample's score
  independent of the company it is scored with; the other two trade that away
  to make the class threshold usable across platforms. See
  `?predict.brs_fit`.
* `predict()` warns when every scored sample falls on the same side of zero,
  which is what a cross-platform threshold mismatch looks like.
* Added `biocViews`, `inst/CITATION`, and a `BiocStyle` vignette; version
  numbering moved to the `0.99.z` series used for Bioconductor submission.

## Other

* `brs_genes` gains `block` and `up_in`: the two row-clusters of Figure S7A
  and the group each is higher in. Refitting on TCGA-THCA recovers this split
  for 70 of 70 genes.
* `data-raw/` documents the provenance of the gene list
  (`brs_genes.md`), re-checks the symbols (`check_symbols.R`), reproduces
  the TCGA validation (`validate_tcga.R`), re-derives the signature from the
  data (`derive_signature.R`, results in `signature_rederivation.md`) and
  records what does and does not transfer to microarrays
  (`microarray_notes.md`).
* `inst/extdata/thca_reference.csv` carries the TCGA-THCA driver-mutation
  status and published BRS, so the validation does not depend on the paper's
  supplementary spreadsheet.

# thyroidBRS 0.1.0

* Initial release.
