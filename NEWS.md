# thyroidBRS 0.1.0

* Initial release.
* `brs_genes`: the 71-gene BRAF-RAS Score (BRS) signature from Agrawal et al.
  (2014), Figure S7A, with outdated gene symbols resolved to current HGNC
  names.
* `brs_fit()`: fit BRAF-like/RAS-like reference centroids from a cohort with
  known driver-mutation labels.
* `predict.brs_fit()`: score any sample's expression against fitted
  centroids, including samples without mutation data.
* `brs_score()`: one-step `brs_fit()` + `predict()` convenience wrapper.
* `validate_brs()`: concordance check against known labels.
