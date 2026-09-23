# Does the BRS transfer to microarray data?

Short answer: the **axis** transfers, the **threshold** does not.

Tested on GSE33630 (Affymetrix HG-U133 Plus 2.0, 49 papillary thyroid
carcinomas, RMA, log2), scored against centroids fitted on TCGA-THCA RNA-seq.
GSE33630 has no BRAF/RAS mutation status in its GEO metadata, so there is no
ground truth on the array side; the checks below are therefore internal.

## The signature is intact on the platform

67 of the 70 signature genes the centroids were fitted on are present on this
matrix (`ITGA3` and `RASGEF1B` were lost in the upstream filtering, not absent
from the array; a third drops out of the fit). Across the 49 tumors:

| | |
| --- | --- |
| Variance of the signature explained by PC1 | 44.9% |
| Genes whose PC1 loading sign matches their Figure S7A block | 66 / 67 |
| Mean correlation within block 1 | +0.284 |
| Mean correlation within block 2 | +0.451 |
| Mean correlation between blocks | -0.336 |

`data-raw/portability.R` reproduces every figure on this page.

The two-block structure that defines the BRS is fully present. The biology is
not the problem.

## The ranking transfers; the zero point does not

Scoring the same 49 tumors three ways, against the same TCGA centroids, and
comparing each to the array's own signature axis (PC1, oriented like the BRS):

| Standardization of the array data | Spearman vs. PC1 | Classes called |
| --- | --- | --- |
| Reference (TCGA per-gene mean/SD, what the package does) | +0.976 | 49 BRAF / 0 RAS |
| Cohort-internal (re-standardize against GSE33630) | +0.979 | 41 BRAF / 8 RAS |
| Within-sample ranks, then standardized | +0.846 | 47 BRAF / 2 RAS |

Every option ranks the tumors nearly identically. They disagree completely on
where zero falls, and that is the whole class call.

## Why cohort-internal standardization is not the fix

It is the obvious repair and it does solve the cross-platform threshold, but
it makes each sample's score depend on the company it keeps. Scoring TCGA
subsets of different composition against the same centroids:

| Cohort scored | n | Reference z-score | Cohort-internal |
| --- | --- | --- | --- |
| All published (272 Braf / 119 Ras) | 391 | 99.0% | 97.4% |
| Only Braf-like | 272 | 99.6% | 85.7% |
| Only Ras-like | 119 | 97.5% | **8.4%** |
| Skewed 90/10 | 200 | 99.5% | 97.0% |
| Balanced 50/50 | 238 | 98.3% | 86.1% |

Frozen reference standardization is insensitive to the class mix of the set
being scored, which is exactly the property a classifier should have. That is
why it stays the default.

## Recommendation

- Same platform and pipeline as the fit: use the package as is.
- Different platform, mutation labels available: fit new centroids on that
  platform. This is the only option that gives a trustworthy class call.
- Different platform, no labels: use `brs_score` as a ranking and do not
  report classes. `predict()` warns when every sample lands on one side,
  which is the usual symptom of this situation.

`standardize = c("reference", "cohort", "rank")` makes the trade-off
selectable rather than implicit; `"reference"` stays the default for the
reason the table above gives.
