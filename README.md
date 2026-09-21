# thyroidBRS

An R package that computes the **BRAF-RAS Score (BRS)** for papillary
thyroid carcinoma (PTC) — a 71-gene expression signature that classifies a
tumor as more BRAF^V600E^-mutant-like or more RAS-mutant-like — as defined
in:

> Agrawal N, Akbani R, Aksoy BA, et al. Integrated Genomic Characterization
> of Papillary Thyroid Carcinoma. *Cell*. 2014;159(3):676-690.
> [doi:10.1016/j.cell.2014.09.050](https://doi.org/10.1016/j.cell.2014.09.050)

## Why this package exists

The original TCGA-THCA study only published a BRS for 391 of its 496
samples, because computing the score required **both** exome sequencing
(to define which tumors are confidently BRAF-mutant vs. RAS-mutant, i.e.
the two reference groups) **and** RNA sequencing (to compute the score
itself). But once the two reference centroids are fixed, scoring an
*additional* tumor only needs its gene expression — no mutation calling
required. At the time of writing, no public package implements this
classifier, and the exact gene list is only published as a heatmap figure
(Figure S7A), not as a machine-readable table.

`thyroidBRS`:

- ships the 71-gene signature, transcribed from Figure S7A and reconciled
  against current HGNC gene symbols (`brs_genes`);
- lets you fit the two reference centroids (`brs_fit()`) from any cohort
  with at least some known BRAF-mutant/RAS-mutant labels;
- scores *any* additional sample's expression against those centroids
  (`predict()`), including samples with no mutation data at all.

Re-running the original TCGA-THCA cohort with this package (using the 391
already-labeled samples as the reference, then scoring the full 505-patient
cohort available today) reaches **98.7% concordance** (386/391) with the
published labels, and extends classification to previously-unclassified
patients with no exome data.

## Installation

```r
# install.packages("remotes")
remotes::install_github("camilasauria/thyroidBRS")
```

## Quick start

```r
library(thyroidBRS)

# expr: log2(TPM + 1) expression matrix, genes (symbols) x samples
# labels: named character vector, values "Braf-like"/"Ras-like", for the
#         subset of samples with a known driver mutation
fit <- brs_fit(expr, labels)
predict(fit, expr)  # scores every sample, including unlabeled ones
```

See `vignette("thyroidBRS")` for a full worked example, including how to
pull TCGA-THCA data with `TCGAbiolinks` and reproduce the validation above.

## Gene list caveat

The 71 genes were transcribed manually from the published heatmap (no
machine-readable list was released with the paper). The gene count matches
exactly (71), and 3 outdated symbols were resolved to their current HGNC
name (`PVRL4` → `NECTIN4`, `FAM176A` → `EVA1A`, `TM7SF4` → `DCSTAMP`). One
gene, `FLJ23867`, could not be resolved to any symbol in current genome
annotations and is excluded (`brs_genes$current_symbol` is `NA` for it) —
the classifier uses the remaining 70 genes.

## License

MIT © Camila Barrios
