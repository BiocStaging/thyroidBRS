#' Fit BRAF-RAS Score (BRS) reference centroids
#'
#' Computes the two reference centroids (BRAF-mutant-like and RAS-mutant-like)
#' needed to score the BRAF-RAS Score (BRS) of Agrawal et al. (2014), from a
#' cohort where at least some samples have a known BRAF-mutant or RAS-mutant
#' label. The original study needed exome sequencing to assign these labels,
#' but that is only required once, to fit the centroids — scoring additional
#' samples afterwards ([predict.brs_fit()]) only needs their gene expression.
#'
#' @param expr A numeric matrix of gene expression, features (genes) as rows
#'   and samples as columns. Row names must be gene symbols; column names
#'   must be sample identifiers. Expression should already be on a
#'   log-like scale (e.g. `log2(TPM + 1)`); see `log2_transform` if it isn't.
#' @param labels A named character vector giving the reference class for
#'   the samples that have one. Names must match a subset of
#'   `colnames(expr)`; values must be `"Braf-like"` or `"Ras-like"` (samples
#'   with any other value, `NA`, or simply absent from `labels` are treated
#'   as unlabeled and excluded from centroid fitting, but can still be
#'   scored with [predict.brs_fit()]).
#' @param genes Character vector of gene symbols to use as the signature.
#'   Defaults to the current, alias-resolved names in [brs_genes] (68 genes
#'   after excluding the ones that never occur in `expr`).
#' @param log2_transform Logical. If `TRUE`, applies `log2(expr + 1)` to
#'   `expr` before fitting (use this if you are passing raw TPM/counts
#'   instead of already-transformed values). Default `FALSE`.
#'
#' @return An object of class `"brs_fit"`, a list with:
#'   \describe{
#'     \item{genes_used}{Signature genes actually present in `expr` and used
#'       to fit the centroids.}
#'     \item{genes_missing}{Requested signature genes not found in `expr`.}
#'     \item{genes_zero_variance}{Genes dropped because they had zero
#'       variance across the reference samples (uninformative for distance).}
#'     \item{gene_mean, gene_sd}{Per-gene mean/SD from the reference samples,
#'       used to z-score any sample later in [predict.brs_fit()].}
#'     \item{centroid_braf, centroid_ras}{Mean z-score profile of the
#'       BRAF-like and RAS-like reference samples.}
#'     \item{n_braf, n_ras}{Number of reference samples in each class.}
#'   }
#'
#' @references Agrawal N, Akbani R, Aksoy BA, et al. "Integrated Genomic
#'   Characterization of Papillary Thyroid Carcinoma." Cell.
#'   2014;159(3):676-690. \doi{10.1016/j.cell.2014.09.050}, Extended
#'   Experimental Procedures section 14.1.
#'
#' @seealso [predict.brs_fit()] to score samples with a fitted object,
#'   [validate_brs()] to check concordance against known labels.
#'
#' @examples
#' set.seed(1)
#' genes <- na.omit(brs_genes$current_symbol)[1:10]
#' samples <- paste0("S", 1:20)
#' expr <- matrix(rnorm(length(genes) * length(samples)), nrow = length(genes),
#'                dimnames = list(genes, samples))
#' labels <- setNames(rep(c("Braf-like", "Ras-like"), each = 8), samples[1:16])
#' fit <- brs_fit(expr, labels, genes = genes)
#' fit
#'
#' @importFrom stats predict
#' @export
brs_fit <- function(expr, labels, genes = NULL, log2_transform = FALSE) {
  .check_expr(expr)

  if (log2_transform) expr <- log2(expr + 1)

  if (is.null(genes)) genes <- stats::na.omit(brs_genes$current_symbol)
  genes <- unique(as.character(genes))

  genes_used <- intersect(genes, rownames(expr))
  genes_missing <- setdiff(genes, rownames(expr))
  if (length(genes_used) == 0) {
    stop("None of the requested signature genes were found in rownames(expr).", call. = FALSE)
  }

  labels <- labels[!is.na(labels)]
  labels <- labels[labels %in% c("Braf-like", "Ras-like")]
  labels <- labels[names(labels) %in% colnames(expr)]
  if (length(labels) == 0) {
    stop("No reference samples found: `labels` must name samples in `colnames(expr)` ",
         "with value \"Braf-like\" or \"Ras-like\".", call. = FALSE)
  }

  ref_expr <- expr[genes_used, names(labels), drop = FALSE]

  gene_mean <- rowMeans(ref_expr)
  gene_sd <- apply(ref_expr, 1, stats::sd)

  genes_zero_variance <- names(gene_sd)[gene_sd == 0 | is.na(gene_sd)]
  if (length(genes_zero_variance) > 0) {
    keep <- !(genes_used %in% genes_zero_variance)
    genes_used <- genes_used[keep]
    gene_mean <- gene_mean[genes_used]
    gene_sd <- gene_sd[genes_used]
    ref_expr <- ref_expr[genes_used, , drop = FALSE]
  }

  z_ref <- (ref_expr - gene_mean) / gene_sd

  is_braf <- labels[colnames(z_ref)] == "Braf-like"
  is_ras <- labels[colnames(z_ref)] == "Ras-like"

  structure(
    list(
      genes_used = genes_used,
      genes_missing = genes_missing,
      genes_zero_variance = genes_zero_variance,
      gene_mean = gene_mean,
      gene_sd = gene_sd,
      centroid_braf = rowMeans(z_ref[, is_braf, drop = FALSE]),
      centroid_ras = rowMeans(z_ref[, is_ras, drop = FALSE]),
      n_braf = sum(is_braf),
      n_ras = sum(is_ras)
    ),
    class = "brs_fit"
  )
}

#' Score samples with a fitted BRS model
#'
#' Computes the BRAF-RAS Score (BRS) and the resulting class for every
#' sample in `newdata`, using the centroids and per-gene mean/SD stored in
#' `object` (from [brs_fit()]) — `newdata` samples do not need to have been
#' part of the reference cohort used to fit those centroids, and do not need
#' any mutation/exome information.
#'
#' @param object A `"brs_fit"` object, as returned by [brs_fit()].
#' @param newdata A numeric matrix of gene expression to score (same gene
#'   symbol row-naming convention as the `expr` originally passed to
#'   [brs_fit()]). Genes in `object$genes_used` missing from `newdata` are
#'   an error, since the distance calculation would then be inconsistent
#'   with how the centroids were defined.
#' @param log2_transform Logical, as in [brs_fit()]. Must match whatever was
#'   used when fitting, since the stored per-gene mean/SD are on that scale.
#' @param ... Unused; present for S3 consistency.
#'
#' @return A data frame with one row per sample in `colnames(newdata)`:
#'   `sample`, `brs_score` (negative = more BRAF-like, positive = more
#'   RAS-like) and `brs_class` (`"Braf-like"` or `"Ras-like"`).
#'
#' @examples
#' set.seed(1)
#' genes <- na.omit(brs_genes$current_symbol)[1:10]
#' samples <- paste0("S", 1:20)
#' expr <- matrix(rnorm(length(genes) * length(samples)), nrow = length(genes),
#'                dimnames = list(genes, samples))
#' labels <- setNames(rep(c("Braf-like", "Ras-like"), each = 8), samples[1:16])
#' fit <- brs_fit(expr, labels, genes = genes)
#' predict(fit, expr[, 17:20])
#'
#' @export
predict.brs_fit <- function(object, newdata, log2_transform = FALSE, ...) {
  .check_expr(newdata)

  if (log2_transform) newdata <- log2(newdata + 1)

  missing_genes <- setdiff(object$genes_used, rownames(newdata))
  if (length(missing_genes) > 0) {
    stop(
      "`newdata` is missing ", length(missing_genes), " signature gene(s) used to fit `object`: ",
      paste(missing_genes, collapse = ", "),
      ". Re-fit with brs_fit() on a gene set newdata actually has, or subset object$genes_used.",
      call. = FALSE
    )
  }

  m <- newdata[object$genes_used, , drop = FALSE]
  z <- (m - object$gene_mean[object$genes_used]) / object$gene_sd[object$genes_used]

  sq_dist <- function(mat, centroid) colSums((mat - centroid)^2)
  score <- sq_dist(z, object$centroid_braf) - sq_dist(z, object$centroid_ras)

  data.frame(
    sample = colnames(m),
    brs_score = unname(score),
    brs_class = ifelse(score < 0, "Braf-like", "Ras-like"),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Fit and score in one step
#'
#' Convenience wrapper equivalent to `predict(brs_fit(expr, labels, ...), newdata)`,
#' for the common case of scoring the same cohort you used to fit the
#' centroids (e.g. to extend a partially-labeled cohort to full coverage).
#'
#' @inheritParams brs_fit
#' @param newdata Expression matrix to score. Defaults to `expr` itself, so
#'   `brs_score(expr, labels)` scores every sample in `expr`, including the
#'   ones without a label.
#'
#' @return See [predict.brs_fit()].
#'
#' @examples
#' set.seed(1)
#' genes <- na.omit(brs_genes$current_symbol)[1:10]
#' samples <- paste0("S", 1:20)
#' expr <- matrix(rnorm(length(genes) * length(samples)), nrow = length(genes),
#'                dimnames = list(genes, samples))
#' labels <- setNames(rep(c("Braf-like", "Ras-like"), each = 8), samples[1:16])
#' brs_score(expr, labels, genes = genes)
#'
#' @export
brs_score <- function(expr, labels, newdata = expr, genes = NULL, log2_transform = FALSE) {
  fit <- brs_fit(expr, labels, genes = genes, log2_transform = log2_transform)
  predict(fit, newdata, log2_transform = log2_transform)
}

#' Check BRS predictions against known labels
#'
#' Compares [predict.brs_fit()] output against a reference label vector —
#' typically the same `labels` used to fit the model, as a sanity check that
#' the fitted centroids reproduce the reference classes (concordance should
#' be very high, since the reference samples define the centroids), but can
#' also be any independent label vector.
#'
#' @param predictions A data frame as returned by [predict.brs_fit()] /
#'   [brs_score()].
#' @param labels A named character vector, as in [brs_fit()]. Only samples
#'   present in both `predictions$sample` and `names(labels)` are compared.
#'
#' @return A list with `n`, `n_concordant`, `pct_concordant`, and the
#'   confusion `table` (published vs. predicted).
#'
#' @examples
#' set.seed(1)
#' genes <- na.omit(brs_genes$current_symbol)[1:10]
#' samples <- paste0("S", 1:20)
#' expr <- matrix(rnorm(length(genes) * length(samples)), nrow = length(genes),
#'                dimnames = list(genes, samples))
#' labels <- setNames(rep(c("Braf-like", "Ras-like"), each = 8), samples[1:16])
#' preds <- brs_score(expr, labels, genes = genes)
#' validate_brs(preds, labels)
#'
#' @export
validate_brs <- function(predictions, labels) {
  common <- intersect(predictions$sample, names(labels))
  if (length(common) == 0) {
    stop("No samples in common between `predictions$sample` and `names(labels)`.", call. = FALSE)
  }

  published <- labels[common]
  predicted <- predictions$brs_class[match(common, predictions$sample)]

  list(
    n = length(common),
    n_concordant = sum(published == predicted),
    pct_concordant = round(100 * mean(published == predicted), 1),
    table = table(published = published, predicted = predicted)
  )
}

#' @export
print.brs_fit <- function(x, ...) {
  cat("<brs_fit>\n")
  cat("  Reference samples: ", x$n_braf, " Braf-like, ", x$n_ras, " Ras-like\n", sep = "")
  cat("  Signature genes used: ", length(x$genes_used), "\n", sep = "")
  if (length(x$genes_missing) > 0) {
    cat("  Missing from expr: ", paste(x$genes_missing, collapse = ", "), "\n", sep = "")
  }
  if (length(x$genes_zero_variance) > 0) {
    cat("  Dropped (zero variance in reference): ",
        paste(x$genes_zero_variance, collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

.check_expr <- function(expr) {
  if (!is.matrix(expr) || !is.numeric(expr)) {
    stop("`expr`/`newdata` must be a numeric matrix (genes x samples).", call. = FALSE)
  }
  if (is.null(rownames(expr)) || is.null(colnames(expr))) {
    stop("`expr`/`newdata` must have rownames (gene symbols) and colnames (sample IDs).", call. = FALSE)
  }
}
