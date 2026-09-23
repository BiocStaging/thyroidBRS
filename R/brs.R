# Provenance: substantial parts of this file were written by Claude Opus 5
# (Anthropic) working from the published description of the method, under the
# direction of and reviewed by the package authors, who are responsible for
# its correctness and maintenance. See the Provenance section of README.md.
# Assisted-by: Claude Opus 5 (Anthropic)

#' Fit BRAF-RAS Score (BRS) reference centroids
#'
#' Computes the two reference centroids (BRAF-mutant and RAS-mutant) needed to
#' score the BRAF-RAS Score (BRS) of Agrawal et al. (2014), from a cohort where
#' at least some samples have a known driver mutation. The original study
#' needed exome sequencing to assign those labels, but that is only required
#' once, to fit the centroids. Scoring additional samples afterwards
#' ([predict.brs_fit()]) only needs their gene expression.
#'
#' @section Which samples define the reference groups:
#' In Agrawal et al. the two reference sets are defined by **driver mutation**:
#' set B is the BRAF-V600E-mutant tumors and set R is the RAS-mutant tumors.
#' They are *not* the BRAF-like/RAS-like classes that the BRS itself produces.
#'
#' This distinction matters when working with the published TCGA-THCA tables.
#' The column `BRAF_RAF_class` (also exposed by
#' `TCGAbiolinks::TCGAquery_subtype()`) is exactly `sign(BRAF_RAF_score)`:
#' it is the *output* of the classifier.
#' Fitting centroids on it is circular: it pulls 43 tumors with other drivers
#' and 1 RAS-mutant tumor into the "BRAF" centroid, and 68 other-driver tumors
#' into the "RAS" centroid. Use the mutation column (`BRAFV600E_RAS`,
#' restricted to `BRAF_V600E` and `RAS`) to fit, and keep `BRAF_RAF_class` /
#' `BRAF_RAF_score` for validation. See `vignette("thyroidBRS")`.
#'
#' @param expr Gene expression, features (genes) as rows and samples as
#'   columns: a numeric matrix, a
#'   [SummarizedExperiment::SummarizedExperiment] or a
#'   `Biobase::ExpressionSet` (see [containers]). Row names must be gene
#'   symbols; column names must be sample identifiers. Expression should
#'   already be on a log-like scale (e.g. `log2(TPM + 1)`); see
#'   `log2_transform` if it isn't.
#' @param labels A named character vector giving the reference group for the
#'   samples that have one, or a single string naming a column of the
#'   object's `colData()`/`pData()`. Names must match a subset of
#'   `colnames(expr)`.
#'   Values naming the BRAF group (`"BRAF_V600E"`, `"BRAF"`, `"BRAF-like"`,
#'   `"Braf-like"`) and the RAS group (`"RAS"`, `"RAS-like"`, `"Ras-like"`)
#'   are recognised, case-insensitively. Samples with any other value, `NA`,
#'   or simply absent from `labels` are treated as unlabeled and excluded from
#'   centroid fitting, but can still be scored with [predict.brs_fit()].
#' @param genes Character vector of gene symbols to use as the signature.
#'   Defaults to the alias-resolved names in [brs_genes] (70 genes; `FLJ23867`
#'   has no current symbol and is excluded).
#' @param log2_transform Logical. If `TRUE`, applies `log2(expr + 1)` to `expr`
#'   before fitting (use this if you are passing raw TPM/counts instead of
#'   already-transformed values). The choice is stored in the returned object
#'   and reused by [predict.brs_fit()], so the two can no longer disagree.
#'   Default `FALSE`.
#' @param assay Name or index of the assay to use when `expr` is a
#'   `SummarizedExperiment`. Ignored otherwise. Defaults to the first assay,
#'   with a warning when there is more than one. For GDC data the first
#'   assay is raw counts, so use `assay = "tpm_unstranded"` with
#'   `log2_transform = TRUE`, or a normalized assay of your own.
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
#'       BRAF-mutant and RAS-mutant reference samples.}
#'     \item{n_braf, n_ras}{Number of reference samples in each group.}
#'     \item{reference_samples}{Names of the samples used to fit, so that
#'       [validate_brs()] can tell resubstitution from real validation.}
#'     \item{log2_transform}{The scale the centroids were fitted on.}
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
#' expr <- matrix(rnorm(length(genes) * length(samples)),
#'     nrow = length(genes),
#'     dimnames = list(genes, samples)
#' )
#' labels <- setNames(rep(c("BRAF_V600E", "RAS"), each = 8), samples[1:16])
#' fit <- brs_fit(expr, labels, genes = genes)
#' fit
#'
#' @importFrom stats predict
#' @export
brs_fit <- function(expr, labels, genes = NULL, log2_transform = FALSE,
                    assay = NULL) {
    .check_flag(log2_transform, "log2_transform")
    labels <- .labels_from_object(expr, labels)
    expr <- .as_expr_matrix(expr, assay, "expr")
    .check_expr(expr)

    if (log2_transform) expr <- .log2p1(expr)

    sel <- .select_genes(expr, genes)
    expr <- sel$expr
    labels <- .normalise_labels(labels, colnames(expr))

    out <- .fit_centroids(expr[sel$used, names(labels), drop = FALSE], labels)
    out$genes_missing <- sel$missing
    out$reference_samples <- names(labels)
    out$log2_transform <- log2_transform

    structure(
        out[c(
            "genes_used", "genes_missing", "genes_zero_variance",
            "gene_mean", "gene_sd", "centroid_braf", "centroid_ras",
            "n_braf", "n_ras", "reference_samples",
            "log2_transform"
        )],
        class = "brs_fit"
    )
}

# Resolve the signature against the matrix and report what is missing.
.select_genes <- function(expr, genes) {
    if (is.null(genes)) {
        genes <- .resolve_signature(rownames(expr))
    } else {
        genes <- unique(stats::na.omit(as.character(genes)))
    }
    expr <- .drop_duplicate_rows(expr, genes)

    used <- intersect(genes, rownames(expr))
    missing <- setdiff(genes, rownames(expr))
    if (length(used) == 0) {
        stop("None of the requested signature genes were found in ",
            "rownames(expr). Are the row names gene symbols?",
            call. = FALSE
        )
    }
    if (length(missing) > 0) {
        warning("Fitting on ", length(used), " of ", length(genes),
            " signature genes; ", length(missing),
            " not found in rownames(expr): ",
            .truncate(missing),
            call. = FALSE
        )
    }
    list(expr = expr, used = used, missing = missing)
}

# Per-gene standardization over the reference samples, then the two centroids.
.fit_centroids <- function(ref_expr, labels) {
    gene_mean <- rowMeans(ref_expr)
    gene_sd <- apply(ref_expr, 1, stats::sd)

    incomplete <- names(gene_sd)[is.na(gene_sd)]
    if (length(incomplete) > 0) {
        warning(length(incomplete), " signature gene(s) have missing values ",
                "in the reference samples and were dropped: ",
                .truncate(incomplete),
                ". A single NA removes the gene for every sample.",
                call. = FALSE)
    }
    flat <- names(gene_sd)[!is.na(gene_sd) & gene_sd == 0]
    drop <- union(flat, incomplete)
    if (length(drop) > 0) {
        keep <- setdiff(rownames(ref_expr), drop)
        if (length(keep) == 0) {
            stop("Every signature gene has zero variance (or is all NA) ",
                "across the reference samples.",
                call. = FALSE
            )
        }
        gene_mean <- gene_mean[keep]
        gene_sd <- gene_sd[keep]
        ref_expr <- ref_expr[keep, , drop = FALSE]
    }

    z_ref <- (ref_expr - gene_mean) / gene_sd
    is_braf <- labels[colnames(z_ref)] == "BRAF"
    is_ras <- labels[colnames(z_ref)] == "RAS"

    list(
        genes_used = rownames(ref_expr),
        genes_zero_variance = flat,
        gene_mean = gene_mean,
        gene_sd = gene_sd,
        centroid_braf = rowMeans(z_ref[, is_braf, drop = FALSE]),
        centroid_ras = rowMeans(z_ref[, is_ras, drop = FALSE]),
        n_braf = sum(is_braf),
        n_ras = sum(is_ras)
    )
}

#' Score samples with a fitted BRS model
#'
#' Computes the BRAF-RAS Score (BRS) and the resulting class for every sample
#' in `newdata`, using the centroids and per-gene mean/SD stored in `object`
#' (from [brs_fit()]). `newdata` samples do not need to have been part of the
#' reference cohort used to fit those centroids, and do not need any
#' mutation/exome information.
#'
#' @section Scoring data from another platform:
#' The centroids carry the per-gene mean and SD of the cohort they were
#' fitted on. That is what makes a sample's score its own: it does not move
#' when you score it alongside a different set of samples. It is also what
#' breaks when `newdata` comes from another platform, because those means and
#' SDs no longer describe it.
#'
#' Measured on GSE33630 (Affymetrix, 49 papillary carcinomas) against
#' centroids fitted on TCGA RNA-seq: all three settings rank the tumors
#' almost identically, at Spearman 0.97, 0.98 and 0.92 against the array's
#' own signature axis, and they disagree completely on where zero falls.
#' `"reference"` calls 49 of 49 BRAF-like; `"cohort"` calls 42 and 7.
#'
#' `"cohort"` and `"rank"` fix the threshold across platforms by paying for
#' it elsewhere: each sample's score then depends on which samples it is
#' scored with. Scoring TCGA subsets of different composition against the
#' same centroids, `"reference"` gets 97.5% to 99.6% right whatever the class
#' mix, while `"cohort"` drops to 8.4% on a cohort that happens to be all
#' RAS-like.
#'
#' So: score on the platform you fitted on and keep the default. Crossing
#' platforms, prefer fitting new centroids there if you have driver-mutation
#' labels. Failing that, use these settings for the *ranking* and treat the
#' class call as unreliable. `predict()` warns when every sample lands on the
#' same side of zero, which is the usual symptom.
#'
#' @section How the score is computed:
#' Following Agrawal et al. (2014), Extended Experimental Procedures 14.1, the
#' score of a tumor `t` is the difference between its distance to the two
#' centroids,
#'
#' \deqn{BRS(t) = \|v(t) - c(B)\|_2 - \|v(t) - c(R)\|_2}
#'
#' where the distance is a normalized Euclidean distance. Here `v(t)` is the
#' sample's expression vector z-scored gene-by-gene against the reference
#' cohort, and the norm is divided by the square root of the number of
#' signature genes, so that scores stay comparable when the gene set shrinks
#' (for example because some signature genes are absent from `newdata`).
#'
#' Negative scores are BRAF-V600E-like and positive scores RAS-like. A score
#' of exactly `0` is not defined by the paper; it is assigned to `"Ras-like"`
#' here, consistent with treating the non-negative half as the RAS side.
#'
#' The BRS values published for TCGA-THCA are additionally rescaled to
#' \eqn{[-1, +1]}, by dividing the negative scores by the most negative value
#' and the positive scores by the largest one. `brs_scaled` reproduces that
#' rescaling. Note that it is a property of the *set of samples being scored*,
#' not of a sample on its own: scoring a different cohort, or a single sample,
#' gives a different (or undefined) rescaling. Use `brs_score` when you need a
#' value that is stable per sample, and `brs_scaled` when you need a value on
#' the same axis as the published figures.
#'
#' @param object A `"brs_fit"` object, as returned by [brs_fit()].
#' @param newdata Gene expression to score, in any of the forms `expr`
#'   accepts (see [containers]) and with the same gene symbol row naming as
#'   the `expr` originally passed to [brs_fit()]. Genes in
#'   `object$genes_used` missing from `newdata` are an error, since the
#'   distance calculation would then be inconsistent with how the centroids
#'   were defined.
#' @param log2_transform Logical, or `NULL` (the default) to reuse whatever
#'   was used when fitting `object`. The stored per-gene mean/SD are on that
#'   scale, so overriding it is almost always a mistake and warns.
#' @param standardize How to put `newdata` on the scale the centroids live
#'   on. `"reference"` (the default) uses the per-gene mean and SD frozen from
#'   the fit. `"cohort"` re-derives them from `newdata` itself, and `"rank"`
#'   replaces each sample's expression by within-sample ranks first. The last
#'   two need at least 3 samples. See the section below before changing it.
#' @param assay Assay to use when `newdata` is a `SummarizedExperiment`; see
#'   [brs_fit()].
#' @param ... Unused; present for S3 consistency.
#'
#' @return A data frame with one row per sample in `colnames(newdata)`:
#'   `sample`, `brs_score` (normalized distance difference; negative = more
#'   BRAF-like, positive = more RAS-like), `brs_scaled` (the same score
#'   rescaled to \eqn{[-1, +1]} across the scored samples, `NA` when it cannot
#'   be determined) and `brs_class` (`"Braf-like"` or `"Ras-like"`).
#'
#' @examples
#' set.seed(1)
#' genes <- na.omit(brs_genes$current_symbol)[1:10]
#' samples <- paste0("S", 1:20)
#' expr <- matrix(rnorm(length(genes) * length(samples)),
#'     nrow = length(genes),
#'     dimnames = list(genes, samples)
#' )
#' labels <- setNames(rep(c("BRAF_V600E", "RAS"), each = 8), samples[1:16])
#' fit <- brs_fit(expr, labels, genes = genes)
#' predict(fit, expr[, 17:20])
#'
#' @export
predict.brs_fit <- function(object, newdata, log2_transform = NULL,
                            standardize = c("reference", "cohort", "rank"),
                            assay = NULL, ...) {
    standardize <- match.arg(standardize)
    newdata <- .as_expr_matrix(newdata, assay, "newdata")
    .check_expr(newdata)

    if (.resolve_scale(object, log2_transform)) newdata <- .log2p1(newdata)

    newdata <- .drop_duplicate_rows(newdata, object$genes_used)
    newdata <- .match_alias_rows(newdata, object$genes_used)

    missing_genes <- setdiff(object$genes_used, rownames(newdata))
    if (length(missing_genes) > 0) {
        stop("`newdata` is missing ", length(missing_genes),
            " signature gene(s) used to fit `object`: ",
            .truncate(missing_genes),
            ". Re-fit with brs_fit() on a gene set newdata actually has, ",
            "or subset object$genes_used.",
            call. = FALSE
        )
    }

    m <- newdata[object$genes_used, , drop = FALSE]
    st <- .standardize(m, object, standardize)
    z <- st$z
    g <- st$genes
    if (standardize == "reference") .warn_if_offscale(z)

    # Normalized Euclidean distance: the L2 norm divided by sqrt(p), so the
    # score does not depend on how many signature genes survived.
    p <- length(g)
    dist_to <- function(centroid) sqrt(colSums((z - centroid[g])^2) / p)
    score <- dist_to(object$centroid_braf) - dist_to(object$centroid_ras)

    .warn_about_scores(score)

    data.frame(
        sample = colnames(m),
        brs_score = unname(score),
        brs_scaled = unname(.rescale_brs(score)),
        brs_class = ifelse(is.na(score), NA_character_,
            ifelse(score < 0, "Braf-like", "Ras-like")
        ),
        row.names = NULL,
        stringsAsFactors = FALSE
    )
}

#' Fit and score in one step
#'
#' Convenience wrapper equivalent to
#' `predict(brs_fit(expr, labels, ...), newdata)`, for the common case of
#' scoring the same cohort you used to fit the centroids (e.g. to extend a
#' partially-labeled cohort to full coverage).
#'
#' @inheritParams brs_fit
#' @param newdata Expression matrix to score. Defaults to `expr` itself, so
#'   `brs_score(expr, labels)` scores every sample in `expr`, including the
#'   ones without a label.
#' @param standardize Passed to [predict.brs_fit()].
#'
#' @return See [predict.brs_fit()].
#'
#' @examples
#' set.seed(1)
#' genes <- na.omit(brs_genes$current_symbol)[1:10]
#' samples <- paste0("S", 1:20)
#' expr <- matrix(rnorm(length(genes) * length(samples)),
#'     nrow = length(genes),
#'     dimnames = list(genes, samples)
#' )
#' labels <- setNames(rep(c("BRAF_V600E", "RAS"), each = 8), samples[1:16])
#' brs_score(expr, labels, genes = genes)
#'
#' @export
brs_score <- function(expr, labels, newdata = expr, genes = NULL,
                      log2_transform = FALSE,
                      standardize = c("reference", "cohort", "rank"),
                      assay = NULL) {
    fit <- brs_fit(expr, labels,
        genes = genes,
        log2_transform = log2_transform, assay = assay
    )
    predict(fit, newdata, standardize = match.arg(standardize), assay = assay)
}

#' Check BRS predictions against known labels
#'
#' Compares [predict.brs_fit()] output against a reference label vector.
#'
#' @section Resubstitution is not validation:
#' Passing the same `labels` that were used to fit the centroids measures how
#' well a nearest-centroid rule reproduces the labels that *defined* those
#' centroids. That statistic is strongly optimistic: on 70 genes and 391
#' samples split 272/119, pure noise with no signal at all still resubstitutes
#' at about 68%. Pass `fit` and this function will tell you how much of the
#' comparison is resubstitution, and warn when all of it is.
#'
#' For a real check, compare against labels that did not enter the fit.
#' For TCGA-THCA, fit on the mutation groups and validate against the
#' published `BRAF_RAF_class` and `BRAF_RAF_score`.
#'
#' @param predictions A data frame as returned by [predict.brs_fit()] /
#'   [brs_score()].
#' @param labels A named character vector, as in [brs_fit()]. Only samples
#'   present in both `predictions$sample` and `names(labels)` are compared.
#' @param fit Optionally the `"brs_fit"` used to produce `predictions`, so the
#'   overlap with its reference samples can be reported.
#'
#' @return A list with `n`, `n_concordant`, `pct_concordant`, `n_resubstituted`
#'   (samples that also defined the centroids, `NA` if `fit` was not given) and
#'   the confusion `table` (published vs. predicted).
#'
#' @examples
#' set.seed(1)
#' genes <- na.omit(brs_genes$current_symbol)[1:10]
#' samples <- paste0("S", 1:20)
#' expr <- matrix(rnorm(length(genes) * length(samples)),
#'     nrow = length(genes),
#'     dimnames = list(genes, samples)
#' )
#' labels <- setNames(rep(c("BRAF_V600E", "RAS"), each = 8), samples[1:16])
#' fit <- brs_fit(expr, labels, genes = genes)
#' preds <- predict(fit, expr)
#' validate_brs(preds, labels, fit = fit)
#'
#' @export
validate_brs <- function(predictions, labels, fit = NULL) {
    groups <- .normalise_labels(labels, predictions$sample,
        what = "predictions$sample", min_per_group = 0L
    )
    labels <- c(BRAF = "Braf-like", RAS = "Ras-like")[groups]
    names(labels) <- names(groups)

    common <- intersect(predictions$sample, names(labels))
    predicted <- predictions$brs_class[match(common, predictions$sample)]
    keep <- !is.na(predicted)
    if (!all(keep)) {
        warning(sum(!keep), " sample(s) had no predicted class and were ",
            "dropped from the comparison.",
            call. = FALSE
        )
        common <- common[keep]
        predicted <- predicted[keep]
    }
    if (length(common) == 0) {
        stop("No samples left to compare between `predictions` and `labels`.",
            call. = FALSE
        )
    }
    published <- unname(labels[common])

    n_resub <- NA_integer_
    if (!is.null(fit)) {
        if (!inherits(fit, "brs_fit")) {
            stop("`fit` must be a \"brs_fit\" object.", call. = FALSE)
        }
        n_resub <- length(intersect(common, fit$reference_samples))
        if (n_resub == length(common)) {
            warning("Every compared sample also defined the centroids: this ",
                "is a resubstitution estimate, not validation. See ",
                "?validate_brs.",
                call. = FALSE
            )
        }
    }

    list(
        n = length(common),
        n_concordant = sum(published == predicted),
        pct_concordant = round(100 * mean(published == predicted), 1),
        n_resubstituted = n_resub,
        table = table(published = published, predicted = predicted)
    )
}

#' @export
print.brs_fit <- function(x, ...) {
    cat("<brs_fit>\n")
    cat("  Reference samples: ", x$n_braf, " BRAF-mutant, ", x$n_ras,
        " RAS-mutant\n",
        sep = ""
    )
    cat("  Signature genes used: ", length(x$genes_used), "\n", sep = "")
    cat("  Fitted on log2(x + 1): ", isTRUE(x$log2_transform), "\n", sep = "")
    if (length(x$genes_missing) > 0) {
        cat("  Missing from expr: ", .truncate(x$genes_missing), "\n", sep = "")
    }
    if (length(x$genes_zero_variance) > 0) {
        cat("  Dropped (zero variance in reference): ",
            .truncate(x$genes_zero_variance), "\n",
            sep = ""
        )
    }
    invisible(x)
}


# Put newdata on the scale the centroids live on. "reference" uses the frozen
# per-gene mean/SD from the fit, which is what makes a sample's score
# independent of the company it is scored with. The other two re-derive the
# location and scale from newdata itself, which is what lets a different
# platform be scored at all, at the cost of that independence.
.standardize <- function(m, object, how) {
    if (how == "reference") {
        return(list(z = (m - object$gene_mean) / object$gene_sd,
                    genes = rownames(m)))
    }

    if (ncol(m) < 3L) {
        stop("`standardize = \"", how, "\"` standardizes against the scored ",
             "cohort itself and needs at least 3 samples in `newdata`; got ",
             ncol(m), ". Use the default \"reference\" for a few samples.",
             call. = FALSE)
    }

    if (how == "rank") {
        rn <- rownames(m)
        cn <- colnames(m)
        m <- apply(m, 2L, function(x) rank(x, na.last = "keep") /
                       sum(!is.na(x)))
        dimnames(m) <- list(rn, cn)
    }

    mu <- rowMeans(m, na.rm = TRUE)
    sdv <- apply(m, 1L, stats::sd, na.rm = TRUE)

    flat <- !is.finite(sdv) | sdv == 0
    if (any(flat)) {
        warning(sum(flat), " signature gene(s) have no variance across ",
                "`newdata` and were dropped from this call: ",
                .truncate(rownames(m)[flat]), call. = FALSE)
        m <- m[!flat, , drop = FALSE]
        mu <- mu[!flat]
        sdv <- sdv[!flat]
    }
    if (nrow(m) == 0L) {
        stop("No signature gene varies across `newdata`.", call. = FALSE)
    }

    list(z = (m - mu) / sdv, genes = rownames(m))
}

# Decide which scale to score on, warning if the caller contradicts the fit.
.resolve_scale <- function(object, log2_transform) {
    stored <- isTRUE(object$log2_transform)
    if (is.null(log2_transform)) {
        return(stored)
    }

    .check_flag(log2_transform, "log2_transform")
    if (!identical(as.logical(log2_transform), stored)) {
        warning("`log2_transform` = ", log2_transform,
            " disagrees with the value used to fit `object` (", stored,
            "). The stored per-gene mean/SD are on the fitted scale, ",
            "so the scores will not be meaningful.",
            call. = FALSE
        )
    }
    log2_transform
}

# Reference standardization subtracts the fit's per-gene mean and divides by
# its SD, so if `newdata` is on a different scale the z-scores blow up. A
# model fitted on log2 values and handed raw TPM gives a median |z| around 27
# and classifies confidently and wrongly, with nothing else to give it away.
.warn_if_offscale <- function(z, limit = 5) {
    typical <- stats::median(abs(z), na.rm = TRUE)
    if (is.finite(typical) && typical > limit) {
        warning("The typical signature gene in `newdata` sits ",
                round(typical, 1), " SDs from the reference mean. `newdata` ",
                "is probably not on the scale `object` was fitted on ",
                "(raw counts or TPM against log2 values, say). Check ",
                "`log2_transform`, or re-fit on the same pipeline.",
                call. = FALSE)
    }
    invisible(typical)
}

# Flag the two situations that usually mean the scores cannot be trusted.
.warn_about_scores <- function(score) {
    scored <- score[!is.na(score)]
    # Below a couple of dozen samples a one-sided cohort is unremarkable, so
    # only flag it once there are enough samples for it to be suspicious.
    if (length(scored) >= 25L && length(unique(sign(scored))) == 1L) {
        warning("All ", length(scored), " scored samples fall on the ",
            if (scored[1L] < 0) "BRAF" else "RAS",
            " side of zero. If the centroids were fitted on a different ",
            "platform or pipeline, the class threshold may not carry ",
            "over even though the ranking does; see ?containers.",
            call. = FALSE
        )
    }
    if (anyNA(score)) {
        warning(sum(is.na(score)), " of ", length(score),
            " samples could not be scored (missing values in `newdata` ",
            "for signature genes); their class is NA.",
            call. = FALSE
        )
    }
    invisible(NULL)
}

## ---- internal helpers -------------------------------------------------

.check_expr <- function(expr) {
    if (!is.matrix(expr) || !is.numeric(expr)) {
        stop("`expr`/`newdata` must be a numeric matrix (genes x samples).",
            call. = FALSE
        )
    }
    if (is.null(rownames(expr)) || is.null(colnames(expr))) {
        stop("`expr`/`newdata` must have rownames (gene symbols) and ",
            "colnames (sample IDs).",
            call. = FALSE
        )
    }
    if (anyDuplicated(colnames(expr))) {
        warning(
            "Duplicated sample identifiers in colnames; the output will ",
            "repeat them: ",
            .truncate(unique(colnames(expr)[duplicated(colnames(expr))])),
            call. = FALSE
        )
    }
}

# A fit made on one annotation vintage should still score a matrix built on
# the other. brs_fit() resolves the spelling once; predict() has to follow,
# or a model fitted on GENCODE v36 (ARNTL) rejects a current matrix (BMAL1).
.match_alias_rows <- function(mat, needed) {
    missing <- setdiff(needed, rownames(mat))
    if (!length(missing)) {
        return(mat)
    }

    partner <- c(
        stats::setNames(brs_genes$original_symbol, brs_genes$current_symbol),
        stats::setNames(brs_genes$current_symbol, brs_genes$original_symbol)
    )
    partner <- partner[!is.na(names(partner)) & !is.na(partner)]

    alt <- partner[missing]
    usable <- !is.na(alt) & alt %in% rownames(mat) & !(alt %in% needed)
    if (!any(usable)) {
        return(mat)
    }

    from <- unname(alt[usable])
    to <- missing[usable]
    rownames(mat)[match(from, rownames(mat))] <- to
    message(
        length(to), " signature gene(s) matched under the other symbol ",
        "spelling: ", paste(from, "->", to, collapse = ", ")
    )
    mat
}

.check_flag <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
        stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
    }
    invisible(TRUE)
}

.log2p1 <- function(x) {
    if (any(x < 0, na.rm = TRUE)) {
        warning("`log2_transform = TRUE` but the matrix has negative values; ",
            "it may already be on a log scale.",
            call. = FALSE
        )
    }
    log2(x + 1)
}

# Accepted synonyms for the two reference groups. The "*-like" labels are
# kept for compatibility, but the paper's groups are defined by driver
# mutation (see ?brs_fit).
.brs_label_map <- c(
    "braf_v600e" = "BRAF", "brafv600e" = "BRAF", "braf" = "BRAF",
    "braf-like" = "BRAF", "braf_like" = "BRAF", "bvl" = "BRAF",
    "ras" = "RAS", "ras-like" = "RAS", "ras_like" = "RAS", "rl" = "RAS"
)

.normalise_labels <- function(labels, samples, what = "colnames(expr)",
                              min_per_group = 2L) {
    if (is.null(names(labels))) {
        stop("`labels` must be a named vector; names are sample identifiers.",
            call. = FALSE
        )
    }
    keys <- tolower(trimws(as.character(labels)))
    out <- unname(.brs_label_map[keys])
    names(out) <- names(labels)

    out <- out[!is.na(out)]
    matched <- names(out) %in% samples
    if (any(!matched)) {
        warning(sum(!matched), " of ", length(out), " labelled sample(s) do ",
                "not appear in ", what, " and were ignored: ",
                .truncate(names(out)[!matched]),
                ". Mismatched identifiers (a barcode against a patient ID, ",
                "say) look exactly like this.", call. = FALSE)
    }
    out <- out[matched]
    if (length(out) == 0) {
        stop("No reference samples found: `labels` must name samples in ",
            what, " with a value identifying the BRAF group (e.g. ",
            "\"BRAF_V600E\") or the RAS group (e.g. \"RAS\").",
            call. = FALSE
        )
    }
    if (anyDuplicated(names(out))) {
        stop("`labels` has duplicated sample names.", call. = FALSE)
    }

    n <- table(factor(out, levels = c("BRAF", "RAS")))
    if (min_per_group > 0L && any(n < min_per_group)) {
        stop("Both reference groups need at least ", min_per_group,
            " samples to fit a centroid; got ", n[["BRAF"]],
            " BRAF-mutant and ", n[["RAS"]], " RAS-mutant.",
            call. = FALSE
        )
    }
    out
}

# Resolve each signature gene against the rows the matrix actually has,
# trying the current HGNC symbol first and the published one second. This is
# needed because the annotation vintage decides which of the two exists: a
# GENCODE v36 matrix has ARNTL, a current one has BMAL1, and neither has both.
.resolve_signature <- function(available) {
    candidates <- Map(
        function(current, original) unique(c(current, original)),
        brs_genes$current_symbol, brs_genes$original_symbol
    )
    resolved <- vapply(candidates, function(cand) {
        hit <- cand[!is.na(cand) & cand %in% available]
        if (length(hit)) hit[1L] else NA_character_
    }, character(1L))

    # A gene that appears under neither spelling is reported under its current
    # name, so that genes_missing stays readable.
    fallback <- ifelse(is.na(brs_genes$current_symbol),
        brs_genes$original_symbol, brs_genes$current_symbol
    )
    out <- ifelse(is.na(resolved), fallback, resolved)
    unique(out[!is.na(out) & !is.na(brs_genes$current_symbol)])
}

.drop_duplicate_rows <- function(mat, genes) {
    dup <- duplicated(rownames(mat))
    if (any(dup & rownames(mat) %in% genes)) {
        hit <- unique(rownames(mat)[dup & rownames(mat) %in% genes])
        warning("Duplicated row names for ", length(hit),
            " signature gene(s); keeping the first occurrence of each: ",
            .truncate(hit),
            call. = FALSE
        )
    }
    mat[!dup, , drop = FALSE]
}

# The paper's rescaling: negatives divided by |min|, positives by max, so the
# set spans [-1, 1] while zero stays the class boundary.
.rescale_brs <- function(score) {
    out <- rep(NA_real_, length(score))
    ok <- !is.na(score)
    # With a single sample there is no set to rescale against: dividing the
    # value by itself would report +-1 for any input at all.
    if (sum(ok) < 2L) return(out)
    neg <- ok & score < 0
    pos <- ok & score > 0
    out[ok & score == 0] <- 0
    if (any(neg)) out[neg] <- score[neg] / abs(min(score[neg]))
    if (any(pos)) out[pos] <- score[pos] / max(score[pos])
    out
}

.truncate <- function(x, n = 10L) {
    if (length(x) <= n) {
        return(paste(x, collapse = ", "))
    }
    paste0(
        paste(x[seq_len(n)], collapse = ", "),
        ", ... (", length(x) - n, " more)"
    )
}
