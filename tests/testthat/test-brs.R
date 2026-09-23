make_synthetic_cohort <- function(n_genes = 20, n_braf = 30, n_ras = 30,
                                  n_unlabeled = 10, seed = 1) {
    set.seed(seed)
    genes <- paste0("GENE", seq_len(n_genes))
    n_total <- n_braf + n_ras + n_unlabeled
    samples <- paste0("S", seq_len(n_total))
    true_class <- c(
        rep("Braf-like", n_braf),
        rep("Ras-like", n_ras + n_unlabeled)
    )

    base_mean <- ifelse(true_class == "Braf-like", 2, 6)
    expr <- matrix(
        stats::rnorm(n_genes * n_total,
            mean = rep(base_mean, each = n_genes),
            sd = 0.3
        ),
        nrow = n_genes, dimnames = list(genes, samples)
    )

    n_ref <- n_braf + n_ras
    labels <- setNames(
        rep(c("BRAF_V600E", "RAS"), c(n_braf, n_ras)),
        samples[seq_len(n_ref)]
    )

    list(
        expr = expr, labels = labels,
        true_class = setNames(true_class, samples)
    )
}

test_that("brs_genes has the expected shape", {
    expect_equal(nrow(brs_genes), 71)
    expect_equal(sum(is.na(brs_genes$current_symbol)), 1)
    expect_equal(
        brs_genes$original_symbol[is.na(brs_genes$current_symbol)],
        "FLJ23867"
    )
    expect_false(anyDuplicated(brs_genes$original_symbol) > 0)
})

test_that("brs_genes records the two blocks of Figure S7A", {
    expect_equal(as.integer(table(brs_genes$block)), c(13L, 58L))
    # Each block is alphabetical in the figure: this is the check that the
    # transcription lost no gene and reordered none. Compared against a radix
    # sort, which orders in the C locale on every platform — is.unsorted()
    # would compare in whatever collation the session happens to run under.
    for (b in c(1L, 2L)) {
        sym <- brs_genes$original_symbol[brs_genes$block == b]
        expect_identical(sym, sort(sym, method = "radix"))
    }
})

test_that("stale symbols are resolved, ARNTL included", {
    map <- setNames(brs_genes$current_symbol, brs_genes$original_symbol)
    expect_equal(unname(map[["ARNTL"]]), "BMAL1")
    expect_equal(unname(map[["PVRL4"]]), "NECTIN4")
    expect_equal(unname(map[["FAM176A"]]), "EVA1A")
    expect_equal(unname(map[["TM7SF4"]]), "DCSTAMP")
})

test_that("brs_fit recovers reference labels with high concordance", {
    cohort <- make_synthetic_cohort()
    fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))

    expect_s3_class(fit, "brs_fit")
    expect_equal(fit$n_braf, 30)
    expect_equal(fit$n_ras, 30)
    expect_length(fit$genes_used, 20)

    preds <- predict(fit, cohort$expr[, names(cohort$labels)])
    expect_equal(mean(preds$brs_class == cohort$true_class[preds$sample]), 1)
})

test_that("predict.brs_fit scores previously-unlabeled samples correctly", {
    cohort <- make_synthetic_cohort()
    fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))

    unlabeled <- setdiff(colnames(cohort$expr), names(cohort$labels))
    preds <- predict(fit, cohort$expr[, unlabeled, drop = FALSE])

    expect_equal(nrow(preds), length(unlabeled))
    expect_equal(preds$brs_class, unname(cohort$true_class[preds$sample]))
})

test_that("brs_score is equivalent to fit + predict", {
    cohort <- make_synthetic_cohort()
    one_shot <- brs_score(cohort$expr, cohort$labels,
        genes = rownames(cohort$expr)
    )
    fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))
    two_step <- predict(fit, cohort$expr)
    expect_equal(one_shot, two_step)
})

test_that("the score follows the published formula", {
    cohort <- make_synthetic_cohort()
    fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))
    preds <- predict(fit, cohort$expr)

    z <- (cohort$expr[fit$genes_used, ] - fit$gene_mean) / fit$gene_sd
    p <- length(fit$genes_used)
    d_b <- sqrt(colSums((z - fit$centroid_braf)^2) / p)
    d_r <- sqrt(colSums((z - fit$centroid_ras)^2) / p)

    expect_equal(preds$brs_score, unname(d_b - d_r))
})

test_that("brs_scaled reproduces the published [-1, 1] rescaling", {
    cohort <- make_synthetic_cohort()
    preds <- brs_score(cohort$expr, cohort$labels,
        genes = rownames(cohort$expr)
    )

    expect_true(all(preds$brs_scaled >= -1 & preds$brs_scaled <= 1))
    # Exactly one sample at each end, as in the published table.
    expect_equal(sum(preds$brs_scaled == -1), 1L)
    expect_equal(sum(preds$brs_scaled == 1), 1L)
    # Rescaling preserves the sign, and therefore the class.
    expect_equal(sign(preds$brs_scaled), sign(preds$brs_score))
})

test_that("validate_brs reports concordance against known labels", {
    cohort <- make_synthetic_cohort()
    fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))
    preds <- predict(fit, cohort$expr)
    result <- suppressWarnings(validate_brs(preds, cohort$labels, fit = fit))

    expect_equal(result$n, length(cohort$labels))
    expect_equal(result$pct_concordant, 100)
})

test_that("validate_brs flags a pure resubstitution estimate", {
    cohort <- make_synthetic_cohort()
    fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))
    preds <- predict(fit, cohort$expr)

    expect_warning(
        validate_brs(preds, cohort$labels, fit = fit),
        "resubstitution"
    )
    res <- suppressWarnings(validate_brs(preds, cohort$labels, fit = fit))
    expect_equal(res$n_resubstituted, length(cohort$labels))
})

test_that("predict.brs_fit errors clearly when newdata is missing a gene", {
    cohort <- make_synthetic_cohort()
    fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))

    incomplete <- cohort$expr[-1, , drop = FALSE]
    expect_error(predict(fit, incomplete), "missing")
})

test_that("zero-variance genes in the reference are dropped, not fatal", {
    cohort <- make_synthetic_cohort()
    cohort$expr["GENE1", names(cohort$labels)] <- 5

    fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))
    expect_true("GENE1" %in% fit$genes_zero_variance)
    expect_false("GENE1" %in% fit$genes_used)
})

test_that("brs_fit errors when no reference labels match", {
    cohort <- make_synthetic_cohort()
    bad <- setNames(cohort$labels, paste0("nope_", names(cohort$labels)))
    expect_error(
        brs_fit(cohort$expr, bad, genes = rownames(cohort$expr)),
        "No reference samples"
    )
})

test_that("brs_fit accepts both the mutation and the '*-like' vocabulary", {
    cohort <- make_synthetic_cohort()
    like <- setNames(
        ifelse(cohort$labels == "BRAF_V600E", "Braf-like", "Ras-like"),
        names(cohort$labels)
    )
    a <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))
    b <- brs_fit(cohort$expr, like, genes = rownames(cohort$expr))
    expect_equal(a$centroid_braf, b$centroid_braf)
    expect_equal(a$n_ras, b$n_ras)
})

## ---- B1: log2_transform can no longer drift --------------------

test_that("predict reuses the log2_transform the model was fitted with", {
    cohort <- make_synthetic_cohort()
    tpm <- 2^cohort$expr - 1

    fit <- brs_fit(tpm, cohort$labels,
        genes = rownames(tpm),
        log2_transform = TRUE
    )
    expect_true(fit$log2_transform)

    preds <- predict(fit, tpm)
    expect_equal(mean(preds$brs_class == cohort$true_class[preds$sample]), 1)
})

test_that("overriding log2_transform against the fitted scale warns", {
    cohort <- make_synthetic_cohort()
    tpm <- 2^cohort$expr - 1
    fit <- brs_fit(tpm, cohort$labels,
        genes = rownames(tpm),
        log2_transform = TRUE
    )

    expect_warning(predict(fit, tpm, log2_transform = FALSE), "disagrees")
})

## ---- B2: a single reference class -------------------------------------

test_that("brs_fit refuses to fit when a reference group is too small", {
    cohort <- make_synthetic_cohort()
    one_class <- cohort$labels[cohort$labels == "BRAF_V600E"]
    expect_error(
        brs_fit(cohort$expr, one_class, genes = rownames(cohort$expr)),
        "at least 2 samples"
    )

    lopsided <- c(
        cohort$labels[cohort$labels == "BRAF_V600E"],
        cohort$labels[cohort$labels == "RAS"][1]
    )
    expect_error(
        brs_fit(cohort$expr, lopsided, genes = rownames(cohort$expr)),
        "at least 2 samples"
    )
})

## ---- B3: missing signature genes --------------------------------------

test_that("brs_fit warns when signature genes are missing from expr", {
    cohort <- make_synthetic_cohort()
    genes <- c(rownames(cohort$expr), "NOT_MEASURED")
    expect_warning(
        brs_fit(cohort$expr, cohort$labels, genes = genes),
        "NOT_MEASURED"
    )
})

## ---- B4/B5: NA y rownames duplicados ----------------------------------

test_that("NA in newdata warns instead of yielding a silent NA class", {
    cohort <- make_synthetic_cohort()
    fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))

    broken <- cohort$expr
    broken["GENE1", "S1"] <- NA
    expect_warning(preds <- predict(fit, broken), "could not be scored")
    expect_true(is.na(preds$brs_class[preds$sample == "S1"]))
    expect_equal(sum(is.na(preds$brs_class)), 1L)
})

test_that("duplicated signature rows warn and resolve to the first", {
    cohort <- make_synthetic_cohort()
    fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))

    dup <- rbind(cohort$expr, cohort$expr["GENE1", , drop = FALSE] * 0)
    expect_warning(preds <- predict(fit, dup), "Duplicated row names")
    expect_equal(preds$brs_score, predict(fit, cohort$expr)$brs_score)
})

## ---- B6: empate exacto ------------------------------------------------

test_that("a score of exactly zero is classified as Ras-like", {
    cohort <- make_synthetic_cohort()
    fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))

    # Exact midpoint between the two centroids.
    midpoint <- (fit$centroid_braf + fit$centroid_ras) / 2
    tie <- as.matrix(midpoint * fit$gene_sd + fit$gene_mean)
    colnames(tie) <- "tie"

    preds <- predict(fit, tie)
    expect_equal(preds$brs_score, 0)
    expect_equal(preds$brs_class, "Ras-like")
    expect_equal(preds$brs_scaled, 0)
})

test_that("the signature resolves against either annotation vintage", {
    lab <- setNames(rep(c("BRAF_V600E", "RAS"), each = 20), paste0("S", 1:40))
    mk <- function(rn) {
        set.seed(3)
        matrix(rnorm(length(rn) * 40, 5), length(rn), 40,
            dimnames = list(rn, paste0("S", 1:40))
        )
    }
    modern <- na.omit(brs_genes$current_symbol) # BMAL1, NECTIN4
    v36 <- brs_genes$original_symbol[!is.na(brs_genes$current_symbol)]

    # GENCODE v36 has ARNTL/PVRL4; current annotation has BMAL1/NECTIN4. All
    # 70 genes must resolve in either case, and in a mixed matrix.
    expect_length(brs_fit(mk(modern), lab)$genes_used, 70L)
    expect_length(brs_fit(mk(v36), lab)$genes_used, 70L)
    mixed <- mk(c(head(modern, 35), tail(v36, 35)))
    expect_length(brs_fit(mixed, lab)$genes_used, 70L)
})

## ---- container support -------------------------------------------------

test_that("brs_fit accepts a SummarizedExperiment and picks the assay", {
    skip_if_not_installed("SummarizedExperiment")
    cohort <- make_synthetic_cohort()

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = 2^cohort$expr - 1, logtpm = cohort$expr),
        colData = data.frame(
            driver = cohort$true_class,
            row.names = colnames(cohort$expr)
        )
    )

    fit <- brs_fit(se, cohort$labels,
        genes = rownames(cohort$expr),
        assay = "logtpm"
    )
    ref <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))
    expect_equal(fit$centroid_braf, ref$centroid_braf)

    # More than one assay and no choice made: first one, with a warning.
    expect_warning(
        brs_fit(se, cohort$labels, genes = rownames(cohort$expr)),
        "none was chosen"
    )
    expect_error(brs_fit(se, cohort$labels,
        genes = rownames(cohort$expr),
        assay = "nope"
    ), "No assay named")
})

test_that("labels can name a colData column", {
    skip_if_not_installed("SummarizedExperiment")
    cohort <- make_synthetic_cohort()
    driver <- ifelse(cohort$true_class == "Braf-like", "BRAF_V600E", "RAS")

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(logtpm = cohort$expr),
        colData = data.frame(
            driver = driver,
            row.names = colnames(cohort$expr)
        )
    )

    fit <- brs_fit(se, "driver", genes = rownames(cohort$expr))
    expect_equal(fit$n_braf + fit$n_ras, ncol(cohort$expr))
    expect_error(
        brs_fit(se, "missing_col", genes = rownames(cohort$expr)),
        "No column"
    )
    expect_error(
        brs_fit(cohort$expr, "driver", genes = rownames(cohort$expr)),
        "named vector"
    )
})

test_that("brs_fit accepts an ExpressionSet", {
    skip_if_not_installed("Biobase")
    cohort <- make_synthetic_cohort()
    es <- Biobase::ExpressionSet(assayData = cohort$expr)

    fit <- brs_fit(es, cohort$labels, genes = rownames(cohort$expr))
    ref <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))
    expect_equal(fit$centroid_ras, ref$centroid_ras)
    expect_equal(predict(fit, es)$brs_score,
                 predict(ref, cohort$expr)$brs_score)
})

test_that("predict warns when every sample lands on one side", {
    cohort <- make_synthetic_cohort()
    fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))

    braf_only <- cohort$expr[, names(cohort$labels)[
        cohort$labels == "BRAF_V600E"
    ], drop = FALSE]
    expect_warning(predict(fit, braf_only), "side of zero")
})

## ---- standardize -------------------------------------------------------

test_that("standardize defaults to reference and is unchanged", {
    cohort <- make_synthetic_cohort()
    fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))

    expect_equal(predict(fit, cohort$expr),
                 predict(fit, cohort$expr, standardize = "reference"))
    expect_error(predict(fit, cohort$expr, standardize = "nope"))
})

# A cohort whose two blocks move in opposite directions, as the real
# signature does: 10 genes up in BRAF, 10 up in RAS.
make_two_block_cohort <- function(seed = 4) {
    set.seed(seed)
    genes <- paste0("GENE", 1:20)
    samples <- paste0("S", 1:60)
    cls <- rep(c("BRAF_V600E", "RAS"), each = 30)
    shift <- rep(c(-1, 1), each = 10)

    expr <- matrix(rnorm(20 * 60, 5), 20, 60, dimnames = list(genes, samples))
    expr <- expr + outer(shift, ifelse(cls == "RAS", 1, -1))
    list(expr = expr, labels = setNames(cls, samples))
}

test_that("every standardization ranks samples the same way", {
    cohort <- make_two_block_cohort()
    fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))

    s <- lapply(c("reference", "cohort", "rank"), function(how) {
        predict(fit, cohort$expr, standardize = how)$brs_score
    })
    # The axis is what carries over; see ?predict.brs_fit.
    expect_gt(cor(s[[1]], s[[2]], method = "spearman"), 0.95)
    expect_gt(cor(s[[1]], s[[3]], method = "spearman"), 0.90)
})

test_that("rank standardization discards a signature-wide shift", {
    # If every signature gene moves the same way between the groups, the
    # difference lives entirely in the overall level and within-sample ranks
    # cannot see it. That is the intended behaviour, and the reason "rank" is
    # the weakest of the three on real data too.
    cohort <- make_synthetic_cohort()
    fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))

    ref <- predict(fit, cohort$expr)$brs_score
    rk <- predict(fit, cohort$expr, standardize = "rank")$brs_score
    expect_lt(abs(cor(ref, rk, method = "spearman")), 0.5)
})

test_that("cohort and rank refuse to score a handful of samples", {
    cohort <- make_synthetic_cohort()
    fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))
    two <- cohort$expr[, 1:2, drop = FALSE]

    expect_error(predict(fit, two, standardize = "cohort"), "at least 3")
    expect_error(predict(fit, two, standardize = "rank"), "at least 3")
    expect_s3_class(predict(fit, two), "data.frame")
})

test_that("reference scoring does not depend on the rest of the cohort", {
    cohort <- make_synthetic_cohort()
    fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))
    one <- "S65"

    alone <- predict(fit, cohort$expr[, one, drop = FALSE])$brs_score
    within <- predict(fit, cohort$expr)
    expect_equal(alone, within$brs_score[within$sample == one])

    # Under "cohort" it does depend on it — that is the trade-off.
    a <- predict(fit, cohort$expr[, 61:70], standardize = "cohort")
    b <- predict(fit, cohort$expr, standardize = "cohort")
    expect_false(isTRUE(all.equal(
        a$brs_score[a$sample == one],
        b$brs_score[b$sample == one]
    )))
})

test_that("genes flat across newdata are dropped with a warning", {
    cohort <- make_synthetic_cohort()
    fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))

    flat <- cohort$expr
    flat["GENE1", ] <- 7
    expect_warning(predict(fit, flat, standardize = "cohort"), "no variance")
})

test_that("brs_score passes standardize through", {
    cohort <- make_synthetic_cohort()
    expect_equal(
        brs_score(cohort$expr, cohort$labels, genes = rownames(cohort$expr),
                  standardize = "cohort"),
        predict(brs_fit(cohort$expr, cohort$labels,
                        genes = rownames(cohort$expr)),
                cohort$expr, standardize = "cohort")
    )
})
