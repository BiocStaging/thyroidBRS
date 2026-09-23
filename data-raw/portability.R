## Reproduces the numbers quoted in microarray_notes.md and in the
## `standardize` documentation, which were measured but had no script:
##
##   - the resubstitution floor on pure noise (~68%)
##   - the circular vs. mutation-defined fit (axes at 0.998, 2 of 505)
##   - the signature's structure on an Affymetrix cohort
##   - the three standardizations against the array's own axis
##   - the class-mix robustness table
##
##   THCA_LOG2TPM=<rnaseq.rds> GSE33630=<array.rds> \
##       Rscript data-raw/portability.R
##
## THCA_LOG2TPM: log2(TPM + 1), gene symbols x TCGA patient barcode.
## GSE33630: log2 RMA matrix of the PTC samples, gene symbols x GSM id.
## Either may be omitted; the sections that need it are then skipped.

library(thyroidBRS)

ref <- read.csv(
    system.file("extdata", "thca_reference.csv", package = "thyroidBRS"),
    stringsAsFactors = FALSE
)
rownames(ref) <- ref$patient

## ---- 1. What resubstitution gives on data with no signal ---------------
message("\n1. Resubstitution floor, 70 genes, 391 samples split 272/119")
floor_pct <- replicate(20, {
    genes <- paste0("G", 1:70)
    samples <- paste0("S", 1:391)
    expr <- matrix(rnorm(70 * 391), 70, 391,
                   dimnames = list(genes, samples))
    labels <- setNames(rep(c("BRAF_V600E", "RAS"), c(272, 119)), samples)
    preds <- brs_score(expr, labels, genes = genes)
    suppressWarnings(validate_brs(preds, labels))$pct_concordant
})
message("   mean ", round(mean(floor_pct), 1), "%  range ",
        paste(round(range(floor_pct), 1), collapse = " - "), "%")

expr_path <- Sys.getenv("THCA_LOG2TPM")
if (nzchar(expr_path) && file.exists(expr_path)) {
    expr <- readRDS(expr_path)

    mut <- setNames(ref$driver_group, ref$patient)
    mut <- mut[mut %in% c("BRAF_V600E", "RAS") & names(mut) %in% colnames(expr)]
    fit <- brs_fit(expr, mut)

    ## ---- 2. Fitting on the published class instead of on mutation ------
    message("\n2. Circular fit vs. mutation-defined fit")
    published <- setNames(ref$published_class, ref$patient)
    published <- published[!is.na(published) & names(published) %in%
                               colnames(expr)]
    circular <- brs_fit(expr, published)

    axis_mut <- fit$centroid_ras - fit$centroid_braf
    axis_cir <- circular$centroid_ras - circular$centroid_braf
    shared <- intersect(names(axis_mut), names(axis_cir))
    message("   correlation between the two axes: ",
            round(cor(axis_mut[shared], axis_cir[shared]), 4))
    message("   tumors classified differently: ",
            sum(predict(fit, expr)$brs_class !=
                    predict(circular, expr)$brs_class),
            " of ", ncol(expr))

    ## ---- 3. Does the score depend on the cohort it is scored with? -----
    message("\n3. Class-mix robustness (same centroids, different cohorts)")
    truth <- setNames(ref$published_class, ref$patient)
    truth <- truth[!is.na(truth) & names(truth) %in% colnames(expr)]
    braf <- names(truth)[truth == "Braf-like"]
    ras <- names(truth)[truth == "Ras-like"]

    cohorts <- list(
        "all published" = c(braf, ras),
        "only Braf-like" = braf,
        "only Ras-like" = ras,
        "skewed 90/10" = c(head(braf, 180), head(ras, 20)),
        "balanced 50/50" = c(head(braf, length(ras)), ras)
    )
    for (nm in names(cohorts)) {
        cols <- cohorts[[nm]]
        acc <- vapply(c("reference", "cohort"), function(how) {
            cls <- suppressWarnings(
                predict(fit, expr[, cols, drop = FALSE], standardize = how)
            )$brs_class
            100 * mean(cls == truth[cols])
        }, numeric(1))
        message(sprintf("   %-16s n=%3d  reference %5.1f%%  cohort %5.1f%%",
                        nm, length(cols), acc[1], acc[2]))
    }

    ## ---- 4. Carrying the centroids to a microarray ---------------------
    array_path <- Sys.getenv("GSE33630")
    if (nzchar(array_path) && file.exists(array_path)) {
        arr <- readRDS(array_path)
        shared <- intersect(fit$genes_used, rownames(arr))
        arr <- arr[shared, ]
        message("\n4. GSE33630: ", length(shared), " of ",
                length(fit$genes_used), " signature genes present")

        array_fit <- brs_fit(expr, mut, genes = shared)

        ## The array's own signature axis, oriented like the BRS: higher
        ## means more RAS-like, so block-1 genes should load positively.
        z <- t(scale(t(arr)))
        pr <- prcomp(t(z), center = FALSE)
        pc <- pr$x[, 1]
        ld <- pr$rotation[, 1]

        up_ras <- c(setNames(brs_genes$up_in, brs_genes$current_symbol),
                    setNames(brs_genes$up_in, brs_genes$original_symbol))
        ras_rows <- up_ras[rownames(arr)] == "RAS"
        if (mean(ld[ras_rows]) < 0) {
            pc <- -pc
            ld <- -ld
        }

        var1 <- pr$sdev[1]^2 / sum(pr$sdev^2)
        message("   PC1 explains ", round(100 * var1, 1), "% of the signature")
        message("   genes loading with their block's sign: ",
                sum(sign(ld) == ifelse(ras_rows, 1, -1)), " of ", length(ld))

        for (how in c("reference", "cohort", "rank")) {
            p <- suppressWarnings(predict(array_fit, arr, standardize = how))
            message(sprintf(
                "   %-10s spearman vs PC1 %.3f | %2d Braf / %2d Ras",
                            how, cor(p$brs_score, pc, method = "spearman"),
                            sum(p$brs_class == "Braf-like"),
                            sum(p$brs_class == "Ras-like")))
        }
    } else {
        message("\n4. GSE33630 not set; skipping the microarray section")
    }
} else {
    message("\nTHCA_LOG2TPM not set; skipping sections 2-4")
}
