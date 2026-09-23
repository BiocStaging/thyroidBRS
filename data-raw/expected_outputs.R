## Assisted-by: Claude Opus 5 (Anthropic). See the Provenance section
## of README.md.
##
## Writes, and on later runs checks, the reference outputs behind the figures
## quoted in README.md. Nothing else does: the test suite runs on synthetic
## data, and the figures come from TCGA and GEO matrices too large to ship.
## A number that silently drifts has, until now, had nothing to drift against.
##
##   THCA_LOG2TPM=<rds> Rscript data-raw/expected_outputs.R          # check
##   THCA_LOG2TPM=<rds> BRS_WRITE_EXPECTED=1 Rscript ...             # rewrite
##
## Writes data-raw/expected/preds_tcga.csv and data-raw/expected/inputs.md5.
## The md5 line matters as much as the predictions: if the input hash moved,
## a changed prediction means the input changed, not the classifier.

library(thyroidBRS)

source(file.path("data-raw", "record_session.R"))

expr_path <- Sys.getenv("THCA_LOG2TPM", "thca_log2tpm.rds")
if (!file.exists(expr_path)) {
    stop("Set THCA_LOG2TPM; see data-raw/build_inputs.R for how to build it.")
}
write_mode <- nzchar(Sys.getenv("BRS_WRITE_EXPECTED", ""))

dir <- file.path("data-raw", "expected")
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
preds_file <- file.path(dir, "preds_tcga.csv")
md5_file <- file.path(dir, "inputs.md5")

expr <- readRDS(expr_path)
ref <- read.csv(
    file.path("inst", "extdata", "thca_reference.csv"),
    stringsAsFactors = FALSE
)

mut <- setNames(ref$driver_group, ref$patient)
mut <- mut[mut %in% c("BRAF_V600E", "RAS") & names(mut) %in% colnames(expr)]

fit <- brs_fit(expr, mut)
preds <- predict(fit, expr)
preds[, c("brs_score", "brs_scaled")] <-
    round(preds[, c("brs_score", "brs_scaled")], 8)

hashes <- data.frame(
    file = c(basename(expr_path), "thca_reference.csv"),
    md5 = c(tools::md5sum(expr_path),
            tools::md5sum(file.path("inst", "extdata",
                                    "thca_reference.csv"))),
    row.names = NULL
)

if (write_mode) {
    write.csv(preds, preds_file, row.names = FALSE)
    write.csv(hashes, md5_file, row.names = FALSE)
    message("wrote ", preds_file, " (", nrow(preds), " samples) and ",
            md5_file)
} else {
    if (!file.exists(preds_file) || !file.exists(md5_file)) {
        stop("No reference output yet; rerun with BRS_WRITE_EXPECTED=1.")
    }
    want <- read.csv(preds_file, stringsAsFactors = FALSE)
    want_md5 <- read.csv(md5_file, stringsAsFactors = FALSE)

    ## Compare the same tumors, not the same row numbers. A GDC refresh or a
    ## change to which aliquot represents a patient alters the sample set,
    ## and comparing positionally would then contrast different tumors and
    ## call it a pass -- in exactly the case this script exists to catch.
    if (!identical(want$sample, preds$sample)) {
        gone <- setdiff(want$sample, preds$sample)
        new_s <- setdiff(preds$sample, want$sample)
        stop("the sample set differs from the record: ",
             length(gone), " missing, ", length(new_s), " new, ",
             "order ", if (setequal(want$sample, preds$sample))
                 "changed" else "aside", ". Rebuild the record if that is ",
             "intended (BRS_WRITE_EXPECTED=1).")
    }

    same_input <- isTRUE(all.equal(hashes$md5, want_md5$md5))
    if (!same_input) {
        message("NOTE: the inputs differ from the ones on record, so any ",
                "difference below is the input's, not the classifier's")
        print(merge(want_md5, hashes, by = "file",
                    suffixes = c("_recorded", "_now")))
    }

    diff_class <- sum(want$brs_class != preds$brs_class)
    max_shift <- max(abs(want$brs_score - preds$brs_score))
    message("compared ", nrow(preds), " samples against the record")
    message("  class changes: ", diff_class)
    message("  largest score shift: ", signif(max_shift, 3))

    if (diff_class > 0L || max_shift > 1e-6) {
        stop("Predictions no longer match data-raw/expected/preds_tcga.csv")
    }
    message("  identical")
}

record_session("expected_outputs")
