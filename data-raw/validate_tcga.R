## Validate thyroidBRS against the published TCGA-THCA BRAF-RAS Score.
##
## The point of this script is to measure agreement with something that did
## NOT go into the fit. The centroids are fitted on driver-mutation status
## (BRAF-V600E vs RAS); the published `BRAF_RAF_score` / `BRAF_RAF_class`
## are only ever used to check the result. Fitting on the published class
## instead is circular — it is the sign of the published score.
##
##   Rscript data-raw/validate_tcga.R
##
## Expects a log2(TPM + 1) matrix of TCGA-THCA primary tumors, gene symbols as
## row names and patient barcodes (TCGA-XX-XXXX) as column names. Build it
## from the GDC "STAR - Counts" files, e.g.
##
##   library(TCGAbiolinks)
##   q <- GDCquery(project = "TCGA-THCA",
##                 data.category = "Transcriptome Profiling",
##                 data.type = "Gene Expression Quantification",
##                 workflow.type = "STAR - Counts")
##   GDCdownload(q)
##   se <- GDCprepare(q)
##   # keep assay "tpm_unstranded", primary tumors (-01), collapse duplicated
##   # symbols, then log2(x + 1).

## Assisted-by: Claude Opus 5 (Anthropic). See the Provenance section
## of README.md.

library(thyroidBRS)

expr_path <- Sys.getenv("THCA_LOG2TPM", "thca_log2tpm.rds")
if (!file.exists(expr_path)) {
    stop("Set THCA_LOG2TPM to a log2(TPM + 1) matrix; see the header.")
}
expr <- readRDS(expr_path)

ref <- read.csv(
    system.file("extdata", "thca_reference.csv", package = "thyroidBRS"),
    stringsAsFactors = FALSE
)
rownames(ref) <- ref$patient

## ---- fit on mutation status -------------------------------------------
mut <- setNames(ref$driver_group, ref$patient)
mut <- mut[mut %in% c("BRAF_V600E", "RAS") & names(mut) %in% colnames(expr)]

fit <- brs_fit(expr, mut)
print(fit)

preds <- predict(fit, expr)
rownames(preds) <- preds$sample

## ---- validate against the published score -----------------------------
pub <- ref[ref$patient %in% preds$sample & !is.na(ref$published_brs), ]
i <- pub$patient
held_out <- !(i %in% fit$reference_samples)

cat("\nAgreement with the published BRS (n = ", nrow(pub), ")\n", sep = "")
cat("  Spearman vs BRAF_RAF_score   : ",
    round(cor(preds[i, "brs_score"], pub$published_brs,
              method = "spearman"), 4), "\n", sep = "")
cat("  Pearson, rescaled            : ",
    round(cor(preds[i, "brs_scaled"], pub$published_brs), 4), "\n", sep = "")
cat("  Class concordance, all       : ",
    round(100 * mean(preds[i, "brs_class"] == pub$published_class), 1),
    "%\n", sep = "")
cat("  Class concordance, held out  : ",
    round(100 * mean(preds[i, "brs_class"][held_out] ==
                     pub$published_class[held_out]), 1),
    "% (n = ", sum(held_out), ")\n", sep = "")
print(table(published = pub$published_class, predicted = preds[i, "brs_class"]))

## ---- cross-validate on the reference set ------------------------------
set.seed(1)
fold <- sample(rep(seq_len(10), length.out = length(mut)))
cv <- character(length(mut))
names(cv) <- names(mut)
for (k in seq_len(10)) {
    f <- brs_fit(expr, mut[fold != k])
    in_fold <- names(mut)[fold == k]
    cv[in_fold] <- predict(f, expr[, in_fold, drop = FALSE])$brs_class
}
truth <- ifelse(mut == "BRAF_V600E", "Braf-like", "Ras-like")
cat("\n10-fold CV against mutation status (n = ", length(mut), "): ",
    round(100 * mean(cv == truth), 1), "%\n", sep = "")

## ---- coverage gained ---------------------------------------------------
new <- setdiff(preds$sample, pub$patient)
cat("\nScored ", nrow(preds), " patients; ", length(new),
    " have no published BRS.\n", sep = "")
print(table(preds[new, "brs_class"]))

## ---- orientation check -------------------------------------------------
## Block 1 of Figure S7A should be higher in RAS, block 2 higher in BRAF.
delta <- fit$centroid_ras - fit$centroid_braf
key <- c(setNames(brs_genes$up_in, brs_genes$current_symbol),
         setNames(brs_genes$up_in, brs_genes$original_symbol))
expected <- key[names(delta)]
cat("\nOrientation matches Figure S7A for ",
    sum(ifelse(delta > 0, "RAS", "BRAF") == expected), " of ", length(delta),
    " genes\n", sep = "")
