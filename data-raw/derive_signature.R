## Re-derive the 71-gene BRS signature from TCGA-THCA, following the
## procedure in Agrawal et al. (2014), Extended Experimental Procedures 14.1:
##
##   "First, we iteratively sub-sampled equal size sets of tumors with either
##    BRAF-V600E (B set) or RAS (R set) activating mutations. At each
##    iteration, we determined the set of differentially expressed genes
##    (q < 0.01) between B and R using the limma R package with voom
##    correction. We then selected the genes that consistently ranked at each
##    iteration among the most significant 100. This gave us a list of 71
##    genes."
##
## This is the strongest available check on the gene list, which otherwise
## rests on a transcription of a heatmap figure.
##
##   THCA_COUNTS=<derive_input.rds> Rscript data-raw/derive_signature.R
##
## THCA_COUNTS is derive_input.rds, the list build_inputs.R writes:
##
##   counts  raw gene counts, gene symbols x patient barcode
##   grp     named character vector of driver status, "BRAF_V600E" / "RAS"
##   keep    gene symbols surviving filterByExpr, the universe to test
##   coding  gene symbols that are protein-coding, for BRS_UNIVERSE=coding
##
## Following an earlier version of this header, which described a different
## structure, left `groups` NULL and the script failed on the first draw.

## Assisted-by: Claude Opus 5 (Anthropic). See the Provenance section
## of README.md.

suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
    library(parallel)
})
library(thyroidBRS)

source(file.path("data-raw", "record_session.R"))

N_ITER <- as.integer(Sys.getenv("BRS_ITER", "200"))
## mclapply() forks, so it does no work in parallel on Windows: set
## BRS_CORES=1 there. The result does not depend on the core count -- each
## iteration seeds from its own index -- only the wall clock does.
N_CORES <- as.integer(Sys.getenv("BRS_CORES", "10"))
if (.Platform$OS.type == "windows" && N_CORES > 1L) {
    message("mclapply cannot fork on Windows; falling back to one core")
    N_CORES <- 1L
}
TOP_N <- 100L
Q_CUT <- 0.01

input <- readRDS(Sys.getenv("THCA_COUNTS", "derive_input.rds"))
stopifnot(
    "THCA_COUNTS must be a list; see the header" = is.list(input),
    "THCA_COUNTS needs `counts`, `grp` and `keep`; see the header" =
        all(c("counts", "grp", "keep") %in% names(input))
)
counts <- input$counts
groups <- input$grp
universe <- input$keep

## The 2014 analysis ran on the TCGA RSEM gene model (~20.5k genes, almost
## all protein-coding). Current STAR counts add ~16k lncRNAs and pseudogenes
## that compete for the top-100 slots and could not have been picked in 2014.
## BRS_UNIVERSE=coding restricts to that older universe.
if (identical(Sys.getenv("BRS_UNIVERSE"), "coding") && !is.null(input$coding)) {
    universe <- intersect(universe, input$coding)
}

message("universe: ", length(universe), " genes; ",
        sum(groups == "BRAF_V600E"), " BRAF-V600E, ",
        sum(groups == "RAS"), " RAS")

braf <- names(groups)[groups == "BRAF_V600E"]
ras <- names(groups)[groups == "RAS"]

## One iteration: equal-size draw from each group, limma-voom, top TOP_N of
## the genes passing q < Q_CUT.
one_iteration <- function(seed, n_per_group) {
    set.seed(seed)
    picked <- c(sample(braf, n_per_group), sample(ras, n_per_group))
    grp <- factor(groups[picked], levels = c("BRAF_V600E", "RAS"))

    dge <- calcNormFactors(DGEList(counts[universe, picked]))
    design <- model.matrix(~grp)
    fit <- eBayes(lmFit(voom(dge, design), design))

    tt <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
    head(rownames(tt)[tt$adj.P.Val < Q_CUT], TOP_N)
}

run <- function(n_per_group, label) {
    message("\n", label, ": ", N_ITER, " iterations of ", n_per_group,
            " vs ", n_per_group)
    hits <- mclapply(seq_len(N_ITER), one_iteration,
                     n_per_group = n_per_group, mc.cores = N_CORES)
    freq <- sort(table(unlist(hits)), decreasing = TRUE) / N_ITER

    ## "consistently ranked at each iteration among the most significant 100"
    always <- names(freq)[freq == 1]
    message("  genes in the top ", TOP_N, " at EVERY iteration: ",
            length(always))
    for (cut in c(0.99, 0.95, 0.90, 0.75, 0.50)) {
        message("  at >= ", 100 * cut, "% of iterations: ",
                sum(freq >= cut))
    }
    list(freq = freq, always = always, hits = hits)
}

## If we simply take the N most consistent genes, how many are the published
## ones? This is the comparison that does not depend on guessing a threshold.
top_n_overlap <- function(freq, published, n = length(published)) {
    top <- names(freq)[seq_len(min(n, length(freq)))]
    hit <- intersect(top, published)
    message("\n  top ", n, " genes by consistency: ", length(hit), " of ",
            length(published), " published (",
            round(100 * length(hit) / length(published)), "%)")
    setdiff(top, published)
}

## The paper does not say how many iterations it ran, and the count drives
## the result directly: "in the top 100 every time" is a weaker filter over
## 10 draws than over 200. Show the whole curve.
iteration_curve <- function(hits, published) {
    message("\n  iterations -> genes kept at 100% (of which published)")
    for (k in c(5, 10, 20, 50, 100, length(hits))) {
        if (k > length(hits)) next
        tab <- table(unlist(hits[seq_len(k)]))
        kept <- names(tab)[tab == k]
        message("    ", formatC(k, width = 4), " -> ",
                formatC(length(kept), width = 4), "  (",
                length(intersect(kept, published)), " published, ",
                round(100 * length(intersect(kept, published)) /
                      max(1, length(kept))), "% precision)")
    }
}

## The published signature, spelled the way this matrix spells it.
published <- ifelse(is.na(brs_genes$current_symbol),
                    brs_genes$original_symbol, brs_genes$current_symbol)
published <- ifelse(published %in% rownames(counts), published,
                    brs_genes$original_symbol)

compare <- function(derived, label) {
    both <- intersect(derived, published)
    message("\n", label, " vs the published 71:")
    message("  derived: ", length(derived), "   published: ", length(published))
    message("  recovered: ", length(both), " (",
            round(100 * length(both) / length(published), 1),
            "% of the published list)")
    missed <- setdiff(published, derived)
    extra <- setdiff(derived, published)
    if (length(missed)) message("  published but not derived: ",
                                paste(missed, collapse = ", "))
    if (length(extra)) message("  derived but not published: ",
                               paste(extra, collapse = ", "))
    invisible(list(both = both, missed = missed, extra = extra))
}

main <- run(min(length(braf), length(ras)), "equal-size sets (all RAS)")
compare(main$always, "Top-100-always")
iteration_curve(main$hits, published)
extra71 <- top_n_overlap(main$freq, published)
message("  not published: ", paste(head(extra71, 25), collapse = ", "))

## Sensitivity: sub-sample both groups, so neither set is fixed.
sens <- run(40L, "sub-sampling both groups (40 vs 40)")
compare(sens$always, "Top-100-always, 40 vs 40")
iteration_curve(sens$hits, published)
top_n_overlap(sens$freq, published)

saveRDS(list(main = main, sens = sens, published = published),
        Sys.getenv("BRS_OUT", "derive_signature_result.rds"))

record_session("derive_signature")
