# Provenance: substantial parts of this file were written by Claude Opus 5
# (Anthropic) working from the published description of the method, under the
# direction of and reviewed by the package authors, who are responsible for
# its correctness and maintenance. See the Provenance section of README.md.
# Assisted-by: Claude Opus 5 (Anthropic)

#' BRAF-RAS Score (BRS) signature genes
#'
#' The 71-gene expression signature used to compute the BRAF-RAS Score (BRS)
#' for papillary thyroid carcinoma (PTC), as published in Agrawal et al.
#' (2014), Figure S7A. No machine-readable gene list was released with the
#' paper, so the symbols were transcribed from the published heatmap and then
#' verified against it row by row. The provenance and the checks that were
#' run are written up in `data-raw/brs_genes.md`, which is in the package
#' repository rather than the installed package:
#' <https://github.com/camilasauria/thyroidBRS/blob/main/data-raw/brs_genes.md>.
#'
#' @format A data frame with 71 rows and 4 columns:
#' \describe{
#'   \item{original_symbol}{Gene symbol exactly as it appears in the 2014
#'     publication's Figure S7A.}
#'   \item{current_symbol}{Current HGNC-approved symbol, resolved against
#'     `org.Hs.eg.db`. Identical to `original_symbol` for genes whose symbol
#'     has not changed. `NA` for `FLJ23867`, which does not map to any symbol
#'     or alias in current annotation and should be treated as permanently
#'     missing from the signature. [brs_fit()] and [predict.brs_fit()]
#'     drop it.}
#'   \item{block}{Which of the two row-clusters of Figure S7A the gene belongs
#'     to: `1` for the 13-gene block, `2` for the 58-gene block. Each block is
#'     listed alphabetically in the figure, which is what makes the
#'     transcription checkable: no gene is out of order within its block.}
#'   \item{up_in}{Which reference group the gene is higher in, `"RAS"` for
#'     block 1 and `"BRAF"` for block 2. Refitting the centroids on TCGA-THCA
#'     recovers this split exactly: all 13 block-1 genes have a higher RAS
#'     centroid and all 57 usable block-2 genes a higher BRAF centroid, with
#'     no exceptions. That is an independent check on the gene list: a
#'     mis-transcribed symbol would land in the wrong block.}
#' }
#'
#' @section Alias resolution:
#' Four symbols in the 2014 figure are no longer the approved HGNC symbol:
#' `FAM176A` (now `EVA1A`), `PVRL4` (`NECTIN4`), `TM7SF4` (`DCSTAMP`) and
#' `ARNTL` (`BMAL1`). Expression matrices built from current GENCODE/HGNC
#' annotation use the new symbols, so matching on `original_symbol` silently
#' loses those genes. Use `current_symbol`, which is what [brs_fit()] does by
#' default. Resolution was checked against `org.Hs.eg.db`; re-run
#' `data-raw/check_symbols.R` when the annotation is updated. `data-raw/` is
#' not shipped in the installed package; it lives in the package repository,
#' <https://github.com/camilasauria/thyroidBRS>.
#'
#' @source Agrawal N, Akbani R, Aksoy BA, et al. "Integrated Genomic
#'   Characterization of Papillary Thyroid Carcinoma." Cell.
#'   2014;159(3):676-690. \doi{10.1016/j.cell.2014.09.050}, Figure S7A
#'   ("mRNA expression of the 71 genes used to derive the BRAF-V600E-RAS
#'   score (BRS) across 391 PTC samples").
#'
#' @examples
#' brs_genes
#' # Genes actually usable today (alias-resolved, FLJ23867 excluded):
#' na.omit(brs_genes$current_symbol)
#' # The two blocks of Figure S7A:
#' table(brs_genes$block)
#'
#' @export
brs_genes <- data.frame(
    original_symbol = c(
        # Block 1 of Figure S7A (13 genes, alphabetical)
        "ANKRD46", "CYB561", "GNA14", "HGD", "KATNAL2", "KCNAB1", "KCNIP3",
        "LGI3", "MLEC", "NQO1", "SFTPC", "SLC4A4", "SORBS2",
        # Block 2 of Figure S7A (58 genes, alphabetical)
        "ABTB2", "AHR", "ANKLE2", "ANXA1", "ANXA2P2", "ARNTL", "ASAP2",
        "BID", "CDC42EP1", "COL8A2", "CREB5", "CTSC", "CYP1B1", "DTX4",
        "DUSP5", "ETHE1", "FAM176A", "FAM20C", "FCHO1", "FLJ23867", "FN1",
        "FSTL3", "GABRB2", "GBP2", "ITGA3", "ITGB8", "KCNN4", "LAMB3",
        "LLGL1", "LY6E", "MDFIC", "MET", "PDE5A", "PDLIM4", "PLCD3",
        "PLEKHA6", "PNPLA5", "PPL", "PRICKLE1", "PROS1", "PTPRE", "PVRL4",
        "RASGEF1B", "RUNX1", "RUNX2", "SDC4", "SEL1L3", "SFTPB", "SLC35F2",
        "SOX4", "SPOCK2", "STK17B", "SYT12", "TACSTD2", "TBC1D2", "TGFBR1",
        "TM7SF4", "TMEM43"
    ),
    block = rep(c(1L, 2L), c(13L, 58L)),
    up_in = rep(c("RAS", "BRAF"), c(13L, 58L)),
    stringsAsFactors = FALSE
)

.brs_alias_to_current <- c(
    ARNTL = "BMAL1",
    FAM176A = "EVA1A",
    PVRL4 = "NECTIN4",
    TM7SF4 = "DCSTAMP"
)

# FLJ23867 resolves to no current symbol or alias, so it is marked NA and
# brs_fit() drops it from the signature rather than looking for it in vain.
.brs_unresolved <- "FLJ23867"

brs_genes$current_symbol <- ifelse(
    brs_genes$original_symbol %in% names(.brs_alias_to_current),
    .brs_alias_to_current[brs_genes$original_symbol],
    brs_genes$original_symbol
)
brs_genes$current_symbol[
    brs_genes$original_symbol %in% .brs_unresolved
] <- NA_character_

brs_genes <- brs_genes[, c(
    "original_symbol", "current_symbol", "block",
    "up_in"
)]
