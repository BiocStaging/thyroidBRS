#' BRAF-RAS Score (BRS) signature genes
#'
#' The 71-gene expression signature used to compute the BRAF-RAS Score (BRS)
#' for papillary thyroid carcinoma (PTC), as published in Agrawal et al.
#' (2014), Figure S7A. Gene symbols were transcribed manually from the
#' published heatmap (no machine-readable gene list was published alongside
#' the figure) and cross-checked against current HGNC aliases.
#'
#' @format A data frame with 71 rows and 2 columns:
#' \describe{
#'   \item{original_symbol}{Gene symbol exactly as it appears in the 2014
#'     publication's Figure S7A.}
#'   \item{current_symbol}{Current HGNC-approved symbol. Identical to
#'     `original_symbol` for genes whose symbol has not changed. `NA` for
#'     `FLJ23867`, which could not be resolved to any symbol present in
#'     current genome annotations (checked against GENCODE-derived gene
#'     symbols and `org.Hs.eg.db` aliases) and should be treated as
#'     permanently missing from the signature — [brs_fit()] and
#'     [predict.brs_fit()] silently drop it.}
#' }
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
#'
#' @export
brs_genes <- data.frame(
  original_symbol = c(
    "ANKRD46", "CYB561", "GNA14", "HGD", "KATNAL2", "KCNAB1", "KCNIP3", "LGI3", "MLEC", "NQO1",
    "SFTPC", "SLC4A4", "SORBS2", "ABTB2", "AHR", "ANKLE2", "ANXA1", "ANXA2P2", "ARNTL", "ASAP2",
    "BID", "CDC42EP1", "COL8A2", "CREB5", "CTSC", "CYP1B1", "DTX4", "DUSP5", "ETHE1", "FAM176A",
    "FAM20C", "FCHO1", "FLJ23867", "FN1", "FSTL3", "GABRB2", "GBP2", "ITGA3", "ITGB8", "KCNN4",
    "LAMB3", "LLGL1", "LY6E", "MDFIC", "MET", "PDE5A", "PDLIM4", "PLCD3", "PLEKHA6", "PNPLA5",
    "PPL", "PRICKLE1", "PROS1", "PTPRE", "PVRL4", "RASGEF1B", "RUNX1", "RUNX2", "SDC4", "SEL1L3",
    "SFTPB", "SLC35F2", "SOX4", "SPOCK2", "STK17B", "SYT12", "TACSTD2", "TBC1D2", "TGFBR1",
    "TM7SF4", "TMEM43"
  ),
  stringsAsFactors = FALSE
)

.brs_alias_to_current <- c(FAM176A = "EVA1A", PVRL4 = "NECTIN4", TM7SF4 = "DCSTAMP")

brs_genes$current_symbol <- ifelse(
  brs_genes$original_symbol %in% names(.brs_alias_to_current),
  .brs_alias_to_current[brs_genes$original_symbol],
  brs_genes$original_symbol
)
brs_genes$current_symbol[brs_genes$original_symbol == "FLJ23867"] <- NA_character_
