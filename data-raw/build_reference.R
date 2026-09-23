## Assisted-by: Claude Opus 5 (Anthropic). See the Provenance section
## of README.md.
##
## Builds inst/extdata/thca_reference.csv from the supplementary table of
## Agrawal et al. (2014), so the reference labels the validation runs on have
## a stated origin rather than appearing in the repository from nowhere.
##
##   THCA_MMC3=<path to mmc3.xlsx> Rscript data-raw/build_reference.R
##
## THCA_MMC3 is supplementary Table S1 of
##
##   Agrawal N, Akbani R, Aksoy BA, et al. Integrated Genomic
##   Characterization of Papillary Thyroid Carcinoma. Cell.
##   2014;159(3):676-690. doi:10.1016/j.cell.2014.09.050
##
## distributed as `1-s2.0-S0092867414012380-mmc3.xlsx`. It is copyrighted and
## is NOT redistributed with this package; download it from the publisher.
##
## Sheet "THCA-TP (496) ", one row per primary tumor. Columns taken:
##
##   sample            TCGA aliquot barcode
##   BRAFV600E_RAS  -> driver_group
##   exome          -> has_exome
##   BRAF_RAF_score -> published_brs
##   BRAF_RAF_class -> published_class
##
## On `driver_group`. The published column has three levels and the third is
## a catch-all: `OTHER` is everything that is neither BRAF-V600E-mutant nor
## RAS-mutant, which mixes two things that matter to keep apart. Of the 210
## OTHER tumors, 94 were never exome sequenced at all, and 116 were. Among
## the 116, those with a published BRS break down as 58 with a fusion driver,
## 10 with another driver mutation, and 43 with no identified driver — the
## "dark matter" of the paper. `has_exome` is carried through so the first
## distinction can be made; the second needs the driver columns of the
## supplement, which this package does not ship.
##
## This matters when reading the validation: the tumors held out of the
## centroid fit are the 111 OTHER tumors that have a published BRS, and all
## of them were exome sequenced. None of the 94 unsequenced tumors is in
## that set.

suppressPackageStartupMessages(library(readxl))

mmc3 <- Sys.getenv("THCA_MMC3", "1-s2.0-S0092867414012380-mmc3.xlsx")
if (!file.exists(mmc3)) {
    stop("Set THCA_MMC3 to the supplementary spreadsheet; see the header.")
}

raw <- as.data.frame(read_excel(
    mmc3,
    sheet = "THCA-TP (496) ",
    col_names = TRUE,
    .name_repair = "minimal"
))

blank_to_na <- function(x) {
    x <- as.character(x)
    x[!nzchar(trimws(x))] <- NA_character_
    x
}

ref <- data.frame(
    sample = as.character(raw$sample),
    patient = substr(as.character(raw$sample), 1, 12),
    driver_group = blank_to_na(raw$BRAFV600E_RAS),
    has_exome = as.integer(suppressWarnings(as.numeric(raw$exome))) == 1L,
    published_brs = suppressWarnings(as.numeric(raw$BRAF_RAF_score)),
    published_class = blank_to_na(raw$BRAF_RAF_class),
    stringsAsFactors = FALSE
)

## Invariants worth failing on rather than discovering later.
stopifnot(
    "expected 496 primary tumors" = nrow(ref) == 496L,
    "patient identifiers should be unique" = !anyDuplicated(ref$patient),
    "driver_group should have no missing values" = !anyNA(ref$driver_group),
    "published_class should be the sign of published_brs" =
        identical(
            ifelse(ref$published_brs < 0, "Braf-like", "Ras-like"),
            ref$published_class
        ),
    "every tumor with a published BRS should have been exome sequenced" =
        all(ref$has_exome[!is.na(ref$published_brs)])
)

out <- file.path("inst", "extdata", "thca_reference.csv")
write.csv(ref, out, row.names = FALSE)

message("wrote ", out, ": ", nrow(ref), " tumors")
message("  driver_group: ",
        paste(names(table(ref$driver_group)), table(ref$driver_group),
              sep = "=", collapse = "  "))
message("  with a published BRS: ", sum(!is.na(ref$published_brs)))
message("  exome sequenced: ", sum(ref$has_exome))
