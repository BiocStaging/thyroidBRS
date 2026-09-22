## Re-check the BRS signature symbols against current annotation.
## Run after an org.Hs.eg.db update; it prints what, if anything, has drifted.
##
##   Rscript data-raw/check_symbols.R

library(org.Hs.eg.db)
devtools::load_all(".", quiet = TRUE)

approved <- AnnotationDbi::keys(org.Hs.eg.db, "SYMBOL")

resolve <- function(symbol) {
    hit <- tryCatch(
        AnnotationDbi::select(org.Hs.eg.db, keys = symbol, keytype = "ALIAS",
                              columns = "SYMBOL"),
        error = function(e) NULL
    )
    if (is.null(hit)) return(NA_character_)
    paste(unique(hit$SYMBOL), collapse = "/")
}

cat("org.Hs.eg.db: ", as.character(packageVersion("org.Hs.eg.db")), "\n",
    "checked on:   ", format(Sys.Date()), "\n\n", sep = "")

## 1. Every current_symbol we ship should be an approved symbol today.
stale <- with(brs_genes, current_symbol[!is.na(current_symbol) &
                                        !current_symbol %in% approved])
if (length(stale)) {
    cat("current_symbol no longer approved:\n")
    for (s in stale) cat(sprintf("  %-12s -> %s\n", s, resolve(s)))
} else {
    cat("All current_symbol values are approved HGNC symbols.\n")
}

## 2. Every original_symbol that is stale should have a mapping (or be known
##    to be unresolvable).
cat("\nOriginal symbols that are not approved symbols today:\n")
for (s in with(brs_genes, original_symbol[!original_symbol %in% approved])) {
    shipped <- brs_genes$current_symbol[brs_genes$original_symbol == s]
    cat(sprintf("  %-12s org.Hs.eg.db: %-10s  shipped: %s\n",
                s, resolve(s), ifelse(is.na(shipped), "<NA>", shipped)))
}

## 3. Shape invariants.
stopifnot(
    nrow(brs_genes) == 71L,
    !anyDuplicated(brs_genes$original_symbol),
    sum(is.na(brs_genes$current_symbol)) == 1L,
    identical(as.integer(table(brs_genes$block)), c(13L, 58L)),
    !is.unsorted(brs_genes$original_symbol[brs_genes$block == 1L]),
    !is.unsorted(brs_genes$original_symbol[brs_genes$block == 2L])
)
cat("\nShape invariants OK.\n")
