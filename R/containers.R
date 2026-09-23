#' Expression container support
#'
#' `brs_fit()`, `predict()` and `brs_score()` accept expression as a plain
#' numeric matrix, a [SummarizedExperiment::SummarizedExperiment] or a
#' `Biobase::ExpressionSet`. Whichever you pass, the signature is matched on
#' row names, so those must be gene symbols.
#'
#' @section Choosing an assay:
#' A `SummarizedExperiment` usually carries several assays (for GDC RNA-seq:
#' `unstranded`, `tpm_unstranded`, `fpkm_unstranded`, ...). Pass the one you
#' want by name through `assay`; the first assay is used if you do not, which
#' for GDC data is raw counts and is almost certainly not what you want.
#' An `ExpressionSet` has a single `exprs()` matrix and ignores `assay`.
#'
#' @section Labels from the object:
#' When `labels` is a bare character string of length one it is looked up as a
#' column of `colData()` (or `pData()`), so `brs_fit(se, "driver_mutation")`
#' works directly on an object that carries its own sample annotation. Sample
#' names come from `colnames()`. A *named* vector of length one is a label for
#' one sample, not a column name.
#'
#' Whichever column you name, it has to hold driver-mutation status. On
#' TCGA-THCA the column giving the published BRAF-like/RAS-like class is the
#' classifier's own output, and fitting on it is circular; see the warning in
#' [brs_fit()].
#'
#' @section Microarray data:
#' Microarray expression works as input: the signature is well represented
#' on common arrays and its structure is preserved there. What does not
#' carry over from one platform to another is the *zero threshold* that
#' separates
#' the two classes, because the per-gene means and SDs stored in a `brs_fit`
#' are those of the platform it was fitted on. Fit centroids on the same
#' platform you are scoring whenever you have driver-mutation labels for it.
#' Scoring array data against RNA-seq centroids still ranks samples correctly,
#' but the class call is not trustworthy; `predict()` warns when every sample
#' lands on the same side, which is the usual symptom.
#'
#' @name containers
#' @keywords internal
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom methods is
NULL


# Pull a genes x samples numeric matrix out of whatever the user passed.
.as_expr_matrix <- function(x, assay = NULL, arg = "expr") {
    if (is.matrix(x)) {
        if (!is.null(assay)) {
            warning("`assay` is ignored when `", arg, "` is a matrix.",
                call. = FALSE
            )
        }
        return(x)
    }

    if (methods::is(x, "SummarizedExperiment")) {
        assays <- SummarizedExperiment::assayNames(x)
        if (is.null(assay)) {
            assay <- if (length(assays)) assays[1L] else 1L
            if (length(assays) > 1L) {
                warning("`", arg, "` has ", length(assays),
                    " assays and none was chosen; using \"", assay,
                    "\". Pass `assay=` to pick another (",
                    .truncate(assays, 6L), ").",
                    call. = FALSE
                )
            }
        } else if (is.character(assay) && !assay %in% assays) {
            stop("No assay named \"", assay, "\" in `", arg, "`; available: ",
                .truncate(assays, 10L),
                call. = FALSE
            )
        }
        m <- SummarizedExperiment::assay(x, assay)
        return(.finish_matrix(m, rownames(x), colnames(x), arg))
    }

    if (methods::is(x, "ExpressionSet")) {
        .need("Biobase")
        if (!is.null(assay)) {
            warning("`assay` is ignored for an ExpressionSet.", call. = FALSE)
        }
        m <- Biobase::exprs(x)
        return(.finish_matrix(m, rownames(m), colnames(m), arg))
    }

    stop("`", arg, "` must be a numeric matrix, a SummarizedExperiment or an ",
        "ExpressionSet; got ", class(x)[1L], ".",
        call. = FALSE
    )
}

.finish_matrix <- function(m, rn, cn, arg) {
    m <- as.matrix(m)
    if (is.null(rownames(m)) && !is.null(rn)) rownames(m) <- rn
    if (is.null(colnames(m)) && !is.null(cn)) colnames(m) <- cn
    if (!is.numeric(m)) {
        stop("The selected assay of `", arg, "` is not numeric.",
            call. = FALSE
        )
    }
    m
}

# `labels` given as a column name: pull it off colData()/pData().
.labels_from_object <- function(x, labels) {
    # A single string means "look this up as a colData()/pData() column", but
    # only when it is bare: a named vector of length one is a label for one
    # sample, and reading it as a column name is how that used to fail.
    if (!is.character(labels) || length(labels) != 1L ||
        !is.null(names(labels))) {
        return(labels)
    }
    if (is.matrix(x)) {
        stop("`labels` must be a named vector when `expr` is a matrix; a ",
            "single string is only read as a column name of a ",
            "SummarizedExperiment or ExpressionSet.",
            call. = FALSE
        )
    }

    pheno <- if (methods::is(x, "SummarizedExperiment")) {
        as.data.frame(SummarizedExperiment::colData(x))
    } else {
        .need("Biobase")
        Biobase::pData(x)
    }

    if (!labels %in% colnames(pheno)) {
        stop("No column \"", labels, "\" in the sample metadata of `expr`; ",
            "available: ", .truncate(colnames(pheno), 10L),
            call. = FALSE
        )
    }
    stats::setNames(as.character(pheno[[labels]]), colnames(x))
}

.need <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Package \"", pkg, "\" is needed for this input type. ",
            "Install it with BiocManager::install(\"", pkg, "\").",
            call. = FALSE
        )
    }
    invisible(TRUE)
}
