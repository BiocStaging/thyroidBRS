## Assisted-by: Claude Opus 5 (Anthropic). See the Provenance section
## of README.md.
##
## Every figure quoted in the README and in data-raw/*.md came out of a
## particular set of package versions. Without them a number that later fails
## to reproduce cannot be told apart from a number that was wrong.
##
## Each script in data-raw/ calls record_session() when it finishes, writing
## data-raw/sessioninfo/<name>.txt. Those files are committed.

record_session <- function(name) {
    dir <- file.path("data-raw", "sessioninfo")
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    path <- file.path(dir, paste0(name, ".txt"))

    con <- file(path, "wt")
    on.exit(close(con))
    writeLines(c(
        paste("script:", name),
        paste("run on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
        ""
    ), con)
    writeLines(capture.output(utils::sessionInfo()), con)

    message("recorded environment in ", path)
    invisible(path)
}
