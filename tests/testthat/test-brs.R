make_synthetic_cohort <- function(n_genes = 20, n_braf = 30, n_ras = 30, n_unlabeled = 10, seed = 1) {
  set.seed(seed)
  genes <- paste0("GENE", seq_len(n_genes))
  n_total <- n_braf + n_ras + n_unlabeled
  samples <- paste0("S", seq_len(n_total))
  true_class <- c(rep("Braf-like", n_braf), rep("Ras-like", n_ras), rep("Ras-like", n_unlabeled))

  base_mean <- ifelse(true_class == "Braf-like", 2, 6)
  expr <- matrix(
    stats::rnorm(n_genes * n_total, mean = rep(base_mean, each = n_genes), sd = 0.3),
    nrow = n_genes, dimnames = list(genes, samples)
  )

  labels <- setNames(true_class[seq_len(n_braf + n_ras)], samples[seq_len(n_braf + n_ras)])

  list(expr = expr, labels = labels, true_class = setNames(true_class, samples))
}

test_that("brs_genes has the expected shape", {
  expect_equal(nrow(brs_genes), 71)
  expect_equal(sum(is.na(brs_genes$current_symbol)), 1)
  expect_equal(brs_genes$original_symbol[is.na(brs_genes$current_symbol)], "FLJ23867")
})

test_that("brs_fit recovers reference labels with high concordance", {
  cohort <- make_synthetic_cohort()
  fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))

  expect_s3_class(fit, "brs_fit")
  expect_equal(fit$n_braf, 30)
  expect_equal(fit$n_ras, 30)
  expect_length(fit$genes_used, 20)

  preds <- predict(fit, cohort$expr[, names(cohort$labels)])
  expect_equal(mean(preds$brs_class == cohort$labels[preds$sample]), 1)
})

test_that("predict.brs_fit scores previously-unlabeled samples correctly", {
  cohort <- make_synthetic_cohort()
  fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))

  unlabeled_samples <- setdiff(colnames(cohort$expr), names(cohort$labels))
  preds <- predict(fit, cohort$expr[, unlabeled_samples, drop = FALSE])

  expect_equal(nrow(preds), length(unlabeled_samples))
  expect_equal(preds$brs_class, unname(cohort$true_class[preds$sample]))
})

test_that("brs_score is equivalent to fit + predict", {
  cohort <- make_synthetic_cohort()
  one_shot <- brs_score(cohort$expr, cohort$labels, genes = rownames(cohort$expr))
  fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))
  two_step <- predict(fit, cohort$expr)
  expect_equal(one_shot, two_step)
})

test_that("validate_brs reports concordance against known labels", {
  cohort <- make_synthetic_cohort()
  preds <- brs_score(cohort$expr, cohort$labels, genes = rownames(cohort$expr))
  result <- validate_brs(preds, cohort$labels)

  expect_equal(result$n, length(cohort$labels))
  expect_equal(result$pct_concordant, 100)
})

test_that("predict.brs_fit errors clearly when newdata is missing a signature gene", {
  cohort <- make_synthetic_cohort()
  fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))

  incomplete <- cohort$expr[-1, , drop = FALSE]
  expect_error(predict(fit, incomplete), "missing")
})

test_that("zero-variance genes in the reference are dropped, not fatal", {
  cohort <- make_synthetic_cohort()
  cohort$expr["GENE1", names(cohort$labels)] <- 5  # constant across reference samples

  fit <- brs_fit(cohort$expr, cohort$labels, genes = rownames(cohort$expr))
  expect_true("GENE1" %in% fit$genes_zero_variance)
  expect_false("GENE1" %in% fit$genes_used)
})

test_that("brs_fit errors when no reference labels match", {
  cohort <- make_synthetic_cohort()
  bad_labels <- setNames(cohort$labels, paste0("nope_", names(cohort$labels)))
  expect_error(brs_fit(cohort$expr, bad_labels, genes = rownames(cohort$expr)), "No reference samples")
})
