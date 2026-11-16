dat <- load_test_data()

test_that("create_literature_indicator flags CpGs per study and overall", {
  cpgs <- c("cg1", "cg2", "cg3", "cg4")
  studies <- list(A = c("cg1", "cg2"), B = c("cg2", "cg9"))

  res <- create_literature_indicator(cpgs, studies)

  expect_equal(res$overall_indicator, c(1L, 1L, 0L, 0L))
  expect_equal(unname(res$by_study[, "A"]), c(1L, 1L, 0L, 0L))
  expect_equal(unname(res$by_study[, "B"]), c(0L, 1L, 0L, 0L))
  expect_equal(unname(res$n_studies), c(1, 2, 0, 0))
  expect_setequal(res$literature_cpgs, c("cg1", "cg2", "cg9"))
})

test_that("filter_cpg_sites keeps CpGs meeting at least min_criteria", {
  res <- quiet(filter_cpg_sites(dat$beta, dat$exposure,
    literature_studies = dat$literature_studies,
    var_percentile = 80, pval_threshold = 0.05, min_criteria = 2, n_cores = 1
  ))

  expect_equal(dim(res$criteria_matrix), c(nrow(dat$beta), 3))
  expect_equal(
    colnames(res$criteria_matrix),
    c("Variance", "Association", "Literature")
  )
  expect_equal(res$criteria_counts, rowSums(res$criteria_matrix))
  expect_setequal(res$selected_cpgs, names(which(res$criteria_counts >= 2)))
  expect_equal(res$n_selected, length(res$selected_cpgs))
  expect_equal(rownames(res$filtered_beta), res$selected_cpgs)
  expect_equal(colnames(res$filtered_beta), colnames(dat$beta))

  # variance criterion: top 20% of CpGs by variance
  expect_equal(sum(res$criteria_matrix[, "Variance"]), 20)
  # association criterion matches the stored p-values
  expect_equal(
    unname(res$criteria_matrix[, "Association"]),
    as.numeric(res$p_values < 0.05)
  )
  # literature criterion matches the literature lists
  lit_cpgs <- unique(unlist(dat$literature_studies))
  expect_equal(
    unname(res$criteria_matrix[, "Literature"]),
    as.numeric(rownames(dat$beta) %in% lit_cpgs)
  )
  expect_equal(res$filters_summary$n_total, res$n_selected)
})

test_that("filter_cpg_sites uses cor.test for a continuous exposure", {
  res <- quiet(filter_cpg_sites(dat$beta, dat$exposure, n_cores = 1))
  expected <- cor.test(dat$beta[5, ], dat$exposure)$p.value
  expect_equal(unname(res$p_values[5]), expected)
})

test_that("filter_cpg_sites uses a t-test for a dichotomous exposure", {
  res <- quiet(filter_cpg_sites(dat$beta, dat$exposure_binary, n_cores = 1))
  expected <- t.test(dat$beta[5, ] ~ dat$exposure_binary)$p.value
  expect_equal(unname(res$p_values[5]), expected)
})

test_that("filter_cpg_sites works without literature studies", {
  res <- quiet(filter_cpg_sites(dat$beta, dat$exposure,
    min_criteria = 1, n_cores = 1
  ))
  expect_null(res$literature_info)
  expect_true(all(res$criteria_matrix[, "Literature"] == 0))
})

test_that("filter_cpg_sites parallel and sequential p-values match", {
  # The parallel branch only runs for more than 1000 CpGs and n_cores > 1
  set.seed(3)
  big <- matrix(runif(1100 * 20), nrow = 1100, dimnames = list(
    sprintf("cg%04d", 1:1100), sprintf("S%02d", 1:20)
  ))
  y <- rnorm(20)

  seq_res <- quiet(filter_cpg_sites(big, y, min_criteria = 1, n_cores = 1))
  par_res <- quiet(filter_cpg_sites(big, y, min_criteria = 1, n_cores = 2))

  expect_equal(par_res$p_values, seq_res$p_values)
})

test_that("analyze_filtering_results prints a report", {
  res <- quiet(filter_cpg_sites(dat$beta, dat$exposure,
    literature_studies = dat$literature_studies, n_cores = 1
  ))
  expect_output(
    analyze_filtering_results(res, dat$beta), "FILTERING ANALYSIS REPORT"
  )
  expect_output(analyze_filtering_results(res, dat$beta), "OVERLAP ANALYSIS")
})

test_that("export_filtering_results writes criteria, CpGs, and summary", {
  res <- quiet(filter_cpg_sites(dat$beta, dat$exposure,
    literature_studies = dat$literature_studies, n_cores = 1
  ))
  out_dir <- withr::local_tempdir()

  quiet(export_filtering_results(res, out_dir))

  criteria <- read.csv(file.path(out_dir, "cpg_filtering_criteria.csv"))
  expect_equal(nrow(criteria), nrow(dat$beta))
  expect_true(
    all(c("StudyA", "StudyB", "StudyC", "N_Studies") %in% names(criteria))
  )

  selected <- readLines(file.path(out_dir, "selected_cpgs.txt"))
  expect_setequal(selected, res$selected_cpgs)

  expect_true(file.exists(file.path(out_dir, "filtering_summary.csv")))
})
