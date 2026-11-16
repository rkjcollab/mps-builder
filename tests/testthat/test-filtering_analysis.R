dat <- load_test_data()
filtered <- quiet(filter_cpg_sites(dat$beta, dat$exposure,
  literature_studies = dat$literature_studies, n_cores = 1
))

test_that("plot_filtering_venn draws without error", {
  withr::local_dir(withr::local_tempdir()) # VennDiagram writes log files here
  withr::local_pdf(NULL)
  expect_no_error(quiet(plot_filtering_venn(filtered)))
})

test_that("sensitivity_analysis_filtering runs every parameter combination", {
  res <- quiet(sensitivity_analysis_filtering(dat$beta, dat$exposure,
    dat$literature_studies,
    var_percentiles = c(70, 90), pval_thresholds = c(0.01, 0.05),
    min_criteria_options = c(1, 2)
  ))

  expect_equal(nrow(res), 8)
  expect_false(anyNA(res$n_selected))
  # loosening min_criteria never selects fewer CpGs
  one <- res[res$min_criteria == 1, "n_selected"]
  two <- res[res$min_criteria == 2, "n_selected"]
  expect_true(all(one >= two))
})

test_that("compare_temporal_filtering compares CpGs across exposure windows", {
  set.seed(4)
  periods <- list(
    trim1 = dat$exposure,
    trim2 = dat$exposure + rnorm(length(dat$exposure), 0, 0.5)
  )
  res <- quiet(
    compare_temporal_filtering(dat$beta, periods, dat$literature_studies)
  )

  expect_named(res$by_period, c("trim1", "trim2"))
  expect_equal(dim(res$overlap_matrix), c(2, 2))
  expect_equal(unname(diag(res$overlap_matrix)), unname(res$summary$N_Selected))
  expect_setequal(
    res$consistent_cpgs,
    intersect(
      res$by_period$trim1$selected_cpgs, res$by_period$trim2$selected_cpgs
    )
  )
})

test_that("analyze_literature_cpgs summarises literature CpGs", {
  res <- quiet(analyze_literature_cpgs(filtered, dat$beta, dat$exposure))

  lit_cpgs <- unique(unlist(dat$literature_studies))
  expect_setequal(res$literature_cpg_details$Literature_CpG, lit_cpgs)
  expect_equal(res$summary$total_literature_cpgs, length(lit_cpgs))
  expect_equal(
    res$literature_cpg_details$Correlation,
    unname(sapply(res$literature_cpg_details$Literature_CpG, function(cpg) {
      cor(dat$beta[cpg, ], dat$exposure)
    }))
  )
})

test_that("analyze_literature_cpgs returns NULL without literature studies", {
  no_lit <- quiet(filter_cpg_sites(dat$beta, dat$exposure, n_cores = 1))
  expect_null(quiet(analyze_literature_cpgs(no_lit, dat$beta, dat$exposure)))
})
