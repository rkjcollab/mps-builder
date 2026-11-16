dat <- load_test_data()

test_that("plot_model_comparison draws continuous and dichotomous plots", {
  results <- list(comparison = list(comparison_table = data.frame(
    Outcome = c("Continuous", "Continuous", "Dichotomous", "Dichotomous"),
    Model = c("Elastic Net", "Random Forest", "Elastic Net", "Random Forest"),
    R2_Mean = c(0.4, 0.3, NA, NA), R2_SD = c(0.1, 0.1, NA, NA),
    AUC_Mean = c(NA, NA, 0.7, 0.6), AUC_SD = c(NA, NA, 0.05, 0.05)
  )))
  withr::local_pdf(NULL)
  expect_no_error(suppressWarnings(plot_model_comparison(results)))
})

test_that("export_results writes CpG, performance, and score files", {
  folds <- create_nested_folds(dat$exposure, k_outer = 3, k_inner = 3)
  en <- quiet(elastic_net_nested_cv(dat$beta, dat$exposure,
    nested_folds = folds, alpha_range = 1
  ))
  rf <- quiet(random_forest_nested_cv_tuned(dat$beta, dat$exposure, folds,
    ntree = 50, n_cores = 1,
    hyperparam_grid = expand.grid(
      mtry = 10, min.node.size = 5, sample.fraction = 0.8
    )
  ))
  results <- list(
    continuous = list(
      elastic_net = en,
      random_forest = rf,
      elastic_net_score =
        create_mps(dat$beta, en, "elastic_net"),
      random_forest_score =
        create_mps(dat$beta, rf, "random_forest")
    ),
    comparison = list(comparison_table = rbind(en$performance, rf$performance))
  )
  out_dir <- withr::local_tempdir()

  quiet(export_results(results, dat$beta, out_dir))

  en_cpgs <- read.csv(file.path(out_dir, "elastic_net_cpgs.csv"))
  expect_equal(en_cpgs$CpG, en$selected_cpgs)
  rf_cpgs <- read.csv(file.path(out_dir, "random_forest_cpgs.csv"))
  expect_equal(nrow(rf_cpgs), length(rf$avg_importance))
  perf <- read.csv(file.path(out_dir, "model_performance.csv"))
  expect_equal(nrow(perf), 2)
  scores <- read.csv(file.path(out_dir, "profile_scores.csv"))
  expect_equal(scores$Sample_ID, colnames(dat$beta))
})
