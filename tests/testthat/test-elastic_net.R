dat <- load_test_data()
folds <- create_nested_folds(dat$exposure, k_outer = 3, k_inner = 3)

test_that("elastic_net_nested_cv runs for continuous exposure and covariates", {
  res <- quiet(elastic_net_nested_cv(dat$beta, dat$exposure,
    covariates = dat$covariates, nested_folds = folds, alpha_range = c(0.5, 1)
  ))

  expect_equal(
    names(res$performance),
    c("Model", "RMSE_Mean", "RMSE_SD", "R2_Mean", "R2_SD", "Avg_Features")
  )
  expect_length(res$outer_fold_results, 3)
  for (fold in res$outer_fold_results) {
    expect_true(fold$best_alpha %in% c(0.5, 1))
    expect_named(fold$covariate_coefficients, names(dat$covariates))
  }

  # CpG coefficients only, one per CpG; covariates are unpenalized
  expect_equal(names(res$coefficients), rownames(dat$beta))
  expect_equal(res$penalty_factors, c(rep(1, nrow(dat$beta)), 0, 0))
  expect_named(res$covariate_coefficients, names(dat$covariates))

  # selection frequency counts outer folds; stable CpGs are in >= 80% of them
  expect_true(all(res$cpg_selection_frequency %in% 1:3))
  expect_true(all(res$cpg_selection_frequency[res$stable_cpgs] >= 3 * 0.8))
  expect_s3_class(res$final_model, "glmnet")
})

test_that("elastic_net_nested_cv runs for a dichotomous exposure", {
  bin_folds <- create_nested_folds(dat$exposure_binary,
    k_outer = 3, k_inner = 3
  )
  res <- quiet(elastic_net_nested_cv(dat$beta, dat$exposure_binary,
    nested_folds = bin_folds, alpha_range = c(0.5, 1)
  ))

  perf <- res$performance
  expect_true(all(c("Accuracy_Mean", "F1_Mean", "AUC_Mean") %in% names(perf)))
  expect_true(perf$Accuracy_Mean >= 0 && perf$Accuracy_Mean <= 1)
  expect_null(res$covariate_coefficients)
})

test_that("elastic_net_nested_cv errors when samples do not match", {
  expect_error(
    elastic_net_nested_cv(dat$beta, dat$exposure[-1], nested_folds = folds),
    "doesn't match"
  )
  expect_error(
    elastic_net_nested_cv(dat$beta, dat$exposure,
      covariates = dat$covariates[-1, ], nested_folds = folds
    ),
    "covariates must match"
  )
})
