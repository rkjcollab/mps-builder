dat <- load_test_data()
folds <- create_nested_folds(dat$exposure, k_outer = 3, k_inner = 3)
# The default grid has 96 combinations; two keep the tests fast
small_grid <- data.frame(
  eta = c(0.1, 0.3), max_depth = 2, subsample = 1, colsample_bytree = 1,
  min_child_weight = 1, gamma = 0
)

test_that("xgboost_nested_cv runs for a continuous exposure", {
  res <- quiet(xgboost_nested_cv(dat$beta, dat$exposure, folds,
    hyperparam_grid = small_grid, nrounds_max = 20, early_stopping_rounds = 5,
    n_cores = 1, importance_threshold = 10, verbose = FALSE
  ))

  perf_names <- names(res$performance)
  expect_true(all(c("RMSE_Mean", "R2_Mean", "Cor_Mean") %in% perf_names))
  expect_length(res$outer_fold_results, 3)
  expect_equal(nrow(res$fold_hyperparameters), 3)
  expect_true(res$consensus_params$eta %in% c(0.1, 0.3))
  expect_true(res$consensus_nrounds >= 1)
  expect_true(length(res$top_cpgs) <= 10)
  expect_true(all(res$top_cpgs %in% rownames(dat$beta)))
  expect_s3_class(res$final_model, "xgb.Booster")
})

test_that("xgboost_nested_cv keeps covariates out of CpG importance", {
  res <- quiet(xgboost_nested_cv(dat$beta, dat$exposure, folds,
    covariates = dat$covariates, hyperparam_grid = small_grid,
    nrounds_max = 20, early_stopping_rounds = 5, n_cores = 1, verbose = FALSE
  ))
  # covariates are features but are kept out of the CpG importance
  expect_false(any(names(dat$covariates) %in% names(res$avg_importance)))
})

test_that("xgboost_nested_cv runs for a dichotomous exposure", {
  bin_folds <- create_nested_folds(dat$exposure_binary,
    k_outer = 3, k_inner = 3
  )
  res <- quiet(xgboost_nested_cv(dat$beta, dat$exposure_binary, bin_folds,
    hyperparam_grid = small_grid, nrounds_max = 20, early_stopping_rounds = 5,
    n_cores = 1, verbose = FALSE
  ))

  expect_true(all(c("Accuracy_Mean", "AUC_Mean") %in% names(res$performance)))
  expect_equal(res$consensus_params$objective, "binary:logistic")
})

test_that("xgboost_nested_cv errors when samples do not match", {
  expect_error(
    xgboost_nested_cv(dat$beta, dat$exposure[-1], folds),
    "doesn't match"
  )
})
