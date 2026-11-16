dat <- load_test_data()
folds <- create_nested_folds(dat$exposure, k_outer = 3, k_inner = 3)
# The default grid uses mtry of 250 to 1000, more than the 100 test CpGs
small_grid <- expand.grid(
  mtry = c(5, 20), min.node.size = 5, sample.fraction = 0.8
)

test_that("random_forest_nested_cv_tuned runs for a continuous exposure", {
  res <- quiet(random_forest_nested_cv_tuned(dat$beta, dat$exposure, folds,
    ntree = 50, hyperparam_grid = small_grid, n_cores = 1
  ))

  expect_equal(
    names(res$performance),
    c("Model", "RMSE_Mean", "RMSE_SD", "R2_Mean", "R2_SD", "Avg_Features")
  )
  expect_length(res$outer_fold_results, 3)
  expect_length(res$all_fold_hyperparameters, 3)
  expect_true(res$best_hyperparameters$mtry %in% c(5, 20))

  # importance is averaged over folds and sorted; top CpGs are capped at 100
  expect_setequal(names(res$avg_importance), rownames(dat$beta))
  expect_false(is.unsorted(rev(res$avg_importance)))
  expect_equal(res$top_cpgs, names(res$avg_importance)[1:100])
  expect_s3_class(res$final_model, "ranger")
})

test_that("random_forest_nested_cv_tuned includes covariates as features", {
  res <- quiet(random_forest_nested_cv_tuned(dat$beta, dat$exposure, folds,
    covariates = dat$covariates, ntree = 50, hyperparam_grid = small_grid,
    n_cores = 1
  ))
  expect_true(all(names(dat$covariates) %in% names(res$avg_importance)))
  expect_equal(
    res$performance$Avg_Features, nrow(dat$beta) + ncol(dat$covariates)
  )
})

test_that("random_forest_nested_cv_tuned runs for a dichotomous exposure", {
  skip(paste(
    "Known issue: checks ranger_model$probability, which ranger objects",
    "do not have (the type is in $treetype), so the if() errors"
  ))
  bin_folds <- create_nested_folds(dat$exposure_binary,
    k_outer = 3, k_inner = 3
  )
  res <- quiet(random_forest_nested_cv_tuned(
    dat$beta, dat$exposure_binary, bin_folds,
    ntree = 50, hyperparam_grid = small_grid, n_cores = 1
  ))

  perf_names <- names(res$performance)
  expect_true(all(c("Accuracy_Mean", "F1_Mean", "AUC_Mean") %in% perf_names))
  for (fold in res$outer_fold_results) {
    expect_s3_class(fold$predictions, "factor")
  }
})
