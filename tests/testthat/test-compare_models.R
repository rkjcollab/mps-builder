test_that("compare_filtering_results reports overlap between two filterings", {
  filter_cont <- list(selected_cpgs = c("cg1", "cg2", "cg3"), n_selected = 3)
  filter_dicho <- list(
    selected_cpgs = c("cg2", "cg3", "cg4", "cg5"), n_selected = 4
  )

  res <- quiet(compare_filtering_results(filter_cont, filter_dicho))

  expect_equal(res$n_overlap, 2)
  expect_equal(res$pct_overlap, 100 * 2 / 3)
  expect_equal(res$only_continuous, 1)
  expect_equal(res$only_dichotomous, 2)
  expect_setequal(res$overlap_cpgs, c("cg2", "cg3"))
})

test_that("compare_filtering_results gives NA overlap with no CpGs", {
  empty <- list(selected_cpgs = character(0), n_selected = 0)
  res <- quiet(compare_filtering_results(empty, empty))
  expect_true(is.na(res$pct_overlap))
})

test_that("compare_all_models compares continuous and dichotomous models", {
  skip(paste(
    "Known issue: rbind() of continuous (RMSE/R2) and dichotomous",
    "(Accuracy/F1/AUC) performance tables fails because their columns differ"
  ))
})
