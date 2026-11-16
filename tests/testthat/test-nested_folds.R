dat <- load_test_data()

test_that("create_nested_folds partitions samples into outer and inner folds", {
  y <- dat$exposure
  folds <- create_nested_folds(y, k_outer = 4, k_inner = 3)

  expect_length(folds, 4)

  # outer test folds cover every sample exactly once
  all_test <- unlist(lapply(folds, `[[`, "outer_test"))
  expect_setequal(all_test, seq_along(y))
  expect_equal(length(all_test), length(y))

  for (f in folds) {
    expect_setequal(c(f$outer_test, f$outer_train), seq_along(y))
    expect_length(intersect(f$outer_test, f$outer_train), 0)

    # inner folds partition the outer training set, in original indices
    expect_length(f$inner_folds, 3)
    inner <- unlist(f$inner_folds)
    expect_setequal(inner, f$outer_train)
    expect_equal(length(inner), length(f$outer_train))
  }
})

test_that("create_nested_folds is reproducible with the same seed", {
  a <- create_nested_folds(dat$exposure, k_outer = 3, k_inner = 3, seed = 10)
  b <- create_nested_folds(dat$exposure, k_outer = 3, k_inner = 3, seed = 10)
  expect_identical(a, b)
})

test_that("create_nested_folds keeps both classes in each outer test fold", {
  folds <- create_nested_folds(dat$exposure_binary, k_outer = 3, k_inner = 3)
  for (f in folds) {
    test_classes <- as.character(dat$exposure_binary[f$outer_test])
    expect_setequal(test_classes, c("High", "Low"))
  }
})

test_that("create_foldid maps inner folds to training-set positions", {
  train_idx <- c(2, 4, 6, 8, 10)
  inner_folds <- list(c(4, 10), c(2, 8), 6)

  expect_equal(create_foldid(inner_folds, train_idx), c(2L, 1L, 3L, 2L, 1L))
})
