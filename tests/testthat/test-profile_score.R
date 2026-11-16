dat <- load_test_data()

# Minimal model results with known weights, so scores can be checked by hand
cpgs <- rownames(dat$beta)
en_results <- list(
  coefficients = setNames(c(0.5, -1, rep(0, length(cpgs) - 2)), cpgs)
)
rf_results <- list(
  top_cpgs = cpgs[1:3],
  avg_importance = setNames(c(3, 1, 0, rep(0, length(cpgs) - 3)), cpgs)
)

test_that("weighted elastic net score is a weighted sum, z-scored", {
  res <- create_mps(dat$beta, en_results,
    method = "elastic_net"
  )

  expected <- 0.5 * dat$beta[1, ] - 1 * dat$beta[2, ]
  expect_equal(res$sample_id, colnames(dat$beta))
  expect_equal(res$risk_score, unname(expected))
  expect_equal(res$risk_score_z, as.numeric(scale(expected)))
})

test_that("unweighted elastic net score is the mean of selected CpGs", {
  res <- create_mps(dat$beta, en_results,
    method = "elastic_net", score_type = "unweighted"
  )
  expect_equal(res$risk_score, unname(colMeans(dat$beta[1:2, ])))
})

test_that("weighted random forest score uses normalized importance", {
  res <- create_mps(dat$beta, rf_results,
    method = "random_forest"
  )
  expected <- 0.75 * dat$beta[1, ] + 0.25 * dat$beta[2, ]
  expect_equal(res$risk_score, unname(expected))
})

test_that("score works on elastic_net_nested_cv output", {
  folds <- create_nested_folds(dat$exposure, k_outer = 3, k_inner = 3)
  en <- quiet(elastic_net_nested_cv(dat$beta, dat$exposure,
    nested_folds = folds, alpha_range = 1
  ))
  res <- create_mps(dat$beta, en, method = "elastic_net")
  expect_equal(nrow(res), ncol(dat$beta))
  expect_equal(mean(res$risk_score_z), 0, tolerance = 1e-8)
})

test_that("score warns and returns zeros when no CpGs are selected", {
  none <- list(coefficients = setNames(rep(0, length(cpgs)), cpgs))
  expect_warning(
    res <- create_mps(dat$beta, none,
      method = "elastic_net"
    ),
    "No CpGs selected"
  )
  expect_equal(res, rep(0, ncol(dat$beta)))
})
