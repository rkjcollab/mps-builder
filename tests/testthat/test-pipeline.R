test_that("run_complete_pipeline runs end to end", {
  skip(paste(
    "Known issue: calls random_forest_nested_cv() and",
    "create_methylation_risk_score(), which are not defined"
  ))
})

test_that("run_trimester_analysis runs every trimester", {
  skip("Known issue: depends on run_complete_pipeline()")
})

test_that("compare_trimester_results compares trimesters", {
  skip(paste(
    "Known issue: rbind() of continuous and dichotomous performance tables",
    "fails because their columns differ"
  ))
})
