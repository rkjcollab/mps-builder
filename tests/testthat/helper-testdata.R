# Load the simulated data files in tests/testthat/fixtures/ (made by
# make_testdata.R in that folder) into the objects the package functions
# expect.
load_test_data <- function() {
  beta <- as.matrix(read.csv(testthat::test_path("fixtures", "beta.csv"),
    row.names = 1, check.names = FALSE
  ))
  pheno <- read.csv(testthat::test_path("fixtures", "pheno.csv"))
  literature <- read.csv(
    testthat::test_path("fixtures", "literature_cpgs.csv")
  )

  covariates <- pheno[, c("maternal_age", "infant_sex")]
  rownames(covariates) <- pheno$sample_id

  list(
    beta = beta,
    exposure = setNames(pheno$exposure, pheno$sample_id),
    exposure_binary = setNames(factor(pheno$exposure_binary), pheno$sample_id),
    covariates = covariates,
    # same construction as the original analysis script
    literature_studies = lapply(split(literature$cpg, literature$study), unique)
  )
}

# Run an expression without its messages and printed output.
quiet <- function(expr) {
  result <- NULL
  utils::capture.output(result <- suppressMessages(expr))
  result
}

# cv.glmnet(parallel = TRUE) and foreach() expect a registered backend; the
# original script registered doParallel. Tests run sequentially.
foreach::registerDoSEQ()
