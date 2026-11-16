# Generate the small simulated data files used by the unit tests.
# Run from the repository root: Rscript tests/testthat/fixtures/make_testdata.R
#
# Files are written to tests/testthat/fixtures/ and mirror the inputs of the
# original analysis script: a beta matrix (CpGs x samples), a sample table with
# exposure and covariates, and a long table of literature CpGs by study.

source("tests/testthat/helper-simulate.R")

out_dir <- "tests/testthat/fixtures"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sim <- simulate_methylation(
  n_samples = 40, n_cpgs = 100, n_signal = 10, seed = 1
)

set.seed(2)
samples <- colnames(sim$beta)
pheno <- data.frame(
  sample_id = samples,
  exposure = round(sim$pheno, 4),
  # median split, as in the original script's example usage
  exposure_binary = ifelse(
    sim$pheno >= stats::median(sim$pheno), "High", "Low"
  ),
  maternal_age = round(rnorm(length(samples), 30, 5), 1),
  infant_sex = rbinom(length(samples), 1, 0.5)
)

# Three studies: two share some signal CpGs, one has only non-signal CpGs
cpgs <- rownames(sim$beta)
literature <- data.frame(
  study = rep(c("StudyA", "StudyB", "StudyC"), times = c(6, 5, 3)),
  cpg = c(
    cpgs[c(1:4, 50, 60)],
    cpgs[c(3:5, 70, 80)],
    cpgs[c(90, 95, 99)]
  )
)

beta <- data.frame(cpg = cpgs, round(sim$beta, 4), check.names = FALSE)

utils::write.csv(beta, file.path(out_dir, "beta.csv"), row.names = FALSE)
utils::write.csv(pheno, file.path(out_dir, "pheno.csv"), row.names = FALSE)
utils::write.csv(literature, file.path(out_dir, "literature_cpgs.csv"),
  row.names = FALSE
)
