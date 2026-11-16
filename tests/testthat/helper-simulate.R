# Simulate a small methylation dataset for tests: CpGs x samples M and beta
# values, a continuous phenotype, and a set of CpGs with a planted effect.
simulate_methylation <- function(n_samples = 30, n_cpgs = 200, n_signal = 20,
                                 seed = 1) {
  set.seed(seed)
  pheno <- rnorm(n_samples)

  # baseline M values per CpG, plus noise
  m_values <- matrix(
    rnorm(n_cpgs * n_samples, mean = rnorm(n_cpgs, 0, 1.5), sd = 0.5),
    nrow = n_cpgs,
    ncol = n_samples,
    dimnames = list(
      sprintf("cg%08d", seq_len(n_cpgs)),
      sprintf("S%02d", seq_len(n_samples))
    )
  )

  # plant signal in the first n_signal CpGs
  signal_rows <- seq_len(n_signal)
  effects <- rnorm(n_signal, 0, 0.4)
  m_values[signal_rows, ] <- m_values[signal_rows, ] + outer(effects, pheno)

  beta <- 2^m_values / (2^m_values + 1)

  list(
    M = m_values,
    beta = beta,
    pheno = setNames(pheno, colnames(m_values)),
    true_cpgs = rownames(m_values)[signal_rows]
  )
}
