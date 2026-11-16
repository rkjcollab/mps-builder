# Visualization and analysis functions.
# Moved from Aim1_Code/methylation_pipeline_v4_shared.R.

#' Create Venn diagram showing CpG filtering criteria overlap
#'
#' @param filter_results Output from filter_cpg_sites()
#' @export
plot_filtering_venn <- function(filter_results) {
  # Extract CpGs meeting each criterion
  variance_cpgs <- rownames(filter_results$criteria_matrix)[
    filter_results$criteria_matrix[, "Variance"] == 1
  ]
  association_cpgs <- rownames(filter_results$criteria_matrix)[
    filter_results$criteria_matrix[, "Association"] == 1
  ]
  literature_cpgs <- rownames(filter_results$criteria_matrix)[
    filter_results$criteria_matrix[, "Literature"] == 1
  ]

  # Generate three-way Venn diagram
  venn <- venn.diagram(
    x = list(
      Variance = variance_cpgs,
      Association = association_cpgs,
      Literature = literature_cpgs
    ),
    filename = NULL,
    fill = c("cornflowerblue", "green", "yellow"),
    alpha = 0.5,
    cex = 1.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.cex = 1.5,
    cat.fontfamily = "serif",
    margin = 0.1
  )

  grid.newpage()
  grid.draw(venn)
}

#' Assess sensitivity of filtering to parameter choices
#'
#' Runs filtering with different combinations of thresholds to evaluate
#' robustness of CpG selection.
#'
#' @param beta_matrix Numeric matrix of beta values
#' @param pm25_data PM2.5 exposure data
#' @param literature_studies Optional literature CpG lists
#' @param var_percentiles Vector of variance percentiles to test
#' @param pval_thresholds Vector of p-value thresholds to test
#' @param min_criteria_options Vector of min_criteria values to test
#' @return Data frame with results for each parameter combination
#' @export
sensitivity_analysis_filtering <- function(
  beta_matrix, pm25_data,
  literature_studies,
  var_percentiles = c(70, 75, 80, 85, 90),
  pval_thresholds = c(0.01, 0.05, 0.1),
  min_criteria_options = c(1, 2, 3)
) {
  # Create grid of all parameter combinations
  results_grid <- expand.grid(
    var_percentile = var_percentiles,
    pval_threshold = pval_thresholds,
    min_criteria = min_criteria_options,
    n_selected = NA,
    pct_literature = NA
  )

  # Progress bar for sensitivity analysis
  pb <- progress_bar$new(
    format = paste0(
      "  Sensitivity analysis [:bar] :percent | ",
      ":current/:total combinations"
    ),
    total = nrow(results_grid),
    clear = FALSE,
    width = 80
  )

  for (i in seq_len(nrow(results_grid))) {
    filter_res <- filter_cpg_sites(
      beta_matrix = beta_matrix,
      pm25_data = pm25_data,
      literature_studies = literature_studies,
      var_percentile = results_grid$var_percentile[i],
      pval_threshold = results_grid$pval_threshold[i],
      min_criteria = results_grid$min_criteria[i]
    )

    results_grid$n_selected[i] <- filter_res$n_selected

    # Calculate percentage of selected CpGs that are from literature
    if (!is.null(filter_res$literature_info)) {
      lit_in_selected <- sum(filter_res$selected_cpgs %in%
        filter_res$literature_info$literature_cpgs)
      if (filter_res$n_selected > 0) {
        results_grid$pct_literature[i] <-
          100 * lit_in_selected / filter_res$n_selected
      } else {
        results_grid$pct_literature[i] <- NA_real_
      }
    }

    pb$tick()
  }

  results_grid
}

#' Compare CpG selection across different PM2.5 exposure windows
#'
#' @param beta_matrix Numeric matrix of beta values
#' @param pm25_data_list Named list of PM2.5 data for different periods
#' @param literature_studies Optional literature CpG lists
#' @param var_percentile Variance filter threshold
#' @param pval_threshold P-value threshold
#' @param min_criteria Minimum criteria for selection
#' @return List with per-period results, consistent CpGs, and overlap matrix
#' @export
compare_temporal_filtering <- function(beta_matrix, pm25_data_list,
                                       literature_studies,
                                       var_percentile = 80,
                                       pval_threshold = 0.05,
                                       min_criteria = 2) {
  temporal_results <- list()

  # Run filtering for each exposure window
  message("\nFiltering CpGs for each exposure window...")
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent | :current/:total periods",
    total = length(pm25_data_list),
    clear = FALSE,
    width = 70
  )

  for (period in names(pm25_data_list)) {
    temporal_results[[period]] <- filter_cpg_sites(
      beta_matrix = beta_matrix,
      pm25_data = pm25_data_list[[period]],
      literature_studies = literature_studies,
      var_percentile = var_percentile,
      pval_threshold = pval_threshold,
      min_criteria = min_criteria
    )
    pb$tick()
  }

  # Find CpGs selected in ALL time periods
  all_cpg_lists <- lapply(temporal_results, function(x) x$selected_cpgs)
  consistent_cpgs <- Reduce(intersect, all_cpg_lists)

  # Create pairwise overlap matrix
  periods <- names(temporal_results)
  overlap_matrix <- matrix(0,
    nrow = length(periods), ncol = length(periods),
    dimnames = list(periods, periods)
  )

  for (i in seq_along(periods)) {
    for (j in seq_along(periods)) {
      overlap_matrix[i, j] <- length(intersect(
        temporal_results[[periods[i]]]$selected_cpgs,
        temporal_results[[periods[j]]]$selected_cpgs
      ))
    }
  }

  # Create summary statistics
  summary_df <- data.frame(
    Period = periods,
    N_Selected = sapply(temporal_results, function(x) x$n_selected),
    N_Variance = sapply(temporal_results, function(x) {
      sum(x$criteria_matrix[, "Variance"])
    }),
    N_Association = sapply(temporal_results, function(x) {
      sum(x$criteria_matrix[, "Association"])
    }),
    N_Literature = sapply(temporal_results, function(x) {
      sum(x$criteria_matrix[, "Literature"])
    })
  )

  list(
    by_period = temporal_results,
    consistent_cpgs = consistent_cpgs,
    n_consistent = length(consistent_cpgs),
    overlap_matrix = overlap_matrix,
    summary = summary_df
  )
}

#' Detailed analysis of literature CpGs in filtering results
#'
#' @param filter_results Output from filter_cpg_sites()
#' @param beta_matrix Numeric matrix of beta values
#' @param pm25_data PM2.5 exposure data
#' @return List with per-CpG details and summary statistics
#' @export
analyze_literature_cpgs <- function(filter_results, beta_matrix, pm25_data) {
  if (is.null(filter_results$literature_info)) {
    message("No literature studies provided")
    return(NULL)
  }

  lit_info <- filter_results$literature_info

  # Group CpGs by number of studies they appear in (replication count)
  cpgs_by_n_studies <- list()
  for (n in 1:4) {
    cpgs_by_n_studies[[paste0("in_", n, "_studies")]] <-
      names(lit_info$n_studies[lit_info$n_studies == n])
  }

  # Get literature CpGs that also met other criteria
  lit_cpgs <- rownames(filter_results$criteria_matrix)[
    filter_results$criteria_matrix[, "Literature"] == 1
  ]

  lit_criteria_summary <- data.frame(
    Literature_CpG = lit_cpgs,
    Passes_Variance = filter_results$criteria_matrix[lit_cpgs, "Variance"],
    Passes_Association =
      filter_results$criteria_matrix[lit_cpgs, "Association"],
    N_Studies = lit_info$n_studies[lit_cpgs],
    Selected = as.integer(lit_cpgs %in% filter_results$selected_cpgs)
  )

  # Calculate effect sizes (correlation) for literature CpGs
  if (is.numeric(pm25_data)) {
    lit_criteria_summary$Correlation <- sapply(lit_cpgs, function(cpg) {
      cor(beta_matrix[cpg, ], pm25_data, use = "complete.obs")
    })
  }

  # Generate summary statistics
  summary_stats <- list(
    total_literature_cpgs = length(lit_cpgs),
    lit_cpgs_selected = sum(lit_criteria_summary$Selected),
    lit_cpgs_also_variance = sum(lit_criteria_summary$Passes_Variance),
    lit_cpgs_also_association = sum(lit_criteria_summary$Passes_Association),
    lit_cpgs_all_three = sum(lit_criteria_summary$Passes_Variance &
      lit_criteria_summary$Passes_Association)
  )

  list(
    cpgs_by_study_count = cpgs_by_n_studies,
    literature_cpg_details = lit_criteria_summary,
    summary = summary_stats
  )
}
