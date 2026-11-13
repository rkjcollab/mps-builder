# Visualization and analysis functions.
# Moved from Aim1_Code/methylation_pipeline_v4_shared.R.

#' Visualize model performance comparison
#'
#' @param results Output from run_complete_pipeline()
#' @export
plot_model_comparison <- function(results) {
  perf_data <- results$comparison$comparison_table

  # Plot R² for continuous outcomes
  if ("R2_Mean" %in% names(perf_data)) {
    cont_data <- perf_data[perf_data$Outcome == "Continuous", ]

    p1 <- ggplot(cont_data, aes(x = Model, y = R2_Mean)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_errorbar(
        aes(
          ymin = R2_Mean - R2_SD,
          ymax = R2_Mean + R2_SD
        ),
        width = 0.2
      ) +
      theme_minimal() +
      labs(
        title = "Model Performance - Continuous PM2.5",
        y = "R\u00b2 (Mean \u00b1 SD)",
        x = "Model"
      ) +
      ylim(0, 1)

    print(p1)
  }

  # Plot AUC for dichotomous outcomes
  if ("AUC_Mean" %in% names(perf_data)) {
    dicho_data <- perf_data[perf_data$Outcome == "Dichotomous", ]

    p2 <- ggplot(dicho_data, aes(x = Model, y = AUC_Mean)) +
      geom_bar(stat = "identity", fill = "darkred") +
      geom_errorbar(
        aes(
          ymin = AUC_Mean - AUC_SD,
          ymax = AUC_Mean + AUC_SD
        ),
        width = 0.2
      ) +
      theme_minimal() +
      labs(
        title = "Model Performance - Dichotomous PM2.5",
        y = "AUC (Mean \u00b1 SD)",
        x = "Model"
      ) +
      ylim(0, 1)

    print(p2)
  }
}

#' Export all results to files
#'
#' @param results Output from run_complete_pipeline()
#' @param beta_matrix Original beta matrix
#' @param output_dir Directory for output files
#' @export
export_results <- function(
  results,
  beta_matrix,
  output_dir = "methylation_results"
) {
  dir.create(output_dir, showWarnings = FALSE)

  # Export elastic net CpGs with selection frequency and coefficients
  en_cpgs <- data.frame(
    CpG = results$continuous$elastic_net$selected_cpgs,
    Selection_Frequency =
      results$continuous$elastic_net$cpg_selection_frequency[
        results$continuous$elastic_net$selected_cpgs
      ],
    Coefficient = results$continuous$elastic_net$coefficients[
      results$continuous$elastic_net$selected_cpgs
    ]
  )

  write.csv(en_cpgs,
    file.path(output_dir, "elastic_net_cpgs.csv"),
    row.names = FALSE
  )

  # Export random forest CpGs with importance scores
  rf_cpgs <- data.frame(
    CpG = names(results$continuous$random_forest$avg_importance),
    Importance = results$continuous$random_forest$avg_importance
  )

  write.csv(rf_cpgs,
    file.path(output_dir, "random_forest_cpgs.csv"),
    row.names = FALSE
  )

  # Export model performance comparison
  write.csv(results$comparison$comparison_table,
    file.path(output_dir, "model_performance.csv"),
    row.names = FALSE
  )

  # Export methylation profile scores for each sample
  profile_scores <- data.frame(
    Sample_ID = colnames(beta_matrix),
    ElasticNet_Score = results$continuous$elastic_net_score,
    RandomForest_Score = results$continuous$random_forest_score
  )

  write.csv(profile_scores,
    file.path(output_dir, "profile_scores.csv"),
    row.names = FALSE
  )

  message(paste("Results exported to", output_dir))
}
