# Step 5: create methylation profile score.
# Moved from Aim1_Code/methylation_pipeline_v4_shared.R.

#' Create methylation profile score (MPS) from model results
#'
#' Calculates a weighted or unweighted sum of methylation values at selected
#' CpG sites to create a summary score representing PM2.5 exposure signature.
#'
#' @param beta_matrix Full beta matrix (CpGs in rows, samples in columns)
#' @param model_results Output from elastic_net_nested_cv() or
#' random_forest_nested_cv()
#' @param method "elastic_net" or "random_forest"
#' @param score_type "weighted" (by coefficients/importance) or "unweighted"
#' (mean)
#' @return Numeric vector of standardized profile scores (one per sample)
#' @export
create_mps <- function(
  beta_matrix, model_results, method = c("elastic_net", "random_forest"),
  score_type = c("weighted", "unweighted")
) {
  method <- match.arg(method)
  score_type <- match.arg(score_type)

  if (method == "elastic_net") {
    # Extract coefficients from elastic net model
    coefs <- model_results$coefficients
    if (is.matrix(coefs)) {
      coefs <- as.vector(coefs)
      names(coefs) <- rownames(model_results$coefficients)
    }

    # Get CpGs with non-zero coefficients (selected by elastic net)
    selected_cpgs <- names(coefs[abs(coefs) > 0])

    if (length(selected_cpgs) == 0) {
      warning("No CpGs selected by elastic net")
      return(rep(0, ncol(beta_matrix)))
    }

    # Subset beta matrix to selected CpGs
    beta_subset <- beta_matrix[selected_cpgs, , drop = FALSE]

    if (score_type == "weighted") {
      # Weighted sum using elastic net coefficients
      weights <- coefs[selected_cpgs]
      risk_score <- colSums(beta_subset * weights, na.rm = TRUE)
    } else {
      # Simple unweighted mean across selected CpGs
      risk_score <- colMeans(beta_subset, na.rm = TRUE)
    }
  } else if (method == "random_forest") {
    # Use top important CpGs from random forest
    selected_cpgs <- model_results$top_cpgs

    # Ensure selected CpGs exist in beta matrix
    selected_cpgs <- intersect(selected_cpgs, rownames(beta_matrix))

    if (length(selected_cpgs) == 0) {
      warning("No selected CpGs found in beta matrix")
      return(rep(0, ncol(beta_matrix)))
    }

    beta_subset <- beta_matrix[selected_cpgs, , drop = FALSE]

    if (score_type == "weighted") {
      # Weight by variable importance scores
      weights <- model_results$avg_importance[selected_cpgs]
      w_sum <- sum(weights)
      if (w_sum == 0) {
        warning(
          "Sum of importance weights is zero; using unweighted mean instead."
        )
        risk_score <- colMeans(beta_subset, na.rm = TRUE)
      } else {
        weights <- weights / w_sum # Normalize weights to sum to 1
        risk_score <- colSums(beta_subset * weights, na.rm = TRUE)
      }
    } else {
      risk_score <- colMeans(beta_subset, na.rm = TRUE)
    }
  }

  # Z-score standardize risk score to z-scores (mean=0, sd=1)
  risk_score_z <- scale(risk_score)[, 1]

  data.frame(
    sample_id        = names(risk_score),
    risk_score       = risk_score,
    risk_score_z     = risk_score_z,
    row.names        = NULL,
    stringsAsFactors = FALSE
  )
}
