# Step 3: elastic net with nested CV.
# Moved from Aim1_Code/methylation_pipeline_v4_shared.R.

#' Elastic Net regression/classification with nested cross-validation
#'
#' Performs elastic net with proper nested CV to avoid overfitting during
#' hyperparameter selection. Inner CV tunes alpha and lambda; outer CV
#' provides unbiased performance estimates. Supports optional covariates
#' that are included but not penalized.
#'
#' @param x CpG beta matrix (CpGs in rows, samples in columns)
#' @param y Response variable (factor for classification, numeric for
#' regression)
#' @param covariates Optional data frame of covariates (samples in rows)
#' @param nested_folds Nested fold structure from create_nested_folds()
#' @param alpha_range Vector of alpha values to search (0=ridge, 1=lasso)
#' @return List with performance metrics, selected CpGs, coefficients, and
#' final model
#' @export
elastic_net_nested_cv <- function(x, y,
                                  covariates = NULL,
                                  nested_folds,
                                  alpha_range = seq(0, 1, 0.1)) {
  # Validate input dimensions
  if (ncol(x) != length(y)) {
    stop(
      "Number of samples in x (", ncol(x),
      ") doesn't match length of y (", length(y), ")"
    )
  }

  if (!is.null(covariates)) {
    if (nrow(covariates) != ncol(x)) {
      stop("Number of rows in covariates must match number of samples in x")
    }
  }

  n_outer <- length(nested_folds)
  outer_results <- list()
  all_selected_cpgs <- list()

  # Combine CpG data with covariates if provided
  # Covariates are appended to feature matrix but will not be penalized
  if (!is.null(covariates)) {
    message("Including unpenalized covariates in elastic net model")
    x_combined <- rbind(x, t(covariates))

    # Penalty factors: 1 = penalized (CpGs), 0 = unpenalized (covariates)
    penalty_factors <- c(rep(1, nrow(x)), rep(0, ncol(covariates)))
    covariate_names <- colnames(covariates)
  } else {
    x_combined <- x
    penalty_factors <- rep(1, nrow(x))
    covariate_names <- NULL
  }

  # Determine evaluation metric based on outcome type
  metric <- ifelse(is.factor(y), "Accuracy", "RMSE")

  # Create progress bar for outer folds
  pb_outer <- progress_bar$new(
    format = paste0(
      "  Elastic Net outer folds [:bar] :percent | ",
      "Fold :current/:total | ETA: :eta"
    ),
    total = n_outer,
    clear = FALSE,
    width = 80
  )

  for (i in 1:n_outer) {
    message(paste("\nOuter fold", i, "of", n_outer))

    # Split data into outer training and test sets
    train_idx <- nested_folds[[i]]$outer_train
    test_idx <- nested_folds[[i]]$outer_test
    inner_folds <- nested_folds[[i]]$inner_folds

    x_train <- x_combined[, train_idx, drop = FALSE]
    x_test <- x_combined[, test_idx, drop = FALSE]
    y_train <- y[train_idx]
    y_test <- y[test_idx]

    # ----------------------------------------
    # INNER CV: Grid search over alpha values
    # Lambda is selected automatically by cv.glmnet
    # ----------------------------------------
    best_alpha <- NA
    best_lambda <- NA
    # Lower is better for both deviance and RMSE
    best_score <- ifelse(metric == "RMSE", Inf, Inf)

    # Progress bar for alpha grid search
    pb_alpha <- progress_bar$new(
      format = "    Tuning alpha [:bar] :percent | Alpha :current/:total",
      total = length(alpha_range),
      clear = FALSE,
      width = 70
    )

    for (alpha in alpha_range) {
      # Run cv.glmnet with inner folds for lambda selection
      # parallel = TRUE enables parallel computation of the CV folds
      cv_model <- cv.glmnet(
        x = t(x_train),
        y = y_train,
        alpha = alpha,
        penalty.factor = penalty_factors,
        nfolds = length(inner_folds),
        foldid = create_foldid(inner_folds, train_idx),
        family = ifelse(is.factor(y_train), "binomial", "gaussian"),
        parallel = TRUE # Enable parallel processing for inner CV
      )

      # Get cross-validated performance at lambda.min
      # cvm is cross-validated mean deviance (lower is better)
      score <- min(cv_model$cvm)
      if (score < best_score) {
        best_score <- score
        best_alpha <- alpha
        best_lambda <- cv_model$lambda.min
      }

      pb_alpha$tick()
    }

    # ----------------------------------------
    # Refit model with optimal hyperparameters on full outer training set
    # ----------------------------------------
    final_model <- glmnet(
      x = t(x_train),
      y = y_train,
      alpha = best_alpha,
      lambda = best_lambda,
      penalty.factor = penalty_factors,
      family = ifelse(is.factor(y_train), "binomial", "gaussian")
    )

    # Extract coefficients (excluding intercept)
    coefs <- as.matrix(coef(final_model, s = best_lambda))

    # Separate CpG coefficients from covariate coefficients
    # Row 1 is intercept, rows 2:(nrow(x)+1) are CpGs
    cpg_coefs <- coefs[2:(nrow(x) + 1), 1]
    selected_features <- rownames(x)[which(abs(cpg_coefs) > 0)]
    all_selected_cpgs[[i]] <- selected_features

    # Extract covariate coefficients if present
    if (!is.null(covariates)) {
      covariate_coefs <- coefs[(nrow(x) + 2):nrow(coefs), 1]
      names(covariate_coefs) <- covariate_names
    } else {
      covariate_coefs <- NULL
    }

    # ----------------------------------------
    # Evaluate on outer test set (unbiased estimate)
    # ----------------------------------------
    predictions <- predict(final_model,
      newx = t(x_test),
      s = best_lambda, type = "response"
    )

    # Calculate performance metrics based on outcome type
    if (is.factor(y)) {
      # Classification metrics
      pred_classes <- ifelse(predictions > 0.5, levels(y)[2], levels(y)[1])
      pred_classes <- factor(pred_classes, levels = levels(y))

      accuracy <- mean(pred_classes == y_test)

      # Calculate precision, recall, and F1 score
      confusion <- table(Predicted = pred_classes, Actual = y_test)
      if (nrow(confusion) == 2 && ncol(confusion) == 2) {
        precision <- confusion[2, 2] / sum(confusion[2, ])
        recall <- confusion[2, 2] / sum(confusion[, 2])
        f1_score <- ifelse(is.na(precision) || is.na(recall) ||
          (precision + recall) == 0, 0,
        2 * (precision * recall) / (precision + recall)
        )
      } else {
        f1_score <- NA
      }

      # Calculate AUC for binary classification
      if (nlevels(y) == 2) {
        auc_value <- auc(roc(y_test, predictions[, 1], quiet = TRUE))
      } else {
        auc_value <- NA
      }

      outer_results[[i]] <- list(
        accuracy = accuracy,
        f1_score = f1_score,
        auc = auc_value,
        predictions = pred_classes,
        actual = y_test,
        n_features = length(selected_features),
        best_alpha = best_alpha,
        best_lambda = best_lambda,
        covariate_coefficients = covariate_coefs
      )
    } else {
      # Regression metrics
      rmse <- sqrt(mean((predictions - y_test)^2))
      mae <- mean(abs(predictions - y_test))
      r_squared <-
        1 - sum((y_test - predictions)^2) / sum((y_test - mean(y_test))^2)

      outer_results[[i]] <- list(
        rmse = rmse,
        mae = mae,
        r_squared = r_squared,
        predictions = as.vector(predictions),
        actual = y_test,
        n_features = length(selected_features),
        best_alpha = best_alpha,
        best_lambda = best_lambda,
        covariate_coefficients = covariate_coefs
      )
    }

    pb_outer$tick()
  }

  # ============================================
  # AGGREGATE RESULTS ACROSS OUTER FOLDS
  # ============================================
  if (is.factor(y)) {
    avg_accuracy <- mean(sapply(outer_results, function(x) x$accuracy))
    avg_f1 <- mean(sapply(outer_results, function(x) x$f1_score), na.rm = TRUE)
    avg_auc <- mean(sapply(outer_results, function(x) x$auc), na.rm = TRUE)

    final_metrics <- data.frame(
      Model = "Elastic Net",
      Accuracy_Mean = avg_accuracy,
      Accuracy_SD = sd(sapply(outer_results, function(x) x$accuracy)),
      F1_Mean = avg_f1,
      F1_SD = sd(sapply(outer_results, function(x) x$f1_score), na.rm = TRUE),
      AUC_Mean = avg_auc,
      AUC_SD = sd(sapply(outer_results, function(x) x$auc), na.rm = TRUE),
      Avg_Features = mean(sapply(outer_results, function(x) x$n_features))
    )
  } else {
    avg_rmse <- mean(sapply(outer_results, function(x) x$rmse))
    avg_r2 <- mean(sapply(outer_results, function(x) x$r_squared))

    final_metrics <- data.frame(
      Model = "Elastic Net",
      RMSE_Mean = avg_rmse,
      RMSE_SD = sd(sapply(outer_results, function(x) x$rmse)),
      R2_Mean = avg_r2,
      R2_SD = sd(sapply(outer_results, function(x) x$r_squared)),
      Avg_Features = mean(sapply(outer_results, function(x) x$n_features))
    )
  }

  # Average covariate coefficients across folds (for interpretation)
  if (!is.null(covariates)) {
    cov_mat <- sapply(
      outer_results, function(z) as.numeric(z$covariate_coefficients)
    )
    avg_covariate_coefs <- rowMeans(cov_mat, na.rm = TRUE)
    names(avg_covariate_coefs) <- covariate_names

    message("\nAverage covariate coefficients across folds:")
    print(avg_covariate_coefs)
  } else {
    avg_covariate_coefs <- NULL
  }

  # Calculate CpG selection frequency across outer folds
  # CpGs selected in more folds are more stable/reliable
  all_cpgs <- unique(unlist(all_selected_cpgs))
  cpg_frequency <- sapply(all_cpgs, function(cpg) {
    sum(sapply(all_selected_cpgs, function(fold_cpgs) cpg %in% fold_cpgs))
  })

  # ============================================
  # TRAIN FINAL MODEL ON ALL DATA
  # Use average hyperparameters from outer folds
  # ============================================
  message("\nTraining final model on all data...")
  final_model <- glmnet(
    x = t(x_combined),
    y = y,
    alpha = mean(sapply(outer_results, function(x) x$best_alpha)),
    lambda = mean(sapply(outer_results, function(x) x$best_lambda)),
    penalty.factor = penalty_factors,
    family = ifelse(is.factor(y), "binomial", "gaussian")
  )

  final_coefs <- coef(final_model)
  final_coefs_matrix <- as.matrix(final_coefs)

  list(
    performance = final_metrics,
    outer_fold_results = outer_results,
    selected_cpgs = names(sort(cpg_frequency, decreasing = TRUE)),
    cpg_selection_frequency = sort(cpg_frequency, decreasing = TRUE),
    # Stable means selected in 80%+ of folds
    stable_cpgs = names(cpg_frequency[cpg_frequency >= n_outer * 0.8]),
    final_model = final_model,
    final_coefs_matrix = final_coefs_matrix,
    # CpG coefficients only (excluding intercept)
    coefficients = final_coefs_matrix[2:(nrow(x) + 1), 1],
    covariate_coefficients = avg_covariate_coefs,
    penalty_factors = penalty_factors
  )
}
