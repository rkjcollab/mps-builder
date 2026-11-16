# Step 4: random forest with nested CV.
# Moved from Aim1_Code/methylation_pipeline_v4_shared.R.

#' Random Forest with nested cross-validation
#'
#' Performs random forest with nested CV for unbiased performance estimation.
#' Supports both ranger (faster, permutation importance) and caret RF.
#' Inner CV tunes mtry parameter; outer CV evaluates performance.
#'
#' @param x CpG beta matrix (CpGs in rows, samples in columns)
#' @param y Response variable (factor for classification, numeric for
#' regression)
#' @param nested_folds Nested fold structure from create_nested_folds()
#' @param covariates Optional data frame of covariates (samples in rows)
#' @param ntree Number of trees in the forest (default 500)
#' @param hyperparam_grid Data frame for tuning hyperparameters: mtry,
#' min.node.size, and sample fraction
#' @param importance_threshold Number of top CpGs to retain based on importance
#' @param n_cores Number of cores for parallel computation (default: half of
#'   available cores, as set in the original analysis script)
#' @return List with performance metrics, importance scores, and final model
#' @export
random_forest_nested_cv_tuned <- function(
  x, y, nested_folds,
  covariates = NULL,
  ntree = 500,
  hyperparam_grid = NULL,
  importance_threshold = 100,
  n_cores = max(1, floor(parallel::detectCores() / 2))
) {
  # Default hyperparameter grid if not provided
  if (is.null(hyperparam_grid)) {
    hyperparam_grid <- expand.grid(
      # TODO:Could rewrite this to be adaptive to p/4, p/3, p/2
      mtry = c(250, 500, 750, 1000),
      min.node.size = c(5, 10, 20),
      sample.fraction = c(0.8, 1.0)
    )
  }

  cat(
    "Testing", nrow(hyperparam_grid),
    "hyperparameter combinations in inner CV\n"
  )

  n_outer <- length(nested_folds)
  outer_results <- list()
  all_importance_scores <- list()
  all_best_params <- list()

  # Combine CpG data with covariates if provided
  if (!is.null(covariates)) {
    message("Including covariates in Random Forest model")
    x_combined <- rbind(x, t(covariates))
  } else {
    x_combined <- x
  }

  # Progress bar for outer folds
  pb_outer <- progress_bar$new(
    format = paste0(
      "  RF outer folds [:bar] :percent | ",
      "Fold :current/:total | ETA: :eta"
    ),
    total = n_outer,
    clear = FALSE,
    width = 80
  )

  for (i in 1:n_outer) {
    message(paste("\n=== Outer fold", i, "of", n_outer, "==="))

    # Split data into outer training and test sets
    train_idx <- nested_folds[[i]]$outer_train
    test_idx <- nested_folds[[i]]$outer_test
    inner_folds <- nested_folds[[i]]$inner_folds

    x_train <- x_combined[, train_idx, drop = FALSE]
    x_test <- x_combined[, test_idx, drop = FALSE]
    y_train <- y[train_idx]
    y_test <- y[test_idx]

    # ----------------------------------------
    # INNER CV: Hyperparameter tuning
    # ----------------------------------------
    message("  Tuning hyperparameters on inner CV...")

    best_params <- NULL
    best_score <- Inf

    # Progress bar for hyperparameter search
    pb_hyper <- progress_bar$new(
      format = "    [:bar] :percent | Combo :current/:total",
      total = nrow(hyperparam_grid),
      clear = FALSE,
      width = 60
    )

    for (j in seq_len(nrow(hyperparam_grid))) {
      inner_scores <- numeric(length(inner_folds))

      # Evaluate this hyperparameter combination on inner folds
      for (k in seq_along(inner_folds)) {
        # Split inner training/validation
        inner_train_idx <- setdiff(train_idx, inner_folds[[k]])
        inner_val_idx <- inner_folds[[k]]

        # Train model with these hyperparameters
        model_temp <- ranger(
          x = t(x_combined[, inner_train_idx, drop = FALSE]),
          y = y[inner_train_idx],
          num.trees = ntree,
          mtry = hyperparam_grid$mtry[j],
          min.node.size = hyperparam_grid$min.node.size[j],
          sample.fraction = hyperparam_grid$sample.fraction[j],
          num.threads = 1, # Don't nest parallelization
          classification = is.factor(y),
          verbose = FALSE
        )

        # Predict on inner validation set
        pred_temp <- predict(model_temp,
          data = t(x_combined[, inner_val_idx, drop = FALSE])
        )

        # Calculate validation score (lower is better)
        if (is.factor(y)) {
          inner_scores[k] <-
            1 - mean(pred_temp$predictions == y[inner_val_idx])
        } else {
          inner_scores[k] <-
            sqrt(mean((pred_temp$predictions - y[inner_val_idx])^2))
        }
      }

      # Average score across inner folds
      avg_score <- mean(inner_scores)

      # Track best hyperparameters
      if (avg_score < best_score) {
        best_score <- avg_score
        best_params <- hyperparam_grid[j, ]
      }

      pb_hyper$tick()
    }

    message(paste(
      "\n  Best params: mtry =", best_params$mtry,
      ", min.node.size =", best_params$min.node.size,
      ", sample.fraction =", best_params$sample.fraction
    ))
    message(paste("  Inner CV score:", round(best_score, 4)))

    all_best_params[[i]] <- best_params

    # ----------------------------------------
    # Train final model for this outer fold with best hyperparameters
    # ----------------------------------------
    ranger_model <- ranger(
      x = t(x_train),
      y = y_train,
      num.trees = ntree,
      mtry = best_params$mtry,
      min.node.size = best_params$min.node.size,
      sample.fraction = best_params$sample.fraction,
      importance = "permutation",
      num.threads = n_cores,
      classification = is.factor(y),
      probability = is.factor(y) && nlevels(y) == 2
    )

    # Generate predictions on test set
    pred_ranger <- predict(ranger_model, data = t(x_test))

    if (is.factor(y)) {
      if (ranger_model$probability) {
        prob_mat <- pred_ranger$predictions
        positive_class <- levels(y)[2]
        prob_vec <- prob_mat[, positive_class]
        predictions <- factor(
          ifelse(prob_vec > 0.5, positive_class, levels(y)[1]),
          levels = levels(y)
        )
      } else {
        predictions <- pred_ranger$predictions
        prob_vec <- NULL
      }
    } else {
      predictions <- pred_ranger$predictions
      prob_vec <- NULL
    }

    # Extract importance
    importance_perm <- ranger_model$variable.importance

    # Store importance
    all_importance_scores[[i]] <- importance_perm

    # RF uses all features
    n_features_used <- nrow(x_combined)

    # ----------------------------------------
    # Calculate performance metrics
    # ----------------------------------------
    if (is.factor(y)) {
      accuracy <- mean(predictions == y_test)

      confusion <- table(Predicted = predictions, Actual = y_test)
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

      if (nlevels(y) == 2 && !is.null(prob_vec)) {
        auc_value <- auc(roc(y_test, prob_vec, quiet = TRUE))
      } else {
        auc_value <- NA
      }

      outer_results[[i]] <- list(
        accuracy = accuracy,
        f1_score = f1_score,
        auc = auc_value,
        predictions = predictions,
        actual = y_test,
        n_features_selected = n_features_used,
        best_mtry = best_params$mtry,
        best_min_node_size = best_params$min.node.size,
        best_sample_fraction = best_params$sample.fraction
      )
    } else {
      rmse <- sqrt(mean((predictions - y_test)^2))
      mae <- mean(abs(predictions - y_test))
      r_squared <-
        1 - sum((y_test - predictions)^2) / sum((y_test - mean(y_test))^2)

      outer_results[[i]] <- list(
        rmse = rmse,
        mae = mae,
        r_squared = r_squared,
        predictions = predictions,
        actual = y_test,
        n_features_selected = n_features_used,
        best_mtry = best_params$mtry,
        best_min_node_size = best_params$min.node.size,
        best_sample_fraction = best_params$sample.fraction
      )
    }

    pb_outer$tick()
  }

  # ============================================
  # AGGREGATE RESULTS
  # ============================================
  if (is.factor(y)) {
    avg_accuracy <- mean(sapply(outer_results, function(x) x$accuracy))
    avg_f1 <- mean(sapply(outer_results, function(x) x$f1_score), na.rm = TRUE)
    avg_auc <- mean(sapply(outer_results, function(x) x$auc), na.rm = TRUE)

    final_metrics <- data.frame(
      Model = "Random Forest (Tuned)",
      Accuracy_Mean = avg_accuracy,
      Accuracy_SD = sd(sapply(outer_results, function(x) x$accuracy)),
      F1_Mean = avg_f1,
      F1_SD = sd(sapply(outer_results, function(x) x$f1_score), na.rm = TRUE),
      AUC_Mean = avg_auc,
      AUC_SD = sd(sapply(outer_results, function(x) x$auc), na.rm = TRUE),
      Avg_Features = mean(sapply(
        outer_results, function(x) x$n_features_selected
      ))
    )
  } else {
    avg_rmse <- mean(sapply(outer_results, function(x) x$rmse))
    avg_r2 <- mean(sapply(outer_results, function(x) x$r_squared))

    final_metrics <- data.frame(
      Model = "Random Forest (Tuned)",
      RMSE_Mean = avg_rmse,
      RMSE_SD = sd(sapply(outer_results, function(x) x$rmse)),
      R2_Mean = avg_r2,
      R2_SD = sd(sapply(outer_results, function(x) x$r_squared)),
      Avg_Features = mean(sapply(
        outer_results, function(x) x$n_features_selected
      ))
    )
  }

  # Average importance scores
  all_features <- unique(unlist(lapply(all_importance_scores, names)))
  avg_importance <- sapply(all_features, function(feat) {
    mean(sapply(all_importance_scores, function(imp) {
      if (feat %in% names(imp)) imp[feat] else 0
    }))
  })

  # ============================================
  # SELECT BEST HYPERPARAMETERS (MOST FREQUENT)
  # ============================================

  # Convert list to data frame
  hyperparam_df <- do.call(rbind, lapply(all_best_params, as.data.frame))

  # Function to get mode (most frequent value)
  get_mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }

  best_mtry <- get_mode(hyperparam_df$mtry)
  best_min_node_size <- get_mode(hyperparam_df$min.node.size)
  best_sample_fraction <- get_mode(hyperparam_df$sample.fraction)

  message("\n=== Best Hyperparameters (Most Frequent Across Folds) ===")
  message(paste(
    "mtry:", best_mtry,
    "- chosen in",
    sum(hyperparam_df$mtry == best_mtry),
    "out of", n_outer, "folds"
  ))
  message(paste(
    "min.node.size:", best_min_node_size,
    "- chosen in",
    sum(hyperparam_df$min.node.size == best_min_node_size),
    "out of", n_outer, "folds"
  ))
  message(paste(
    "sample.fraction:", best_sample_fraction,
    "- chosen in",
    sum(hyperparam_df$sample.fraction == best_sample_fraction),
    "out of", n_outer, "folds"
  ))

  # ============================================
  # TRAIN FINAL MODEL ON ALL DATA
  # Using MOST FREQUENT hyperparameters
  # ============================================
  message("\nTraining final model on all data with best hyperparameters...")
  final_model <- ranger(
    x = t(x_combined),
    y = y,
    num.trees = ntree,
    mtry = best_mtry,
    min.node.size = best_min_node_size,
    sample.fraction = best_sample_fraction,
    importance = "permutation",
    num.threads = n_cores,
    classification = is.factor(y),
    probability = is.factor(y) && nlevels(y) == 2
  )

  list(
    performance = final_metrics,
    outer_fold_results = outer_results,
    avg_importance = sort(avg_importance, decreasing = TRUE),
    top_cpgs = names(
      sort(avg_importance, decreasing = TRUE)[seq_len(min(
        100,
        length(avg_importance)
      ))]
    ),
    final_model = final_model,
    best_hyperparameters = list(
      mtry = best_mtry,
      min.node.size = best_min_node_size,
      sample.fraction = best_sample_fraction
    ),
    all_fold_hyperparameters = all_best_params,
    hyperparam_frequencies = list(
      mtry = table(hyperparam_df$mtry),
      min.node.size = table(hyperparam_df$min.node.size),
      sample.fraction = table(hyperparam_df$sample.fraction)
    )
  )
}
