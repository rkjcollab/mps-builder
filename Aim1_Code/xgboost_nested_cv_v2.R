# XGBoost with nested cross-validation
# Matches the structure of elastic_net_nested_cv() and random_forest_nested_cv()
# Uses 5x5 nested CV with inner-loop hyperparameter tuning
# Parallelization via xgboost's nthread parameter
#
# Arguments:
#   X                      CpG beta matrix (CpGs in rows, samples in columns)
#   y                      Response variable (factor for classification, numeric for regression)
#   nested_folds           Nested fold structure from create_nested_folds()
#   covariates             Optional data frame of covariates (samples in rows)
#   hyperparam_grid        Optional data.frame of hyperparameter combos to search
#   nrounds_max            Maximum boosting rounds (early stopping determines actual)
#   early_stopping_rounds  Stop if no improvement after this many rounds
#   n_cores                Number of threads for xgboost (default 4)
#   importance_threshold   Number of top CpGs to flag by gain importance
#   verbose                Print progress messages (default TRUE)

xgboost_nested_cv <- function(X, y,
                              nested_folds,
                              covariates = NULL,
                              hyperparam_grid = NULL,
                              nrounds_max = 1000,
                              early_stopping_rounds = 50,
                              n_cores = 4,
                              importance_threshold = 100,
                              verbose = TRUE) {
  
  # Load required packages
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required. Install with: install.packages('xgboost')")
  }
  library(xgboost)
  
  # Input validation
  if (ncol(X) != length(y)) {
    stop("Number of samples in X (", ncol(X),
         ") doesn't match length of y (", length(y), ")")
  }
  
  if (!is.null(covariates)) {
    if (nrow(covariates) != ncol(X)) {
      stop("Number of rows in covariates must match number of samples in X")
    }
  }
  
  # Determine task type
  is_classification <- is.factor(y)
  if (is_classification) {
    y_numeric <- as.numeric(y) - 1
    objective <- "binary:logistic"
    eval_metric <- "logloss"
    positive_class <- levels(y)[2]
    message("Classification mode: positive class = '", positive_class, "'")
  } else {
    y_numeric <- as.numeric(y)
    objective <- "reg:squarederror"
    eval_metric <- "rmse"
    message("Regression mode")
  }
  
  # Combine CpGs with covariates
  if (!is.null(covariates)) {
    message("Including covariates in XGBoost model")
    cov_matrix <- as.matrix(covariates)
    if (!is.numeric(cov_matrix)) {
      cov_numeric <- apply(covariates, 2, function(col) {
        if (is.character(col) || is.factor(col)) {
          as.numeric(as.factor(col))
        } else {
          as.numeric(col)
        }
      })
      cov_matrix <- as.matrix(cov_numeric)
    }
    X_combined <- rbind(X, t(cov_matrix))
    feature_names <- c(rownames(X), colnames(covariates))
    cpg_names <- rownames(X)
  } else {
    X_combined <- X
    feature_names <- rownames(X)
    cpg_names <- rownames(X)
  }
  
  # Default hyperparameter grid
  if (is.null(hyperparam_grid)) {
    hyperparam_grid <- expand.grid(
      eta              = c(0.01, 0.05, 0.1),
      max_depth        = c(3, 6),
      subsample        = c(0.7, 1.0),
      colsample_bytree = c(0.5, 0.8),
      min_child_weight = c(1, 10),
      gamma            = c(0, 1),
      stringsAsFactors = FALSE
    )
  }
  
  n_combos <- nrow(hyperparam_grid)
  message("Hyperparameter grid: ", n_combos, " combinations to evaluate in inner CV")
  
  # Storage for outer fold results
  n_outer <- length(nested_folds)
  outer_results <- list()
  all_importance_scores <- list()
  all_best_params <- list()
  all_best_nrounds <- numeric(n_outer)
  
  # Progress bar
  if (verbose) {
    pb_outer <- progress::progress_bar$new(
      format = "  XGBoost outer folds [:bar] :percent | Fold :current/:total | ETA: :eta",
      total = n_outer,
      clear = FALSE,
      width = 80
    )
  }
  
  # ~ ~ ~ OUTER LOOP: Performance estimation ~ ~ ~
  for (i in 1:n_outer) {
    message(paste("\n*** Outer fold", i, "of", n_outer, "***"))
    
    # Split into outer train/test
    train_idx  <- nested_folds[[i]]$outer_train
    test_idx   <- nested_folds[[i]]$outer_test
    inner_folds <- nested_folds[[i]]$inner_folds
    
    X_train <- X_combined[, train_idx, drop = FALSE]
    X_test  <- X_combined[, test_idx, drop = FALSE]
    y_train <- y_numeric[train_idx]
    y_test  <- y_numeric[test_idx]
    
    # Create DMatrix objects (xgboost optimized data structure)
    # xgboost expects samples in rows, features in columns
    dtrain <- xgb.DMatrix(data = t(X_train), label = y_train)
    dtest  <- xgb.DMatrix(data = t(X_test),  label = y_test)
    
    # ~ ~ ~ INNER LOOP: Hyperparameter tuning ~ ~ ~
    # Uses manual inner CV with xgb.train instead of xgb.cv,
    # which can segfault with custom folds in some xgboost versions
    best_inner_score <- Inf
    best_params <- NULL
    best_nrounds <- 100
    n_inner <- length(inner_folds)
    
    for (j in 1:n_combos) {
      
      if (verbose && j %% 10 == 0) {
        message("  Inner tuning: combo ", j, "/", n_combos)
      }
      
      # Current hyperparameter set
      params_j <- list(
        objective        = objective,
        eval_metric      = eval_metric,
        eta              = hyperparam_grid$eta[j],
        max_depth        = hyperparam_grid$max_depth[j],
        subsample        = hyperparam_grid$subsample[j],
        colsample_bytree = hyperparam_grid$colsample_bytree[j],
        min_child_weight = hyperparam_grid$min_child_weight[j],
        gamma            = hyperparam_grid$gamma[j],
        nthread          = n_cores
      )
      
      # Manual inner CV: train on (k-1) inner folds, validate on held-out fold
      # Track best nrounds per inner fold via watchlist early stopping
      inner_fold_scores <- numeric(n_inner)
      inner_fold_nrounds <- numeric(n_inner)
      combo_failed <- FALSE
      
      for (k in 1:n_inner) {
        inner_val_idx <- inner_folds[[k]]
        inner_train_idx <- setdiff(train_idx, inner_val_idx)
        
        X_inner_train <- X_combined[, inner_train_idx, drop = FALSE]
        X_inner_val   <- X_combined[, inner_val_idx, drop = FALSE]
        y_inner_train <- y_numeric[inner_train_idx]
        y_inner_val   <- y_numeric[inner_val_idx]
        
        d_inner_train <- xgb.DMatrix(data = t(X_inner_train), label = y_inner_train)
        d_inner_val   <- xgb.DMatrix(data = t(X_inner_val),   label = y_inner_val)
        
        # Train with evals list for early stopping
        inner_model <- tryCatch({
          xgb.train(
            params    = params_j,
            data      = d_inner_train,
            nrounds   = nrounds_max,
            evals     = list(val = d_inner_val),
            early_stopping_rounds = early_stopping_rounds,
            verbose   = 0
          )
        }, error = function(e) {
          message("  Warning: xgb.train failed for combo ", j, " inner fold ", k, ": ", e$message)
          return(NULL)
        })
        
        if (is.null(inner_model)) {
          combo_failed <- TRUE
          break
        }
        
        # Extract best score and iteration from this inner fold
        # Uses xgb.attr() which stores values as strings in the C++ object
        b_iter <- as.numeric(xgb.attr(inner_model, "best_iteration"))
        b_score <- as.numeric(xgb.attr(inner_model, "best_score"))
        inner_fold_nrounds[k] <- b_iter
        inner_fold_scores[k] <- b_score
      }
      
      if (combo_failed) next
      
      # Average validation score across inner folds
      mean_inner_score <- mean(inner_fold_scores)
      mean_inner_nrounds <- round(mean(inner_fold_nrounds))
      
      if (mean_inner_score < best_inner_score) {
        best_inner_score <- mean_inner_score
        best_params <- params_j
        best_nrounds <- mean_inner_nrounds
      }
    }
    
    message("  Best inner CV ", eval_metric, ": ", round(best_inner_score, 4),
            " | nrounds: ", best_nrounds,
            " | eta: ", best_params$eta,
            " | max_depth: ", best_params$max_depth)
    
    # Fit final model on outer training set with best hyperparameters
    final_model_fold <- xgb.train(
      params    = best_params,
      data      = dtrain,
      nrounds   = best_nrounds,
      verbose   = 0
    )
    
    # Predictions on outer test set
    pred_raw <- predict(final_model_fold, dtest)
    
    if (is_classification) {
      pred_class <- factor(
        ifelse(pred_raw > 0.5, positive_class, levels(y)[1]),
        levels = levels(y)
      )
      y_test_factor <- factor(
        ifelse(y_test == 1, positive_class, levels(y)[1]),
        levels = levels(y)
      )
      prob_vec <- pred_raw
    } else {
      predictions <- pred_raw
    }
    
    # Compute fold performance metrics
    if (is_classification) {
      acc <- mean(pred_class == y_test_factor)
      
      auc_val <- tryCatch({
        roc_obj <- pROC::roc(y_test, prob_vec, quiet = TRUE)
        as.numeric(pROC::auc(roc_obj))
      }, error = function(e) NA_real_)
      
      cm <- table(Predicted = pred_class, Actual = y_test_factor)
      
      fold_metrics <- list(
        Accuracy = acc,
        AUC = auc_val,
        ConfusionMatrix = cm
      )
    } else {
      residuals <- y_test - predictions
      fold_metrics <- list(
        RMSE = sqrt(mean(residuals^2)),
        MAE  = mean(abs(residuals)),
        R2   = 1 - sum(residuals^2) / sum((y_test - mean(y_test))^2),
        Cor  = cor(predictions, y_test)
      )
    }
    
    # Feature importance (gain-based)
    importance_matrix <- xgb.importance(
      feature_names = feature_names,
      model = final_model_fold
    )
    
    all_importance_scores[[i]] <- importance_matrix
    
    # Store fold results
    outer_results[[i]] <- list(
      metrics       = fold_metrics,
      predictions   = if (is_classification) pred_class else predictions,
      probabilities = if (is_classification) prob_vec else NULL,
      actual        = y_test,
      actual_factor = if (is_classification) y_test_factor else NULL,
      best_params   = best_params,
      best_nrounds  = best_nrounds,
      best_inner_score = best_inner_score,
      model         = final_model_fold
    )
    
    all_best_params[[i]] <- data.frame(
      eta              = best_params$eta,
      max_depth        = best_params$max_depth,
      subsample        = best_params$subsample,
      colsample_bytree = best_params$colsample_bytree,
      min_child_weight = best_params$min_child_weight,
      gamma            = best_params$gamma,
      nrounds          = best_nrounds,
      inner_score      = best_inner_score,
      stringsAsFactors = FALSE
    )
    
    all_best_nrounds[i] <- best_nrounds
    
    if (verbose) pb_outer$tick()
  }
  
  # ~ ~ ~ AGGREGATE RESULTS ACROSS OUTER FOLDS ~ ~ ~
  message("\n*** Aggregating results across outer folds ***")
  
  if (is_classification) {
    final_metrics <- data.frame(
      Accuracy_Mean = mean(sapply(outer_results, function(x) x$metrics$Accuracy)),
      Accuracy_SD   = sd(sapply(outer_results, function(x) x$metrics$Accuracy)),
      AUC_Mean      = mean(sapply(outer_results, function(x) x$metrics$AUC), na.rm = TRUE),
      AUC_SD        = sd(sapply(outer_results, function(x) x$metrics$AUC), na.rm = TRUE)
    )
  } else {
    final_metrics <- data.frame(
      RMSE_Mean = mean(sapply(outer_results, function(x) x$metrics$RMSE)),
      RMSE_SD   = sd(sapply(outer_results, function(x) x$metrics$RMSE)),
      MAE_Mean  = mean(sapply(outer_results, function(x) x$metrics$MAE)),
      MAE_SD    = sd(sapply(outer_results, function(x) x$metrics$MAE)),
      R2_Mean   = mean(sapply(outer_results, function(x) x$metrics$R2)),
      R2_SD     = sd(sapply(outer_results, function(x) x$metrics$R2)),
      Cor_Mean  = mean(sapply(outer_results, function(x) x$metrics$Cor)),
      Cor_SD    = sd(sapply(outer_results, function(x) x$metrics$Cor))
    )
  }
  
  # Aggregate feature importance across folds
  all_gains <- list()
  for (i in 1:n_outer) {
    imp_i <- all_importance_scores[[i]]
    gains_i <- setNames(imp_i$Gain, imp_i$Feature)
    all_gains[[i]] <- gains_i
  }
  
  all_feature_names <- unique(unlist(lapply(all_gains, names)))
  
  gain_matrix <- matrix(0, nrow = length(all_feature_names), ncol = n_outer,
                        dimnames = list(all_feature_names, paste0("Fold", 1:n_outer)))
  for (i in 1:n_outer) {
    matched <- match(names(all_gains[[i]]), all_feature_names)
    gain_matrix[matched, i] <- all_gains[[i]]
  }
  
  avg_gain <- rowMeans(gain_matrix)
  avg_gain <- sort(avg_gain, decreasing = TRUE)
  
  appearance_count <- rowSums(gain_matrix > 0)
  
  cpg_gain <- avg_gain[names(avg_gain) %in% cpg_names]
  
  if (!is.null(covariates)) {
    cov_gain <- avg_gain[names(avg_gain) %in% colnames(covariates)]
    message("Covariate importance (avg gain): ")
    for (nm in names(cov_gain)) {
      message("  ", nm, ": ", round(cov_gain[nm], 6))
    }
  } else {
    cov_gain <- NULL
  }
  
  top_cpgs <- head(names(cpg_gain), importance_threshold)
  
  stable_cpgs <- names(which(appearance_count[cpg_names[cpg_names %in% names(appearance_count)]] == n_outer))
  
  message("Top ", importance_threshold, " CpGs identified by average gain")
  message("CpGs appearing in all ", n_outer, " folds: ", length(stable_cpgs))
  
  # Aggregate hyperparameter selections
  params_df <- do.call(rbind, all_best_params)
  
  get_mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  
  consensus_params <- list(
    objective        = objective,
    eval_metric      = eval_metric,
    eta              = get_mode(params_df$eta),
    max_depth        = get_mode(params_df$max_depth),
    subsample        = get_mode(params_df$subsample),
    colsample_bytree = get_mode(params_df$colsample_bytree),
    min_child_weight = get_mode(params_df$min_child_weight),
    gamma            = get_mode(params_df$gamma),
    nthread          = n_cores
  )
  consensus_nrounds <- round(median(all_best_nrounds))
  
  message("\nConsensus hyperparameters (mode across folds):")
  message("  eta: ", consensus_params$eta,
          " | max_depth: ", consensus_params$max_depth,
          " | subsample: ", consensus_params$subsample,
          " | colsample_bytree: ", consensus_params$colsample_bytree)
  message("  min_child_weight: ", consensus_params$min_child_weight,
          " | gamma: ", consensus_params$gamma,
          " | nrounds: ", consensus_nrounds)
  
  # Train final model on all data
  message("\nTraining final model on all data with consensus hyperparameters...")
  
  dfull <- xgb.DMatrix(data = t(X_combined), label = y_numeric)
  
  final_model <- xgb.train(
    params  = consensus_params,
    data    = dfull,
    nrounds = consensus_nrounds,
    verbose = 0
  )
  
  final_importance <- xgb.importance(
    feature_names = feature_names,
    model = final_model
  )
  
  # Return results
  message("\n*** XGBoost nested CV complete ***")
  if (is_classification) {
    message("Mean Accuracy: ", round(final_metrics$Accuracy_Mean, 4),
            " +/- ", round(final_metrics$Accuracy_SD, 4))
    message("Mean AUC: ", round(final_metrics$AUC_Mean, 4),
            " +/- ", round(final_metrics$AUC_SD, 4))
  } else {
    message("Mean RMSE: ", round(final_metrics$RMSE_Mean, 4),
            " +/- ", round(final_metrics$RMSE_SD, 4))
    message("Mean R^2: ", round(final_metrics$R2_Mean, 4),
            " +/- ", round(final_metrics$R2_SD, 4))
    message("Mean Cor: ", round(final_metrics$Cor_Mean, 4),
            " +/- ", round(final_metrics$Cor_SD, 4))
  }
  
  return(list(
    performance = final_metrics,
    outer_fold_results = outer_results,
    avg_importance = cpg_gain,
    importance_matrix = gain_matrix,
    top_cpgs = top_cpgs,
    stable_cpgs = stable_cpgs,
    cpg_appearance_count = appearance_count[intersect(names(appearance_count), cpg_names)],
    covariate_importance = cov_gain,
    fold_hyperparameters = params_df,
    consensus_params = consensus_params,
    consensus_nrounds = consensus_nrounds,
    final_model = final_model,
    final_importance = final_importance
  ))
}
