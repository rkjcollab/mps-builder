# Step 2: nested cross-validation setup.
# Moved from Aim1_Code/methylation_pipeline_v4_shared.R.

#' Create nested cross-validation fold structure
#'
#' Generates indices for nested cross-validation, where outer folds are used
#' for unbiased performance estimation and inner folds for hyperparameter
#' tuning. Uses stratified sampling for classification problems.
#'
#' @param y Response variable (factor for classification, numeric for
#' regression)
#' @param k_outer Number of outer folds (default 5)
#' @param k_inner Number of inner folds (default 5)
#' @param seed Random seed for reproducibility
#' @return List of nested fold structures with outer_test, outer_train, and
#' inner_folds
#' @export
create_nested_folds <- function(y, k_outer = 5, k_inner = 5, seed = 123) {
  set.seed(seed)

  # Create stratified folds for classification to maintain class balance
  # Use regular k-fold for regression
  if (is.factor(y)) {
    outer_folds <- createFolds(
      y,
      k = k_outer, list = TRUE, returnTrain = FALSE
    )
  } else {
    outer_folds <- createFolds(
      seq_along(y),
      k = k_outer, list = TRUE, returnTrain = FALSE
    )
  }

  # For each outer fold, create inner folds from the training subset
  nested_folds <- list()
  for (i in 1:k_outer) {
    train_indices <- setdiff(seq_along(y), outer_folds[[i]])
    y_train <- y[train_indices]

    # Create inner folds (stratified for classification)
    if (is.factor(y_train)) {
      inner_folds <- createFolds(
        y_train,
        k = k_inner, list = TRUE, returnTrain = FALSE
      )
    } else {
      inner_folds <- createFolds(seq_along(y_train),
        k = k_inner,
        list = TRUE, returnTrain = FALSE
      )
    }

    # Convert inner fold indices back to original data indices
    inner_folds <- lapply(inner_folds, function(idx) train_indices[idx])

    nested_folds[[i]] <- list(
      outer_test = outer_folds[[i]],
      outer_train = train_indices,
      inner_folds = inner_folds
    )
  }

  nested_folds
}

#' Convert nested fold structure to glmnet foldid format
#'
#' @param inner_folds List of inner fold indices (in original data coordinates)
#' @param train_idx Vector of training indices for the current outer fold
#' @return Integer vector where each element indicates which inner fold that
#' sample belongs to
#' @keywords internal
create_foldid <- function(inner_folds, train_idx) {
  n_total <- length(train_idx)
  foldid <- rep(NA_integer_, n_total)

  for (i in seq_along(inner_folds)) {
    # Indices in original data coordinates
    idx_orig <- inner_folds[[i]]
    # Convert to training set coordinates (1:n_train)
    idx_train <- match(idx_orig, train_idx)
    foldid[idx_train] <- i
  }
  foldid
}
