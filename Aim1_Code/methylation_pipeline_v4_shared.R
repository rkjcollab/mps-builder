## Set working directory ##
setwd(paste0(Sys.getenv("RKJCOLLAB"), "/Project_Folder"))

library(caret)
library(glmnet)
library(randomForest)
library(xgboost)
library(limma)
library(matrixStats)
library(dplyr)
library(readr)
library(readxl)
library(pROC)
library(doParallel)
library(ranger)
library(VennDiagram)
library(grid)
library(progress)
library(purrr)

# ============================================
# PARALLEL PROCESSING CONFIGURATION
# ============================================
# Enable parallel processing using half of available cores
# This balances performance with system responsiveness
n_cores <- max(1, floor(detectCores() / 2))
registerDoParallel(cores = n_cores)

# ============================================
# LOADING IN DATA
# ============================================

# beta_raw <- read_csv("InvestigatorPackageCompressedCBMethylation/CB_m_norm_batch_plate.csv")
# saveRDS(beta_matrix, file = "InvestigatorPackageCompressedCBMethylation/beta_raw.rds")
beta_raw <- readRDS("InvestigatorPackageCompressedCBMethylation/beta_raw.rds")
sample_ids <- beta_raw %>% pull(pid)
beta_raw <- beta_raw %>% select(-pid)
beta_matrix <- t(as.matrix(beta_raw))
colnames(beta_matrix) <- sample_ids
nrow(beta_matrix) 
ncol(beta_matrix)

# Convert M-values to beta-values
beta_matrix <- 2^beta_matrix / (2^beta_matrix + 1)
range(beta_matrix)
beta_matrix <- beta_matrix[order(rownames(beta_matrix)), ]
beta_matrix <- beta_matrix[, order(colnames(beta_matrix))]

all_cpgs <- rownames(beta_matrix)

prev_ewas <- read_excel("previous_EWAS_CpGs.xlsx")
literature_studies <- split(prev_ewas$cpg, prev_ewas$study)
literature_studies <- lapply(literature_studies, unique)

hs_data <- read_csv("p0001_HohsfieldKirk_150210.csv")

hs_cell_proportions <- read_csv(file = "InvestigatorPackageCompressedCBMethylation/DecompressedCBMethylation/HS_CB_cell_proportions_IDOL.csv")
hs_cell_proportions <- hs_cell_proportions %>%
  mutate(CD8T = ifelse(CD8T == 0, min(CD8T[CD8T > 0])/2, CD8T),
         Bcell = ifelse(Bcell == 0, min(Bcell[Bcell > 0])/2, Bcell),
         Mono = ifelse(Mono == 0, min(Mono[Mono > 0])/2, Mono),
         nRBC = ifelse(nRBC == 0, min(nRBC[nRBC > 0])/2, nRBC)) %>%
  arrange(Sample_Name)

hs_data <- hs_data %>%
  left_join((hs_cell_proportions %>% select(-Sample_Plate, -rgName, -Gran)), 
            by = c("childID" = "Sample_Name")) %>%
  filter(as.character(childID) %in% colnames(beta_matrix)) #,
         # !is.na(pm25_fullpreg_d)) %>%
  mutate(pre_preg_bmi = ifelse(is.na(pre_preg_bmi), median(pre_preg_bmi, na.rm = TRUE), pre_preg_bmi))


beta_matrix <- beta_matrix[, as.character(hs_data$childID)]

# Subset methylation matrix to EPIC-compatible CpGs
beta_matrix_epic <- beta_matrix[rownames(beta_matrix) %in% epic_cpgs, ]

cat("Original CpGs:", nrow(beta_matrix), "\n")
cat("EPIC-compatible CpGs:", nrow(beta_matrix_epic), "\n")

# Subset beta matrix for children with ancestry PC data
beta_matrix_epic_aPC <- beta_matrix_epic[, colnames(beta_matrix_epic) %in% as.character(hs_data_pcs_childID$childID)]

# Rerun nested CV with same fold structure for comparability
results_epic_compatible <- run_nested_cv(
  meth_matrix_epic, 
  pm25_contin,  # or pm25_residualized
  nested_folds,  # same folds as original
  ...
)

pm25_contin <- hs_data$pm25_fullpreg_d
names(pm25_contin) <- hs_data$childID
pm25_contin <- pm25_contin[colnames(beta_matrix)]

pm25_data_new <- read.csv(file = "PM25_pregnancy_exposures.csv")
pm25_contin_new <- pm25_data_new$pm25_full_pregnancy
names(pm25_contin_new) <- pm25_data_new$childID
pm25_contin_new <- pm25_contin_new[colnames(beta_matrix)]

pm25_data_highres <- read.csv(file = "PM25_pregnancy_exposures_highres.csv")
pm25_contin_highres <- pm25_data_highres$pm25_full
names(pm25_contin_highres) <- pm25_data_highres$childID
pm25_contin_highres <- pm25_contin_highres[colnames(beta_matrix)]

# ============================================
# STEP 0: TESTING
# ============================================
set.seed(123)
sample_idx <- sample(1:ncol(beta_matrix), size = 100)
cpg_idx <- sample(1:nrow(beta_matrix), size = 5000)

beta_subset <- beta_matrix[cpg_idx, sample_idx]

test_lit_ind

pm25_subset <- hs_data_subset$pm25_fullpreg_d
names(pm25_subset) <- hs_data_subset$childID
pm25_subset <- pm25_subset[colnames(beta_subset)]

test_cpg_summary <- data.frame(
  CpG = rownames(test_filtered_cpgs_full$criteria_matrix),
  Variance = test_filtered_cpgs_full$criteria_matrix[, "Variance"],
  Association = test_filtered_cpgs_full$criteria_matrix[, "Association"],
  Literature = test_filtered_cpgs_full$criteria_matrix[, "Literature"],
  Total_Criteria = test_filtered_cpgs_full$criteria_counts,
  Selected = rownames(test_filtered_cpgs_full$criteria_matrix) %in% 
    test_filtered_cpgs_full$selected_cpgs
)

covariates <- hs_data %>%
  mutate(GA_at_birth_in_years = GA_at_birth_in_days/365.25,
         ChildSeasonBirth = factor(ChildSeasonBirth),
         gestsmoking = factor(gestsmoking),
         infant_sex = factor(infant_sex),
         mat_edu = factor(epi5)) %>%
  select(temp_fullpreg_d, ChildSeasonBirth, GA_at_birth_in_years, 
         gestsmoking, pre_preg_bmi, infant_sex, mat_edu, hs1_medianincome, 
         CD8T, CD4T, NK, Bcell, Mono, nRBC) %>%
  data.frame()
rownames(covariates) <- as.character(hs_data$childID)

complete_idx <- complete.cases(covariates)

cpg_filtering_sensitivity_analysis <- sensitivity_analysis_filtering(beta_matrix = beta_matrix, 
                                                                     pm25_data = pm25_contin, 
                                                                     literature_studies = literature_studies)

cpg_filtering_sensitivity_analysis_dichot <- sensitivity_analysis_filtering(beta_matrix = beta_matrix,
                                                                            pm25_data = pm25_dichot,
                                                                            literature_studies = literature_studies)

# Define parameter grid
param_grid <- expand.grid(
  var_percentile = c(70, 75, 80, 85, 90),
  pval_threshold = c(0.01, 0.05, 0.10),
  min_criteria = c(1, 2)
)

# Run all filtering steps
filtered_results_contin <- pmap(param_grid, function(var_percentile, pval_threshold, min_criteria) {
  filter_cpg_sites(
    beta_matrix = beta_matrix,
    pm25_data = pm25_contin,
    literature_studies = literature_studies,
    var_percentile = var_percentile,
    pval_threshold = pval_threshold,
    min_criteria = min_criteria
  )
})

# Name the list
names(filtered_results_contin) <- paste0(
  param_grid$var_percentile, "_",
  sprintf("%03d", param_grid$pval_threshold * 100), "_",
  param_grid$min_criteria
)

nested_folds_contin <- create_nested_folds(
  y = pm25_contin,
  k_outer = k_outer,
  k_inner = k_inner
)

# Check covariate associations
cor.test(covariates$temp_fullpreg_d, pm25_contin)
cor.test(covariates$GA_at_birth_in_years, pm25_contin)
cor.test(covariates$pre_preg_bmi, pm25_contin)

cor.test(covariates$CD8T, pm25_contin)
cor.test(covariates$CD4T, pm25_contin)
cor.test(covariates$NK, pm25_contin)
cor.test(covariates$Bcell, pm25_contin)
cor.test(covariates$Mono, pm25_contin)
cor.test(covariates$nRBC, pm25_contin)

summary(aov(pm25_contin ~ covariates$ChildSeasonBirth))
summary(aov(pm25_contin ~ covariates$gestsmoking))
summary(aov(pm25_contin ~ covariates$infant_sex))
summary(aov(pm25_contin ~ covariates$mat_edu))
summary(aov(pm25_contin ~ covariates$hs1_medianincome))


# Run elastic net on all filtered results
en_results_contin_min1 <- map(filtered_results_contin[param_grid$min_criteria == 1 & (param_grid$var_percentile %in% c(80, 85, 90))], function(filtered) {
  elastic_net_nested_cv(
    X = filtered$filtered_beta,
    y = pm25_contin,
    covariates = covariates,
    nested_folds = nested_folds_contin,
    alpha_range = alpha_range
  )
})


en_results_contin_min2 <- map(filtered_results_contin[param_grid$min_criteria == 2], function(filtered) {
  elastic_net_nested_cv(
    X = filtered$filtered_beta,
    y = pm25_contin,
    covariates = covariates,
    nested_folds = nested_folds_contin,
    alpha_range = alpha_range
  )
})


en_results_contin_all <- c(en_results_contin_min1, en_results_contin_min2)

# Extract performance metrics from all models
performance_comparison <- map_dfr(names(en_results_contin_all), function(name) {
  # Parse the name to get parameters
  parts <- strsplit(name, "_")[[1]]
  
  data.frame(
    config = name,
    var_percentile = as.numeric(parts[1]),
    pval_threshold = as.numeric(parts[2]) / 100,  # Convert 001 -> 0.01
    min_criteria = as.numeric(parts[3]),
    en_results_contin_all[[name]]$performance,
    n_cpgs_selected = length(en_results_contin_all[[name]]$selected_cpgs),
    n_stable_cpgs = length(en_results_contin_all[[name]]$stable_cpgs),
    avg_features_per_fold = en_results_contin_all[[name]]$performance$Avg_Features
  )
})


# Baseline with nested CV (same folds as elastic net)
baseline_cv_results <- list()

for(i in 1:length(nested_folds_contin)) {
  train_idx <- nested_folds_contin[[i]]$outer_train
  test_idx <- nested_folds_contin[[i]]$outer_test
  
  # Train covariates-only model
  baseline_model <- lm(pm25_contin[train_idx] ~ ., 
                       data = covariates[train_idx, , drop = FALSE])
  
  # Predict on test set
  predictions <- predict(baseline_model, 
                         newdata = covariates[test_idx, , drop = FALSE])
  
  # Calculate metrics
  actual <- pm25_contin[test_idx]
  rmse <- sqrt(mean((predictions - actual)^2))
  r2 <- 1 - sum((actual - predictions)^2) / sum((actual - mean(actual))^2)
  
  baseline_cv_results[[i]] <- list(rmse = rmse, r2 = r2)
}

# Average across folds
baseline_cv_rmse <- mean(sapply(baseline_cv_results, function(x) x$rmse))
baseline_cv_r2 <- mean(sapply(baseline_cv_results, function(x) x$r2))

cat("\n=== Baseline with Nested CV (Covariates Only) ===\n")
cat("CV RMSE:", baseline_cv_rmse, "\n")
cat("CV R²:", baseline_cv_r2, "\n\n")

cat("=== Your Best Model vs CV Baseline ===\n")
cat("Best model RMSE:", en_results_contin_min2$`70_005_2`$performance$RMSE_Mean, "\n")
cat("Baseline CV RMSE:", baseline_cv_rmse, "\n")
cat("Improvement:", round((baseline_cv_rmse - en_results_contin_min2$`70_005_2`$performance$RMSE_Mean) / baseline_cv_rmse * 100, 1), "%\n")
cat("\nBest model R²:", en_results_contin_min2$`70_005_2`$performance$R2_Mean, "\n")
cat("Baseline CV R²:", baseline_cv_r2, "\n")
cat("Added R² from CpGs:", round(en_results_contin_min2$`70_005_2`$performance$R2_Mean - baseline_cv_r2, 3), "\n")


# Calculate profile score
en_contin_risk_score <- create_methylation_risk_score(
  beta_matrix = filtered_results_contin[["70_005_2"]]$filtered_beta,
  model_results = en_results_contin_min2[["70_005_2"]],
  method = "elastic_net",
  score_type = "weighted"
)

# Validate
cor.test(en_contin_risk_score, pm25_contin)




# ============================================
# STEP 1: CpG SITE MULTI-CRITERIA FILTERING
# ============================================

#' Create indicator variables for literature-based CpG filtering
#' 
#' @param all_cpgs Character vector of all CpG names in the dataset
#' @param literature_studies Named list where each element is a character vector 
#'        of CpG names from a published study (e.g., Gruzieva, Wang, Parikh, Broseus)
#' @return List containing: overall_indicator (binary), by_study (matrix of indicators),
#'         n_studies (count per CpG), literature_cpgs (all unique literature CpGs)
create_literature_indicator <- function(all_cpgs, literature_studies) {
  
  # Combine all CpGs from literature into a single unique set
  all_literature_cpgs <- unique(unlist(literature_studies))
  
  # Create binary indicator for each study (1 = CpG present in study)
  study_indicators <- sapply(literature_studies, function(study_cpgs) {
    as.integer(all_cpgs %in% study_cpgs)
  })
  
  # Create overall indicator (1 = CpG appears in at least one study)
  in_any_study <- as.integer(all_cpgs %in% all_literature_cpgs)
  
  # Count number of studies each CpG appears in (replication score)
  n_studies <- rowSums(study_indicators)
  
  return(list(
    overall_indicator = in_any_study,
    by_study = study_indicators,
    n_studies = n_studies,
    literature_cpgs = all_literature_cpgs
  ))
}

#' Filter CpG sites using three criteria: variance, association, and literature support
#' 
#' This function implements a multi-criteria filtering approach to reduce the 
#' dimensionality of methylation data before machine learning. CpGs must meet
#' at least min_criteria of the three filtering criteria to be retained.
#' 
#' @param beta_matrix Numeric matrix of beta values (CpGs in rows, samples in columns)
#' @param pm25_data Numeric vector (continuous) or factor (dichotomous) of PM2.5 exposure
#' @param literature_studies Optional named list of CpG vectors from published studies
#' @param var_percentile Percentile threshold for variance filter (default 80 = top 20%)
#' @param pval_threshold P-value threshold for association filter (default 0.05)
#' @param min_criteria Minimum number of criteria a CpG must meet (default 2)
#' @return List containing filtered beta matrix, selected CpGs, and filtering statistics
filter_cpg_sites <- function(beta_matrix, pm25_data,
                             literature_studies = NULL,
                             var_percentile = 80,
                             pval_threshold = 0.05,
                             min_criteria = 2) {
  
  message("Starting multi-criteria CpG filtering...")
  message(paste("Input CpGs:", nrow(beta_matrix)))
  
  all_cpgs <- rownames(beta_matrix)
  n_cpgs <- nrow(beta_matrix)
  
  # Initialize matrix to track which criteria each CpG meets (rows=CpGs, cols=criteria)
  criteria_matrix <- matrix(0, nrow = n_cpgs, ncol = 3,
                            dimnames = list(all_cpgs, 
                                            c("Variance", "Association", "Literature")))
  
  # ============================================
  # CRITERION 1: Variance-based filtering
  # Rationale: CpGs with low variance across samples are unlikely to be 
  # informative for predicting exposure differences
  # ============================================
  message("\n--- Criterion 1: Variance-based filtering ---")
  cpg_variance <- rowVars(beta_matrix, na.rm = TRUE)
  var_threshold <- quantile(cpg_variance, var_percentile/100, na.rm = TRUE)
  criteria_matrix[, "Variance"] <- as.integer(cpg_variance > var_threshold)
  n_variance <- sum(criteria_matrix[, "Variance"])
  message(paste("CpGs passing variance filter (top", 100-var_percentile, "%): ", 
                n_variance, " (", round(100*n_variance/n_cpgs, 1), "%)", sep=""))
  
  # ============================================
  # CRITERION 2: Association-based filtering
  # Rationale: Retain CpGs that show preliminary association with PM2.5
  # Uses correlation test for continuous outcomes, t-test for dichotomous
  # ============================================
  message("\n--- Criterion 2: Association-based filtering ---")
  
  # Create progress bar for p-value calculations
  pb <- progress_bar$new(
    format = "  Computing p-values [:bar] :percent | ETA: :eta | :current/:total CpGs",
    total = n_cpgs,
    clear = FALSE,
    width = 80
  )
  
  # Calculate association p-values with progress tracking
  # Using parallel processing via foreach for large datasets
  if (n_cpgs > 1000 && n_cores > 1) {
    # Parallel computation for large datasets
    message(paste("  Using parallel computation with", n_cores, "cores..."))
    
    # Split CpGs into chunks for parallel processing
    chunk_size <- ceiling(n_cpgs / n_cores)
    chunks <- split(1:n_cpgs, ceiling(seq_along(1:n_cpgs) / chunk_size))
    
    pvals <- foreach(chunk = chunks, .combine = c, .packages = c("stats")) %dopar% {
      chunk_pvals <- sapply(chunk, function(i) {
        cpg_values <- beta_matrix[i, ]
        if(is.numeric(pm25_data)) {
          tryCatch(cor.test(cpg_values, pm25_data)$p.value,  # Equivalent to lm() test for p-value unless covariates
                   error = function(e) 1)                    # are present in the model. Then, use lm().
        } else {
          tryCatch(t.test(cpg_values ~ pm25_data)$p.value,
                   error = function(e) 1)
        }
      })
      chunk_pvals
    }
    names(pvals) <- all_cpgs
    pb$terminate()
    
  } else {
    # Sequential computation with progress bar for smaller datasets
    pvals <- numeric(n_cpgs)
    names(pvals) <- all_cpgs
    
    for(i in 1:n_cpgs) {
      cpg_values <- beta_matrix[i, ]
      if(is.numeric(pm25_data)) {
        pvals[i] <- tryCatch(cor.test(cpg_values, pm25_data)$p.value,
                             error = function(e) 1)
      } else {
        pvals[i] <- tryCatch(t.test(cpg_values ~ pm25_data)$p.value,
                             error = function(e) 1)
      }
      pb$tick()
    }
  }
  
  criteria_matrix[, "Association"] <- as.integer(pvals < pval_threshold)
  n_association <- sum(criteria_matrix[, "Association"])
  message(paste("CpGs passing association filter (p < ", pval_threshold, "): ",
                n_association, " (", round(100*n_association/n_cpgs, 1), "%)", sep=""))
  
  # ============================================
  # CRITERION 3: Literature-based filtering
  # Rationale: Prioritize CpGs previously associated with air pollution
  # in published epigenome-wide association studies
  # ============================================
  message("\n--- Criterion 3: Literature-based filtering ---")
  if(!is.null(literature_studies)) {
    lit_indicators <- create_literature_indicator(all_cpgs, literature_studies)
    criteria_matrix[, "Literature"] <- lit_indicators$overall_indicator
    n_literature <- sum(criteria_matrix[, "Literature"])
    message(paste("CpGs in external studies: ", n_literature, 
                  " (", round(100*n_literature/n_cpgs, 1), "%)", sep=""))
    
    # Report breakdown by individual study
    for(study in names(literature_studies)) {
      n_study <- sum(lit_indicators$by_study[, study])
      message(paste("  -", study, ":", n_study, "CpGs"))
    }
  } else {
    message("No literature studies provided - criterion not used")
    criteria_matrix[, "Literature"] <- 0
    n_literature <- 0
    lit_indicators <- NULL
  }
  
  # ============================================
  # COMBINE CRITERIA
  # Select CpGs meeting at least min_criteria of the three filters
  # ============================================
  message(paste("\n--- Combining criteria (at least", min_criteria, "of 3) ---"))
  
  # Count how many criteria each CpG meets
  criteria_counts <- rowSums(criteria_matrix)
  
  # Select CpGs meeting minimum threshold
  selected <- which(criteria_counts >= min_criteria)
  selected_cpgs <- all_cpgs[selected]
  
  # Create breakdown of which criteria combinations CpGs meet
  combination_breakdown <- table(apply(criteria_matrix, 1, function(x) {
    paste(c("Var", "Assoc", "Lit")[as.logical(x)], collapse="+")
  }))
  
  # Report selection statistics
  message("\nCpG selection breakdown:")
  message(paste("  - All 3 criteria:", sum(criteria_counts == 3)))
  message(paste("  - Exactly 2 criteria:", sum(criteria_counts == 2)))
  if(min_criteria <= 1) {
    message(paste("  - Exactly 1 criterion:", sum(criteria_counts == 1)))
  }
  message(paste("  - No criteria:", sum(criteria_counts == 0)))
  
  # Report specific two-way and three-way combinations
  message("\nSpecific combinations (2+ criteria):")
  for(combo in names(combination_breakdown)) {
    n_combo <- combination_breakdown[combo]
    if(grepl("\\+", combo) && n_combo > 0) {
      message(paste("  -", combo, ":", n_combo, "CpGs"))
    }
  }
  
  message(paste("\nTotal CpGs selected:", length(selected_cpgs), 
                "(", round(100*length(selected_cpgs)/n_cpgs, 1), "%)"))
  
  # Compute pairwise and three-way overlap statistics
  overlap_var_assoc <- sum(criteria_matrix[, "Variance"] & 
                             criteria_matrix[, "Association"])
  overlap_var_lit <- sum(criteria_matrix[, "Variance"] & 
                           criteria_matrix[, "Literature"])
  overlap_assoc_lit <- sum(criteria_matrix[, "Association"] & 
                             criteria_matrix[, "Literature"])
  overlap_all <- sum(criteria_counts == 3)
  
  return(list(
    filtered_beta = beta_matrix[selected, ],
    selected_cpgs = selected_cpgs,
    n_selected = length(selected_cpgs),
    criteria_matrix = criteria_matrix,
    criteria_counts = criteria_counts,
    p_values = pvals,
    variance_values = cpg_variance,
    literature_info = lit_indicators,
    filters_summary = list(
      n_variance_only = sum(criteria_counts == 1 & criteria_matrix[, "Variance"]),
      n_association_only = sum(criteria_counts == 1 & criteria_matrix[, "Association"]),
      n_literature_only = sum(criteria_counts == 1 & criteria_matrix[, "Literature"]),
      n_var_assoc = overlap_var_assoc,
      n_var_lit = overlap_var_lit,
      n_assoc_lit = overlap_assoc_lit,
      n_all_three = overlap_all,
      n_at_least_two = sum(criteria_counts >= 2),
      n_total = length(selected_cpgs)
    ),
    combination_breakdown = combination_breakdown
  ))
}


## Additional (optional) functions for analyzing filtering results ##
#' Generate detailed report of CpG filtering results
#' 
#' @param filter_results Output from filter_cpg_sites()
#' @param beta_matrix Optional original beta matrix for variance statistics
analyze_filtering_results <- function(filter_results, beta_matrix = NULL) {
  
  cat("\n========== FILTERING ANALYSIS REPORT ==========\n\n")
  
  # Overall selection statistics
  cat("OVERALL SELECTION:\n")
  cat(sprintf("  Total CpGs selected: %d (%.1f%%)\n", 
              filter_results$n_selected,
              100 * filter_results$n_selected / nrow(filter_results$criteria_matrix)))
  
  # Individual criteria performance
  cat("\nINDIVIDUAL CRITERIA:\n")
  for(crit in colnames(filter_results$criteria_matrix)) {
    n_pass <- sum(filter_results$criteria_matrix[, crit])
    pct_pass <- 100 * n_pass / nrow(filter_results$criteria_matrix)
    cat(sprintf("  %s: %d CpGs (%.1f%%)\n", crit, n_pass, pct_pass))
  }
  
  # Pairwise and three-way overlap statistics
  cat("\nOVERLAP ANALYSIS:\n")
  cat(sprintf("  Variance ∩ Association: %d\n", filter_results$filters_summary$n_var_assoc))
  cat(sprintf("  Variance ∩ Literature: %d\n", filter_results$filters_summary$n_var_lit))
  cat(sprintf("  Association ∩ Literature: %d\n", filter_results$filters_summary$n_assoc_lit))
  cat(sprintf("  All three criteria: %d\n", filter_results$filters_summary$n_all_three))
  
  # Distribution of criteria counts across all CpGs
  cat("\nCRITERIA COUNT DISTRIBUTION:\n")
  counts_table <- table(filter_results$criteria_counts)
  for(i in names(counts_table)) {
    cat(sprintf("  %s criteria: %d CpGs\n", i, counts_table[i]))
  }
  
  # Variance statistics for selected CpGs (if beta matrix provided)
  if(!is.null(beta_matrix)) {
    selected_var <- filter_results$variance_values[filter_results$selected_cpgs]
    if (length(selected_var) > 0) {
      cat("\nSELECTED CPGs VARIANCE STATISTICS:\n")
      cat(sprintf("  Min variance: %.4f\n", min(selected_var)))
      cat(sprintf("  Median variance: %.4f\n", median(selected_var)))
      cat(sprintf("  Mean variance: %.4f\n", mean(selected_var)))
      cat(sprintf("  Max variance: %.4f\n", max(selected_var)))
    } else {
      cat("\nSELECTED CPGs VARIANCE STATISTICS:\n")
      cat("  No CpGs selected; variance statistics not available.\n")
    }
  }
  
  cat("\n================================================\n")
}

#' Export filtering results to files for downstream analysis
#' 
#' @param filter_results Output from filter_cpg_sites()
#' @param output_dir Directory to save results (created if doesn't exist)
export_filtering_results <- function(filter_results, output_dir = "filtering_results") {
  
  dir.create(output_dir, showWarnings = FALSE)
  
  # Create comprehensive data frame with all filtering information
  export_df <- data.frame(
    CpG = rownames(filter_results$criteria_matrix),
    filter_results$criteria_matrix,
    Criteria_Count = filter_results$criteria_counts,
    Selected = as.integer(rownames(filter_results$criteria_matrix) %in% 
                            filter_results$selected_cpgs),
    P_Value = filter_results$p_values,
    Variance = filter_results$variance_values
  )
  
  # Add per-study indicators if literature data available
  if(!is.null(filter_results$literature_info)) {
    export_df <- cbind(export_df, filter_results$literature_info$by_study)
    export_df$N_Studies <- filter_results$literature_info$n_studies
  }
  
  # Save detailed per-CpG results
  write.csv(export_df, 
            file.path(output_dir, "cpg_filtering_criteria.csv"),
            row.names = FALSE)
  
  # Save list of selected CpG names only
  write.table(filter_results$selected_cpgs,
              file.path(output_dir, "selected_cpgs.txt"),
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  # Save summary statistics
  summary_df <- data.frame(
    Metric = names(unlist(filter_results$filters_summary)),
    Value = unlist(filter_results$filters_summary)
  )
  write.csv(summary_df,
            file.path(output_dir, "filtering_summary.csv"),
            row.names = FALSE)
  
  message(paste("\nResults exported to:", output_dir))
}

# ============================================
# STEP 2: NESTED CROSS-VALIDATION SETUP
# ============================================

#' Create nested cross-validation fold structure
#' 
#' Generates indices for nested cross-validation, where outer folds
#' are used for unbiased performance estimation and inner folds for hyperparameter
#' tuning. Uses stratified sampling for classification problems.
#' 
#' @param y Response variable (factor for classification, numeric for regression)
#' @param k_outer Number of outer folds (default 5)
#' @param k_inner Number of inner folds (default 5)
#' @param seed Random seed for reproducibility
#' @return List of nested fold structures with outer_test, outer_train, and inner_folds
create_nested_folds <- function(y, k_outer = 5, k_inner = 5, seed = 123) {
  set.seed(seed)
  
  # Create stratified folds for classification to maintain class balance
  # Use regular k-fold for regression
  if(is.factor(y)) {
    outer_folds <- createFolds(y, k = k_outer, list = TRUE, returnTrain = FALSE)
  } else {
    outer_folds <- createFolds(seq_along(y), k = k_outer, list = TRUE, returnTrain = FALSE)
  }
  
  # For each outer fold, create inner folds from the training subset
  nested_folds <- list()
  for(i in 1:k_outer) {
    train_indices <- setdiff(seq_along(y), outer_folds[[i]])
    y_train <- y[train_indices]
    
    # Create inner folds (stratified for classification)
    if(is.factor(y_train)) {
      inner_folds <- createFolds(y_train, k = k_inner, list = TRUE, returnTrain = FALSE)
    } else {
      inner_folds <- createFolds(seq_along(y_train), k = k_inner, 
                                 list = TRUE, returnTrain = FALSE)
    }
    
    # Convert inner fold indices back to original data indices
    inner_folds <- lapply(inner_folds, function(idx) train_indices[idx])
    
    nested_folds[[i]] <- list(
      outer_test = outer_folds[[i]],
      outer_train = train_indices,
      inner_folds = inner_folds
    )
  }
  
  return(nested_folds)
}

# ============================================
# STEP 3: ELASTIC NET WITH NESTED CV
# ============================================

#' Convert nested fold structure to glmnet foldid format
#' 
#' @param inner_folds List of inner fold indices (in original data coordinates)
#' @param train_idx Vector of training indices for the current outer fold
#' @return Integer vector where each element indicates which inner fold that sample belongs to
create_foldid <- function(inner_folds, train_idx) {
  n_total <- length(train_idx)
  foldid <- rep(NA_integer_, n_total)
  
  for(i in seq_along(inner_folds)) {
    idx_orig  <- inner_folds[[i]]            # Indices in original data coordinates
    idx_train <- match(idx_orig, train_idx)  # Convert to training set coordinates (1:n_train)
    foldid[idx_train] <- i
  }
  return(foldid)
}

#' Elastic Net regression/classification with nested cross-validation
#' 
#' Performs elastic net with proper nested CV to avoid overfitting during
#' hyperparameter selection. Inner CV tunes alpha and lambda; outer CV
#' provides unbiased performance estimates. Supports optional covariates
#' that are included but not penalized.
#' 
#' @param X CpG beta matrix (CpGs in rows, samples in columns)
#' @param y Response variable (factor for classification, numeric for regression)
#' @param covariates Optional data frame of covariates (samples in rows)
#' @param nested_folds Nested fold structure from create_nested_folds()
#' @param alpha_range Vector of alpha values to search (0=ridge, 1=lasso)
#' @return List with performance metrics, selected CpGs, coefficients, and final model
elastic_net_nested_cv <- function(X, y, 
                                  covariates = NULL,
                                  nested_folds, 
                                  alpha_range = seq(0, 1, 0.1)) {
  
  # Validate input dimensions
  if(ncol(X) != length(y)) {
    stop("Number of samples in X (", ncol(X), 
         ") doesn't match length of y (", length(y), ")")
  }
  
  if(!is.null(covariates)) {
    if(nrow(covariates) != ncol(X)) {
      stop("Number of rows in covariates must match number of samples in X")
    }
  }
  
  n_outer <- length(nested_folds)
  outer_results <- list()
  all_selected_cpgs <- list()
  
  # Combine CpG data with covariates if provided
  # Covariates are appended to feature matrix but will not be penalized
  if(!is.null(covariates)) {
    message("Including unpenalized covariates in elastic net model")
    X_combined <- rbind(X, t(covariates))
    
    # Penalty factors: 1 = penalized (CpGs), 0 = unpenalized (covariates)
    penalty_factors <- c(rep(1, nrow(X)), rep(0, ncol(covariates)))
    cpg_indices <- 1:nrow(X)
    covariate_names <- colnames(covariates)
  } else {
    X_combined <- X
    penalty_factors <- rep(1, nrow(X))
    cpg_indices <- 1:nrow(X)
    covariate_names <- NULL
  }
  
  # Determine evaluation metric based on outcome type
  metric <- ifelse(is.factor(y), "Accuracy", "RMSE")
  
  # Create progress bar for outer folds
  pb_outer <- progress_bar$new(
    format = "  Elastic Net outer folds [:bar] :percent | Fold :current/:total | ETA: :eta",
    total = n_outer,
    clear = FALSE,
    width = 80
  )
  
  for(i in 1:n_outer) {
    message(paste("\nOuter fold", i, "of", n_outer))
    
    # Split data into outer training and test sets
    train_idx <- nested_folds[[i]]$outer_train
    test_idx <- nested_folds[[i]]$outer_test
    inner_folds <- nested_folds[[i]]$inner_folds
    
    X_train <- X_combined[, train_idx, drop = FALSE]
    X_test <- X_combined[, test_idx, drop = FALSE]
    y_train <- y[train_idx]
    y_test <- y[test_idx]
    
    # ----------------------------------------
    # INNER CV: Grid search over alpha values
    # Lambda is selected automatically by cv.glmnet
    # ----------------------------------------
    best_alpha <- NA
    best_lambda <- NA
    best_score <- ifelse(metric == "RMSE", Inf, Inf)  # Lower is better for both deviance and RMSE
    
    # Progress bar for alpha grid search
    pb_alpha <- progress_bar$new(
      format = "    Tuning alpha [:bar] :percent | Alpha :current/:total",
      total = length(alpha_range),
      clear = FALSE,
      width = 70
    )
    
    for(alpha in alpha_range) {
      # Run cv.glmnet with inner folds for lambda selection
      # parallel = TRUE enables parallel computation of the CV folds
      cv_model <- cv.glmnet(
        x = t(X_train),
        y = y_train,
        alpha = alpha,
        penalty.factor = penalty_factors,
        nfolds = length(inner_folds),
        foldid = create_foldid(inner_folds, train_idx),
        family = ifelse(is.factor(y_train), "binomial", "gaussian"),
        parallel = TRUE  # Enable parallel processing for inner CV
      )
      
      # Get cross-validated performance at lambda.min
      # cvm is cross-validated mean deviance (lower is better)
      score <- min(cv_model$cvm)
      if(score < best_score) {
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
      x = t(X_train),
      y = y_train,
      alpha = best_alpha,
      lambda = best_lambda,
      penalty.factor = penalty_factors,
      family = ifelse(is.factor(y_train), "binomial", "gaussian")
    )
    
    # Extract coefficients (excluding intercept)
    coefs <- as.matrix(coef(final_model, s = best_lambda))
    
    # Separate CpG coefficients from covariate coefficients
    # Row 1 is intercept, rows 2:(nrow(X)+1) are CpGs
    cpg_coefs <- coefs[2:(nrow(X)+1), 1]
    selected_features <- rownames(X)[which(abs(cpg_coefs) > 0)]
    all_selected_cpgs[[i]] <- selected_features
    
    # Extract covariate coefficients if present
    if(!is.null(covariates)) {
      covariate_coefs <- coefs[(nrow(X)+2):nrow(coefs), 1]
      names(covariate_coefs) <- covariate_names
    } else {
      covariate_coefs <- NULL
    }
    
    # ----------------------------------------
    # Evaluate on outer test set (unbiased estimate)
    # ----------------------------------------
    predictions <- predict(final_model, newx = t(X_test), 
                           s = best_lambda, type = "response")
    
    # Calculate performance metrics based on outcome type
    if(is.factor(y)) {
      # Classification metrics
      pred_classes <- ifelse(predictions > 0.5, levels(y)[2], levels(y)[1])
      pred_classes <- factor(pred_classes, levels = levels(y))
      
      accuracy <- mean(pred_classes == y_test)
      
      # Calculate precision, recall, and F1 score
      confusion <- table(Predicted = pred_classes, Actual = y_test)
      if(nrow(confusion) == 2 && ncol(confusion) == 2) {
        precision <- confusion[2,2] / sum(confusion[2,])
        recall <- confusion[2,2] / sum(confusion[,2])
        f1_score <- ifelse(is.na(precision) || is.na(recall) || 
                             (precision + recall) == 0, 0,
                           2 * (precision * recall) / (precision + recall))
      } else {
        f1_score <- NA
      }
      
      # Calculate AUC for binary classification
      if(nlevels(y) == 2) {
        auc_value <- auc(roc(y_test, predictions[,1], quiet = TRUE))
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
      r_squared <- 1 - sum((y_test - predictions)^2) / sum((y_test - mean(y_test))^2)
      
      # Adjusted R-squared accounts for number of predictors
      n <- length(y_test)
      p <- length(selected_features) + ifelse(!is.null(covariates), ncol(covariates), 0)
      # adj_r_squared <- 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))
      
      outer_results[[i]] <- list(
        rmse = rmse,
        mae = mae,
        r_squared = r_squared,
        # adj_r_squared = adj_r_squared,
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
  if(is.factor(y)) {
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
    # avg_adj_r2 <- mean(sapply(outer_results, function(x) x$adj_r_squared))
    
    final_metrics <- data.frame(
      Model = "Elastic Net",
      RMSE_Mean = avg_rmse,
      RMSE_SD = sd(sapply(outer_results, function(x) x$rmse)),
      R2_Mean = avg_r2,
      R2_SD = sd(sapply(outer_results, function(x) x$r_squared)),
      # Adj_R2_Mean = avg_adj_r2,
      # Adj_R2_SD = sd(sapply(outer_results, function(x) x$adj_r_squared)),
      Avg_Features = mean(sapply(outer_results, function(x) x$n_features))
    )
  }
  
  # Average covariate coefficients across folds (for interpretation)
  if (!is.null(covariates)) {
    cov_mat <- sapply(outer_results, function(z) as.numeric(z$covariate_coefficients))
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
    x = t(X_combined),
    y = y,
    alpha = mean(sapply(outer_results, function(x) x$best_alpha)),
    lambda = mean(sapply(outer_results, function(x) x$best_lambda)),
    penalty.factor = penalty_factors,
    family = ifelse(is.factor(y), "binomial", "gaussian")
  )
  
  final_coefs <- coef(final_model)
  final_coefs_matrix <- as.matrix(final_coefs)
  
  return(list(
    performance = final_metrics,
    outer_fold_results = outer_results,
    selected_cpgs = names(sort(cpg_frequency, decreasing = TRUE)),
    cpg_selection_frequency = sort(cpg_frequency, decreasing = TRUE),
    stable_cpgs = names(cpg_frequency[cpg_frequency >= n_outer * 0.8]),  # Selected in 80%+ of folds
    final_model = final_model,
    final_coefs_matrix = final_coefs_matrix,
    coefficients = final_coefs_matrix[2:(nrow(X)+1), 1],  # CpG coefficients only (excluding intercept)
    covariate_coefficients = avg_covariate_coefs,
    penalty_factors = penalty_factors
  ))
}

# ============================================
# STEP 4: RANDOM FOREST WITH NESTED CV
# ============================================

#' Random Forest with nested cross-validation
#' 
#' Performs random forest with nested CV for unbiased performance estimation.
#' Supports both ranger (faster, permutation importance) and caret RF.
#' Inner CV tunes mtry parameter; outer CV evaluates performance.
#' 
#' @param X CpG beta matrix (CpGs in rows, samples in columns)
#' @param y Response variable (factor for classification, numeric for regression)
#' @param nested_folds Nested fold structure from create_nested_folds()
#' @param covariates Optional data frame of covariates (samples in rows)
#' @param ntree Number of trees in the forest (default 500)
#' @param hyperparam_grid Data frame for tuning hyperparameters: mtry, min.node.size, and sample fraction
#' @param importance_threshold Number of top CpGs to retain based on importance
#' @return List with performance metrics, importance scores, and final model
random_forest_nested_cv_tuned <- function(X, y, nested_folds, 
                                          covariates = NULL,
                                          ntree = 500,
                                          hyperparam_grid = NULL,
                                          importance_threshold = 100) {
  
  # Default hyperparameter grid if not provided
  if(is.null(hyperparam_grid)) {
    hyperparam_grid <- expand.grid(
      mtry = c(250, 500, 750, 1000),  # Could rewrite this to be adaptive to p/4, p/3, p/2
      min.node.size = c(5, 10, 20), 
      sample.fraction = c(0.8, 1.0) 
    )
  }
  
  cat("Testing", nrow(hyperparam_grid), "hyperparameter combinations in inner CV\n")
  
  n_outer <- length(nested_folds)
  outer_results <- list()
  all_importance_scores <- list()
  all_best_params <- list()
  
  # Combine CpG data with covariates if provided
  if(!is.null(covariates)) {
    message("Including covariates in Random Forest model")
    X_combined <- rbind(X, t(covariates))
    feature_names <- c(rownames(X), colnames(covariates))
    cpg_names <- rownames(X)
  } else {
    X_combined <- X
    feature_names <- rownames(X)
    cpg_names <- rownames(X)
  }
  
  metric <- ifelse(is.factor(y), "Accuracy", "RMSE")
  
  # Progress bar for outer folds
  pb_outer <- progress_bar$new(
    format = "  RF outer folds [:bar] :percent | Fold :current/:total | ETA: :eta",
    total = n_outer,
    clear = FALSE,
    width = 80
  )
  
  for(i in 1:n_outer) {
    message(paste("\n=== Outer fold", i, "of", n_outer, "==="))
    
    # Split data into outer training and test sets
    train_idx <- nested_folds[[i]]$outer_train
    test_idx <- nested_folds[[i]]$outer_test
    inner_folds <- nested_folds[[i]]$inner_folds
    
    X_train <- X_combined[, train_idx, drop = FALSE]
    X_test <- X_combined[, test_idx, drop = FALSE]
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
    
    for(j in 1:nrow(hyperparam_grid)) {
      inner_scores <- numeric(length(inner_folds))
      
      # Evaluate this hyperparameter combination on inner folds
      for(k in 1:length(inner_folds)) {
        # Split inner training/validation
        inner_train_idx <- setdiff(train_idx, inner_folds[[k]])
        inner_val_idx <- inner_folds[[k]]
        
        # Train model with these hyperparameters
        model_temp <- ranger(
          x = t(X_combined[, inner_train_idx, drop = FALSE]),
          y = y[inner_train_idx],
          num.trees = ntree,
          mtry = hyperparam_grid$mtry[j],
          min.node.size = hyperparam_grid$min.node.size[j],
          sample.fraction = hyperparam_grid$sample.fraction[j],
          num.threads = 1,  # Don't nest parallelization
          classification = is.factor(y),
          verbose = FALSE
        )
        
        # Predict on inner validation set
        pred_temp <- predict(model_temp, 
                             data = t(X_combined[, inner_val_idx, drop = FALSE]))
        
        # Calculate validation score (lower is better)
        if(is.factor(y)) {
          inner_scores[k] <- 1 - mean(pred_temp$predictions == y[inner_val_idx])
        } else {
          inner_scores[k] <- sqrt(mean((pred_temp$predictions - y[inner_val_idx])^2))
        }
      }
      
      # Average score across inner folds
      avg_score <- mean(inner_scores)
      
      # Track best hyperparameters
      if(avg_score < best_score) {
        best_score <- avg_score
        best_params <- hyperparam_grid[j, ]
      }
      
      pb_hyper$tick()
    }
    
    message(paste("\n  Best params: mtry =", best_params$mtry,
                  ", min.node.size =", best_params$min.node.size,
                  ", sample.fraction =", best_params$sample.fraction))
    message(paste("  Inner CV score:", round(best_score, 4)))
    
    all_best_params[[i]] <- best_params
    
    # ----------------------------------------
    # Train final model for this outer fold with best hyperparameters
    # ----------------------------------------
    ranger_model <- ranger(
      x = t(X_train),
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
    pred_ranger <- predict(ranger_model, data = t(X_test))
    
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
    n_features_used = nrow(X_combined)
    
    # ----------------------------------------
    # Calculate performance metrics
    # ----------------------------------------
    if (is.factor(y)) {
      accuracy <- mean(predictions == y_test)
      
      confusion <- table(Predicted = predictions, Actual = y_test)
      if (nrow(confusion) == 2 && ncol(confusion) == 2) {
        precision <- confusion[2, 2] / sum(confusion[2, ])
        recall    <- confusion[2, 2] / sum(confusion[, 2])
        f1_score  <- ifelse(is.na(precision) || is.na(recall) || 
                              (precision + recall) == 0, 0,
                            2 * (precision * recall) / (precision + recall))
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
      mae  <- mean(abs(predictions - y_test))
      r_squared <- 1 - sum((y_test - predictions)^2) / sum((y_test - mean(y_test))^2)
      
      n <- length(y_test)
      p <- n_features_used
      # adj_r_squared <- 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))
      
      outer_results[[i]] <- list(
        rmse = rmse,
        mae = mae,
        r_squared = r_squared,
        # adj_r_squared = adj_r_squared,
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
  if(is.factor(y)) {
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
      Avg_Features = mean(sapply(outer_results, function(x) x$n_features_selected))
    )
  } else {
    avg_rmse <- mean(sapply(outer_results, function(x) x$rmse))
    avg_r2 <- mean(sapply(outer_results, function(x) x$r_squared))
    # avg_adj_r2 <- mean(sapply(outer_results, function(x) x$adj_r_squared))
    
    final_metrics <- data.frame(
      Model = "Random Forest (Tuned)",
      RMSE_Mean = avg_rmse,
      RMSE_SD = sd(sapply(outer_results, function(x) x$rmse)),
      R2_Mean = avg_r2,
      R2_SD = sd(sapply(outer_results, function(x) x$r_squared)),
      # Adj_R2_Mean = avg_adj_r2,
      # Adj_R2_SD = sd(sapply(outer_results, function(x) x$adj_r_squared)),
      Avg_Features = mean(sapply(outer_results, function(x) x$n_features_selected))
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
  message(paste("mtry:", best_mtry, 
                "- chosen in", sum(hyperparam_df$mtry == best_mtry), "out of", n_outer, "folds"))
  message(paste("min.node.size:", best_min_node_size,
                "- chosen in", sum(hyperparam_df$min.node.size == best_min_node_size), "out of", n_outer, "folds"))
  message(paste("sample.fraction:", best_sample_fraction,
                "- chosen in", sum(hyperparam_df$sample.fraction == best_sample_fraction), "out of", n_outer, "folds"))
  
  # ============================================
  # TRAIN FINAL MODEL ON ALL DATA
  # Using MOST FREQUENT hyperparameters
  # ============================================
  message("\nTraining final model on all data with best hyperparameters...")
  final_model <- ranger(
    x = t(X_combined),
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
  
  return(list(
    performance = final_metrics,
    outer_fold_results = outer_results,
    avg_importance = sort(avg_importance, decreasing = TRUE),
    top_cpgs = names(sort(avg_importance, decreasing = TRUE)[1:min(100, length(avg_importance))]),
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
  ))
}

# ============================================
# STEP 5: CREATE METHYLATION PROFILE SCORE
# ============================================

#' Create methylation profile score (MPS) from model results
#' 
#' Calculates a weighted or unweighted sum of methylation values at selected
#' CpG sites to create a summary score representing PM2.5 exposure signature.
#' 
#' @param beta_matrix Full beta matrix (CpGs in rows, samples in columns)
#' @param model_results Output from elastic_net_nested_cv() or random_forest_nested_cv()
#' @param method "elastic_net" or "random_forest"
#' @param score_type "weighted" (by coefficients/importance) or "unweighted" (mean)
#' @return Numeric vector of standardized profile scores (one per sample)
create_methylation_profile_score <- function(beta_matrix, model_results, 
                                          method = c("elastic_net", "random_forest"),
                                          score_type = c("weighted", "unweighted")) {
  method <- match.arg(method)
  score_type <- match.arg(score_type)
  
  if(method == "elastic_net") {
    # Extract coefficients from elastic net model
    coefs <- model_results$coefficients
    if(is.matrix(coefs)) {
      coefs <- as.vector(coefs)
      names(coefs) <- rownames(model_results$coefficients)
    }
    
    # Get CpGs with non-zero coefficients (selected by elastic net)
    selected_cpgs <- names(coefs[abs(coefs) > 0])
    
    if(length(selected_cpgs) == 0) {
      warning("No CpGs selected by elastic net")
      return(rep(0, ncol(beta_matrix)))
    }
    
    # Subset beta matrix to selected CpGs
    beta_subset <- beta_matrix[selected_cpgs, , drop = FALSE]
    
    if(score_type == "weighted") {
      # Weighted sum using elastic net coefficients
      weights <- coefs[selected_cpgs]
      risk_score <- colSums(beta_subset * weights, na.rm = TRUE)
    } else {
      # Simple unweighted mean across selected CpGs
      risk_score <- colMeans(beta_subset, na.rm = TRUE)
    }
    
  } else if(method == "random_forest") {
    # Use top important CpGs from random forest
    selected_cpgs <- model_results$top_cpgs
    
    # Ensure selected CpGs exist in beta matrix
    selected_cpgs <- intersect(selected_cpgs, rownames(beta_matrix))
    
    if(length(selected_cpgs) == 0) {
      warning("No selected CpGs found in beta matrix")
      return(rep(0, ncol(beta_matrix)))
    }
    
    beta_subset <- beta_matrix[selected_cpgs, , drop = FALSE]
    
    if(score_type == "weighted") {
      # Weight by variable importance scores
      weights <- model_results$avg_importance[selected_cpgs]
      w_sum <- sum(weights)
      if (w_sum == 0) {
        warning("Sum of importance weights is zero; using unweighted mean instead.")
        risk_score <- colMeans(beta_subset, na.rm = TRUE)
      } else {
        weights <- weights / w_sum  # Normalize weights to sum to 1
        risk_score <- colSums(beta_subset * weights, na.rm = TRUE)
      }
    } else {
      risk_score <- colMeans(beta_subset, na.rm = TRUE)
    }
  }
  
  # Z-score standardize risk score to z-scores (mean=0, sd=1)
  risk_score_z <- scale(risk_score)[, 1]
  
  return(data.frame(
    sample_id        = names(risk_score),
    risk_score       = risk_score,
    risk_score_z     = risk_score_z,
    row.names        = NULL,
    stringsAsFactors = FALSE
  ))
}

# ============================================
# STEP 6: COMPREHENSIVE MODEL COMPARISON
# ============================================

#' Compare performance across all models
#' 
#' @param results_continuous List with elastic_net and random_forest results for continuous outcome
#' @param results_dichotomous List with elastic_net and random_forest results for dichotomous outcome
#' @return List with comparison table, best models, and recommendations
compare_all_models <- function(results_continuous, results_dichotomous) {
  
  # Combine performance metrics into single table
  comparison <- rbind(
    cbind(Outcome = "Continuous", results_continuous$elastic_net$performance),
    cbind(Outcome = "Continuous", results_continuous$random_forest$performance),
    cbind(Outcome = "Dichotomous", results_dichotomous$elastic_net$performance),
    cbind(Outcome = "Dichotomous", results_dichotomous$random_forest$performance)
  )
  
  # Calculate composite scores that balance performance and stability
  # Higher R2 with lower SD is preferred for continuous
  if("R2_Mean" %in% names(comparison)) {
    comparison$Score_Continuous <- comparison$R2_Mean - 
      0.5 * comparison$R2_SD
  }
  
  # Higher AUC with lower SD is preferred for dichotomous
  if("AUC_Mean" %in% names(comparison)) {
    comparison$Score_Dichotomous <- comparison$AUC_Mean - 
      0.5 * comparison$AUC_SD
  }
  
  # Identify best model for each outcome type
  cont_models <- comparison[comparison$Outcome == "Continuous", ]
  dicho_models <- comparison[comparison$Outcome == "Dichotomous", ]
  
  best_continuous <- cont_models[which.max(cont_models$Score_Continuous), ]
  best_dichotomous <- dicho_models[which.max(dicho_models$Score_Dichotomous), ]
  
  return(list(
    comparison_table = comparison,
    best_continuous_model = best_continuous,
    best_dichotomous_model = best_dichotomous,
    recommendations = paste(
      "For continuous PM2.5:", best_continuous$Model,
      "\nFor dichotomous PM2.5:", best_dichotomous$Model
    )
  ))
}

#' Compare CpG filtering results between continuous and dichotomous outcomes
#' 
#' @param filter_cont Filtering results for continuous outcome
#' @param filter_dicho Filtering results for dichotomous outcome
#' @return List with overlap statistics and shared CpGs
compare_filtering_results <- function(filter_cont, filter_dicho) {
  
  # Calculate overlap between the two filtering approaches
  overlap_cpgs <- intersect(filter_cont$selected_cpgs, filter_dicho$selected_cpgs)
  only_cont <- setdiff(filter_cont$selected_cpgs, filter_dicho$selected_cpgs)
  only_dicho <- setdiff(filter_dicho$selected_cpgs, filter_cont$selected_cpgs)
  
  # Calculate overlap percentage relative to smaller set
  min_sel <- min(filter_cont$n_selected, filter_dicho$n_selected)
  pct_overlap <- if (min_sel > 0) {
    100 * length(overlap_cpgs) / min_sel
  } else {
    NA_real_
  }
  
  comparison <- list(
    n_continuous = filter_cont$n_selected,
    n_dichotomous = filter_dicho$n_selected,
    n_overlap = length(overlap_cpgs),
    pct_overlap = pct_overlap,
    only_continuous = length(only_cont),
    only_dichotomous = length(only_dicho),
    overlap_cpgs = overlap_cpgs
  )
  
  message("\n--- Filtering Comparison (Continuous vs Dichotomous) ---")
  message(paste("CpGs selected for continuous:", comparison$n_continuous))
  message(paste("CpGs selected for dichotomous:", comparison$n_dichotomous))
  message(paste("Overlapping CpGs:", comparison$n_overlap, 
                "(", ifelse(is.na(comparison$pct_overlap),
                            "NA", round(comparison$pct_overlap, 1)), "% overlap)"))
  
  return(comparison)
}

# ============================================
# MAIN PIPELINE WRAPPER
# ============================================

#' Run complete methylation analysis pipeline
#' 
#' Executes the full analysis workflow: CpG filtering, nested CV with elastic net
#' and random forest, risk score creation, and model comparison for both
#' continuous and dichotomous PM2.5 outcomes.
#' 
#' @param beta_matrix Numeric matrix of beta values (CpGs in rows, samples in columns)
#' @param pm25_continuous Numeric vector of continuous PM2.5 exposure values
#' @param pm25_dichotomous Factor with two levels for dichotomous exposure
#' @param covariates Optional data frame of covariates (samples in rows)
#' @param literature_studies Optional named list of CpG vectors from published studies
#' @param k_outer Number of outer CV folds (default 5)
#' @param k_inner Number of inner CV folds (default 5)
#' @param var_percentile Variance filter threshold (default 80 = top 20%)
#' @param pval_threshold P-value threshold for association filter (default 0.05)
#' @param min_criteria Minimum criteria for CpG selection (default 2)
#' @param alpha_range Alpha values to search for elastic net (default 0 to 1 by 0.05)
#' @param seed Random seed for reproducibility
#' @return List containing all results for both outcomes and model comparisons
run_complete_pipeline <- function(beta_matrix, 
                                  pm25_continuous, 
                                  pm25_dichotomous,
                                  covariates = NULL,
                                  literature_studies = NULL,
                                  k_outer = 5,
                                  k_inner = 5,
                                  var_percentile = 80,
                                  pval_threshold = 0.05,
                                  min_criteria = 2,
                                  alpha_range = seq(0, 1, 0.05),
                                  seed = 123) {
  
  set.seed(seed)
  results <- list()
  
  # Validate covariate dimensions
  if(!is.null(covariates)) {
    if(nrow(covariates) != ncol(beta_matrix)) {
      stop("Number of rows in covariates must match number of samples (columns) in beta_matrix")
    }
    message(paste("Including", ncol(covariates), "covariates in models"))
    message(paste("Covariates:", paste(colnames(covariates), collapse = ", ")))
  }
  
  # Ensure dichotomous outcome is a factor
  if(!is.null(pm25_dichotomous) && !is.factor(pm25_dichotomous)) {
    pm25_dichotomous <- as.factor(pm25_dichotomous)
  }
  
  # ============================================
  # PROCESS CONTINUOUS OUTCOME
  # ============================================
  message("\n=== PROCESSING CONTINUOUS OUTCOME ===")
  
  # Filter CpGs using multi-criteria approach
  filtered_contin <- filter_cpg_sites(beta_matrix, pm25_continuous,
                                      literature_studies = literature_studies,
                                      var_percentile = var_percentile,
                                      pval_threshold = pval_threshold,
                                      min_criteria = min_criteria)
  
  results$continuous_filtering <- filtered_contin
  
  # Create nested CV fold structure
  nested_folds_contin <- create_nested_folds(pm25_continuous, 
                                             k_outer = k_outer, 
                                             k_inner = k_inner,
                                             seed = seed)
  
  # Run elastic net with nested CV
  message("\nRunning Elastic Net with nested CV...")
  en_results_contin <- elastic_net_nested_cv(X = filtered_contin$filtered_beta, 
                                             y = pm25_continuous,
                                             covariates = covariates,
                                             nested_folds = nested_folds_contin,
                                             alpha_range = alpha_range)
  
  # Run random forest with nested CV
  message("\nRunning Random Forest with nested CV...")
  rf_results_contin <- random_forest_nested_cv(X = filtered_contin$filtered_beta, 
                                               y = pm25_continuous,
                                               covariates = covariates,
                                               nested_folds = nested_folds_contin,
                                               importance_method = "all")
  
  # Create methylation risk scores
  en_score_contin <- create_methylation_risk_score(beta_matrix, en_results_contin,
                                                   method = "elastic_net",
                                                   score_type = "weighted")
  
  rf_score_contin <- create_methylation_risk_score(beta_matrix, rf_results_contin,
                                                   method = "random_forest",
                                                   score_type = "weighted")
  
  results$continuous <- list(
    filtered_cpgs = filtered_contin$selected_cpgs,
    n_cpgs_filtered = filtered_contin$n_selected,
    filtering_details = filtered_contin,
    elastic_net = en_results_contin,
    random_forest = rf_results_contin,
    elastic_net_score = en_score_contin,
    random_forest_score = rf_score_contin
  )
  
  # ============================================
  # PROCESS DICHOTOMOUS OUTCOME
  # ============================================
  message("\n=== PROCESSING DICHOTOMOUS OUTCOME ===")
  
  # Filter CpGs using multi-criteria approach
  filtered_dichot <- filter_cpg_sites(beta_matrix, pm25_dichotomous,
                                      literature_studies = literature_studies,
                                      var_percentile = var_percentile,
                                      pval_threshold = pval_threshold,
                                      min_criteria = min_criteria)
  
  results$dichotomous_filtering <- filtered_dichot
  
  # Create nested CV fold structure (stratified for classification)
  nested_folds_dichot <- create_nested_folds(pm25_dichotomous, 
                                             k_outer = k_outer, 
                                             k_inner = k_inner,
                                             seed = seed)
  
  # Run elastic net with nested CV
  message("\nRunning Elastic Net with nested CV...")
  en_results_dichot <- elastic_net_nested_cv(X = filtered_dichot$filtered_beta, 
                                             y = pm25_dichotomous,
                                             covariates = covariates,
                                             nested_folds = nested_folds_dichot,
                                             alpha_range = alpha_range)
  
  # Run random forest with nested CV
  message("\nRunning Random Forest with nested CV...")
  rf_results_dichot <- random_forest_nested_cv(X = filtered_dichot$filtered_beta, 
                                               y = pm25_dichotomous,
                                               covariates = covariates,
                                               nested_folds = nested_folds_dichot,
                                               importance_method = "all")
  
  # Create methylation risk scores
  en_score_dichot <- create_methylation_risk_score(beta_matrix, en_results_dichot,
                                                   method = "elastic_net",
                                                   score_type = "weighted")
  
  rf_score_dichot <- create_methylation_risk_score(beta_matrix, rf_results_dichot,
                                                   method = "random_forest",
                                                   score_type = "weighted")
  
  results$dichotomous <- list(
    filtered_cpgs = filtered_dichot$selected_cpgs,
    n_cpgs_filtered = filtered_dichot$n_selected,
    filtering_details = filtered_dichot,
    elastic_net = en_results_dichot,
    random_forest = rf_results_dichot,
    elastic_net_score = en_score_dichot,
    random_forest_score = rf_score_dichot
  )
  
  # ============================================
  # COMPARE ALL MODELS
  # ============================================
  message("\n=== COMPARING ALL MODELS ===")
  
  results$comparison <- compare_all_models(
    list(elastic_net = results$continuous$elastic_net,
         random_forest = results$continuous$random_forest),
    list(elastic_net = results$dichotomous$elastic_net,
         random_forest = results$dichotomous$random_forest)
  )
  
  # Compare filtering between continuous and dichotomous approaches
  results$filtering_comparison <- compare_filtering_results(
    filtered_contin, filtered_dichot
  )
  
  # Print final summary
  message("\n=== PIPELINE COMPLETE ===")
  message("\nModel Performance Summary:")
  print(results$comparison$comparison_table)
  message(paste("\n", results$comparison$recommendations))
  
  return(results)
}

# ============================================
# VISUALIZATION AND ANALYSIS FUNCTIONS
# ============================================

#' Create Venn diagram showing CpG filtering criteria overlap
#' 
#' @param filter_results Output from filter_cpg_sites()
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
sensitivity_analysis_filtering <- function(beta_matrix, pm25_data, 
                                           literature_studies,
                                           var_percentiles = c(70, 75, 80, 85, 90),
                                           pval_thresholds = c(0.01, 0.05, 0.1),
                                           min_criteria_options = c(1, 2, 3)) {
  
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
    format = "  Sensitivity analysis [:bar] :percent | :current/:total combinations",
    total = nrow(results_grid),
    clear = FALSE,
    width = 80
  )
  
  for(i in 1:nrow(results_grid)) {
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
    if(!is.null(filter_res$literature_info)) {
      lit_in_selected <- sum(filter_res$selected_cpgs %in% 
                               filter_res$literature_info$literature_cpgs)
      if (filter_res$n_selected > 0) {
        results_grid$pct_literature[i] <- 100 * lit_in_selected / filter_res$n_selected
      } else {
        results_grid$pct_literature[i] <- NA_real_
      }
    }
    
    pb$tick()
  }
  
  return(results_grid)
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
  
  for(period in names(pm25_data_list)) {
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
  overlap_matrix <- matrix(0, nrow = length(periods), ncol = length(periods),
                           dimnames = list(periods, periods))
  
  for(i in 1:length(periods)) {
    for(j in 1:length(periods)) {
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
    N_Variance = sapply(temporal_results, function(x) 
      sum(x$criteria_matrix[, "Variance"])),
    N_Association = sapply(temporal_results, function(x) 
      sum(x$criteria_matrix[, "Association"])),
    N_Literature = sapply(temporal_results, function(x) 
      sum(x$criteria_matrix[, "Literature"]))
  )
  
  return(list(
    by_period = temporal_results,
    consistent_cpgs = consistent_cpgs,
    n_consistent = length(consistent_cpgs),
    overlap_matrix = overlap_matrix,
    summary = summary_df
  ))
}

#' Detailed analysis of literature CpGs in filtering results
#' 
#' @param filter_results Output from filter_cpg_sites()
#' @param beta_matrix Numeric matrix of beta values
#' @param pm25_data PM2.5 exposure data
#' @return List with per-CpG details and summary statistics
analyze_literature_cpgs <- function(filter_results, beta_matrix, pm25_data) {
  
  if(is.null(filter_results$literature_info)) {
    message("No literature studies provided")
    return(NULL)
  }
  
  lit_info <- filter_results$literature_info
  
  # Group CpGs by number of studies they appear in (replication count)
  cpgs_by_n_studies <- list()
  for(n in 1:4) {
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
    Passes_Association = filter_results$criteria_matrix[lit_cpgs, "Association"],
    N_Studies = lit_info$n_studies[lit_cpgs],
    Selected = as.integer(lit_cpgs %in% filter_results$selected_cpgs)
  )
  
  # Calculate effect sizes (correlation) for literature CpGs
  if(is.numeric(pm25_data)) {
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
  
  return(list(
    cpgs_by_study_count = cpgs_by_n_studies,
    literature_cpg_details = lit_criteria_summary,
    summary = summary_stats
  ))
}

# ============================================
# TRIMESTER ANALYSIS FUNCTIONS
# ============================================

#' Run complete pipeline across all pregnancy trimesters
#' 
#' @param beta_matrix Numeric matrix of beta values
#' @param pm25_data Data frame with trimester-specific PM2.5 columns
#' @param trimesters Vector of trimester names to analyze
#' @param k_outer Number of outer CV folds
#' @param k_inner Number of inner CV folds
#' @param var_percentile Variance filter threshold
#' @param pval_threshold P-value threshold
#' @param seed Random seed base (incremented for each trimester)
#' @return List with per-trimester results and cross-trimester comparison
run_trimester_analysis <- function(beta_matrix, pm25_data, 
                                   trimesters = c("trim1", "trim2", "trim3", "overall"),
                                   k_outer = 5, k_inner = 5,
                                   var_percentile = 80,
                                   pval_threshold = 0.05,
                                   seed = 123) {
  
  trimester_results <- list()
  
  for(trim in trimesters) {
    message(paste("\n\n===== ANALYZING", toupper(trim), "====="))
    
    # Extract continuous PM2.5 for this trimester and create dichotomous version
    pm25_cont <- pm25_data[[paste0(trim, "_avg")]]
    pm25_median <- median(pm25_cont, na.rm = TRUE)
    pm25_dicho <- factor(ifelse(pm25_cont >= pm25_median, "High", "Low"))
    
    # Run full pipeline with different seed for each trimester
    trimester_results[[trim]] <- run_complete_pipeline(
      beta_matrix = beta_matrix,
      pm25_continuous = pm25_cont,
      pm25_dichotomous = pm25_dicho,
      k_outer = k_outer,
      k_inner = k_inner,
      var_percentile = var_percentile,
      pval_threshold = pval_threshold,
      seed = seed + which(trimesters == trim)
    )
  }
  
  # Compare results across trimesters
  trimester_comparison <- compare_trimester_results(trimester_results)
  
  return(list(
    by_trimester = trimester_results,
    trimester_comparison = trimester_comparison
  ))
}

#' Compare model performance and CpG selection across trimesters
#' 
#' @param trimester_results Output from run_trimester_analysis()
#' @return List with performance table, best models, and consistent CpGs
compare_trimester_results <- function(trimester_results) {
  
  all_metrics <- list()
  
  # Extract performance metrics for each trimester and model
  for(trim in names(trimester_results)) {
    trim_data <- trimester_results[[trim]]
    
    # Continuous outcome models
    cont_en <- trim_data$continuous$elastic_net$performance
    cont_en$Trimester <- trim
    cont_en$Outcome <- "Continuous"
    
    cont_rf <- trim_data$continuous$random_forest$performance
    cont_rf$Trimester <- trim
    cont_rf$Outcome <- "Continuous"
    
    # Dichotomous outcome models
    dicho_en <- trim_data$dichotomous$elastic_net$performance
    dicho_en$Trimester <- trim
    dicho_en$Outcome <- "Dichotomous"
    
    dicho_rf <- trim_data$dichotomous$random_forest$performance
    dicho_rf$Trimester <- trim
    dicho_rf$Outcome <- "Dichotomous"
    
    all_metrics[[trim]] <- rbind(cont_en, cont_rf, dicho_en, dicho_rf)
  }
  
  comparison_table <- do.call(rbind, all_metrics)
  
  # Identify best overall model for each outcome type
  best_continuous <- NULL
  best_dichotomous <- NULL
  
  if("R2_Mean" %in% names(comparison_table)) {
    cont_subset <- comparison_table[comparison_table$Outcome == "Continuous", ]
    best_continuous <- cont_subset[which.max(cont_subset$R2_Mean), ]
  }
  
  if("AUC_Mean" %in% names(comparison_table)) {
    dicho_subset <- comparison_table[comparison_table$Outcome == "Dichotomous", ]
    best_dichotomous <- dicho_subset[which.max(dicho_subset$AUC_Mean), ]
  }
  
  # Find CpGs consistently selected across all trimesters
  all_cpgs_en <- list()
  all_cpgs_rf <- list()
  
  for(trim in names(trimester_results)) {
    all_cpgs_en[[trim]] <- trimester_results[[trim]]$continuous$elastic_net$stable_cpgs
    # Take top 50 RF CpGs for comparison
    all_cpgs_rf[[trim]] <- trimester_results[[trim]]$continuous$random_forest$top_cpgs[1:50]
  }
  
  consistent_cpgs_en <- Reduce(intersect, all_cpgs_en)
  consistent_cpgs_rf <- Reduce(intersect, all_cpgs_rf)
  
  return(list(
    performance_table = comparison_table,
    best_continuous = best_continuous,
    best_dichotomous = best_dichotomous,
    consistent_cpgs_elastic_net = consistent_cpgs_en,
    consistent_cpgs_random_forest = consistent_cpgs_rf,
    n_consistent_en = length(consistent_cpgs_en),
    n_consistent_rf = length(consistent_cpgs_rf)
  ))
}

#' Visualize model performance comparison
#' 
#' @param results Output from run_complete_pipeline()
plot_model_comparison <- function(results) {
  library(ggplot2)
  library(tidyr)
  
  perf_data <- results$comparison$comparison_table
  
  # Plot R² for continuous outcomes
  if("R2_Mean" %in% names(perf_data)) {
    cont_data <- perf_data[perf_data$Outcome == "Continuous", ]
    
    p1 <- ggplot(cont_data, aes(x = Model, y = R2_Mean)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_errorbar(aes(ymin = R2_Mean - R2_SD, 
                        ymax = R2_Mean + R2_SD),
                    width = 0.2) +
      theme_minimal() +
      labs(title = "Model Performance - Continuous PM2.5",
           y = "R² (Mean ± SD)",
           x = "Model") +
      ylim(0, 1)
    
    print(p1)
  }
  
  # Plot AUC for dichotomous outcomes
  if("AUC_Mean" %in% names(perf_data)) {
    dicho_data <- perf_data[perf_data$Outcome == "Dichotomous", ]
    
    p2 <- ggplot(dicho_data, aes(x = Model, y = AUC_Mean)) +
      geom_bar(stat = "identity", fill = "darkred") +
      geom_errorbar(aes(ymin = AUC_Mean - AUC_SD, 
                        ymax = AUC_Mean + AUC_SD),
                    width = 0.2) +
      theme_minimal() +
      labs(title = "Model Performance - Dichotomous PM2.5",
           y = "AUC (Mean ± SD)",
           x = "Model") +
      ylim(0, 1)
    
    print(p2)
  }
}

#' Export all results to files
#' 
#' @param results Output from run_complete_pipeline()
#' @param beta_matrix Original beta matrix
#' @param output_dir Directory for output files
export_results <- function(results, beta_matrix, output_dir = "methylation_results") {
  
  dir.create(output_dir, showWarnings = FALSE)
  
  # Export elastic net CpGs with selection frequency and coefficients
  en_cpgs <- data.frame(
    CpG = results$continuous$elastic_net$selected_cpgs,
    Selection_Frequency = results$continuous$elastic_net$cpg_selection_frequency[
      results$continuous$elastic_net$selected_cpgs],
    Coefficient = results$continuous$elastic_net$coefficients[
      results$continuous$elastic_net$selected_cpgs]
  )
  
  write.csv(en_cpgs, 
            file.path(output_dir, "elastic_net_cpgs.csv"),
            row.names = FALSE)
  
  # Export random forest CpGs with importance scores
  rf_cpgs <- data.frame(
    CpG = names(results$continuous$random_forest$avg_importance),
    Importance = results$continuous$random_forest$avg_importance
  )
  
  write.csv(rf_cpgs,
            file.path(output_dir, "random_forest_cpgs.csv"),
            row.names = FALSE)
  
  # Export model performance comparison
  write.csv(results$comparison$comparison_table,
            file.path(output_dir, "model_performance.csv"),
            row.names = FALSE)
  
  # Export methylation profile scores for each sample
  profile_scores <- data.frame(
    Sample_ID = colnames(beta_matrix),
    ElasticNet_Score = results$continuous$elastic_net_score,
    RandomForest_Score = results$continuous$random_forest_score
  )
  
  write.csv(profile_scores,
            file.path(output_dir, "profile_scores.csv"),
            row.names = FALSE)
  
  message(paste("Results exported to", output_dir))
}

# ============================================
# EXAMPLE USAGE
# ============================================
# 
# # Load your data
# beta_matrix <- as.matrix(read.csv("beta_values.csv", row.names = 1))
# pm25_data <- read.csv("pm25_exposures.csv")
# 
# # Load literature CpGs (example structure)
# literature_studies <- list(
#   Gruzieva = read.table("gruzieva_cpgs.txt", stringsAsFactors = FALSE)$V1,
#   Wang = read.table("wang_cpgs.txt", stringsAsFactors = FALSE)$V1,
#   Parikh = read.table("parikh_cpgs.txt", stringsAsFactors = FALSE)$V1,
#   Broseus = read.table("broseus_cpgs.txt", stringsAsFactors = FALSE)$V1
# )
# 
# # Create dichotomous variable based on median split
# pm25_median <- median(pm25_data$overall_pregnancy)  # alternatively, use model derived threshold
# pm25_dichotomous <- factor(ifelse(pm25_data$overall_pregnancy >= pm25_median, 
#                                   "High", "Low"))
# 
# # Optional: prepare covariates (samples in rows, covariates in columns)
# covariates <- data.frame(
#   maternal_age = pm25_data$maternal_age,
#   gestational_age = pm25_data$gestational_age,
#   sex = as.numeric(pm25_data$sex == "Male")
# )
# 
# # Run complete pipeline
# results <- run_complete_pipeline(
#   beta_matrix = beta_matrix,
#   pm25_continuous = pm25_data$trimester1_avg,
#   pm25_dichotomous = pm25_dichotomous,
#   covariates = covariates,  # Optional
#   literature_studies = literature_studies,
#   k_outer = 5,
#   k_inner = 5,
#   var_percentile = 80,
#   pval_threshold = 0.05,
#   min_criteria = 2,
#   seed = 123
# )
# 
# # Analyze and export results
# analyze_filtering_results(results$continuous$filtering_details, beta_matrix)
# export_filtering_results(results$continuous$filtering_details, "continuous_filtering")
# export_results(results, beta_matrix, "methylation_results")
# 
# # Visualize results
# plot_filtering_venn(results$continuous$filtering_details)
# plot_model_comparison(results)