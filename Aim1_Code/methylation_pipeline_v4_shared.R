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
library(mpsbuilder) # functions previously defined in this script

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

# create_literature_indicator() moved to R/filter_cpgs.R

# filter_cpg_sites() moved to R/filter_cpgs.R


## Additional (optional) functions for analyzing filtering results ##
# analyze_filtering_results() moved to R/filter_cpgs.R

# export_filtering_results() moved to R/filter_cpgs.R

# ============================================
# STEP 2: NESTED CROSS-VALIDATION SETUP
# ============================================

# create_nested_folds() moved to R/nested_folds.R

# ============================================
# STEP 3: ELASTIC NET WITH NESTED CV
# ============================================

# create_foldid() moved to R/nested_folds.R

# elastic_net_nested_cv() moved to R/elastic_net.R

# ============================================
# STEP 4: RANDOM FOREST WITH NESTED CV
# ============================================

# random_forest_nested_cv_tuned() moved to R/random_forest.R

# ============================================
# STEP 5: CREATE METHYLATION PROFILE SCORE
# ============================================

# create_methylation_profile_score() moved to R/profile_score.R

# ============================================
# STEP 6: COMPREHENSIVE MODEL COMPARISON
# ============================================

# compare_all_models() moved to R/compare_models.R

# compare_filtering_results() moved to R/compare_models.R

# ============================================
# MAIN PIPELINE WRAPPER
# ============================================

# run_complete_pipeline() moved to R/pipeline.R

# ============================================
# VISUALIZATION AND ANALYSIS FUNCTIONS
# ============================================

# plot_filtering_venn() moved to R/filtering_analysis.R

# sensitivity_analysis_filtering() moved to R/filtering_analysis.R

# compare_temporal_filtering() moved to R/filtering_analysis.R

# analyze_literature_cpgs() moved to R/filtering_analysis.R

# ============================================
# TRIMESTER ANALYSIS FUNCTIONS
# ============================================

# run_trimester_analysis() moved to R/trimester_analysis.R

# compare_trimester_results() moved to R/trimester_analysis.R

# plot_model_comparison() moved to R/results_output.R

# export_results() moved to R/results_output.R

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