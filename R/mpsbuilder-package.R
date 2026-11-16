#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom caret createFolds
#' @importFrom foreach %dopar% foreach
#' @importFrom ggplot2 aes geom_bar geom_errorbar ggplot labs theme_minimal
#'   ylim
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom grid grid.draw grid.newpage
#' @importFrom matrixStats rowVars
#' @importFrom pROC auc roc
#' @importFrom progress progress_bar
#' @importFrom ranger ranger
#' @importFrom stats coef cor cor.test median predict quantile sd setNames
#'   t.test
#' @importFrom utils head write.csv write.table
#' @importFrom VennDiagram venn.diagram
#' @importFrom xgboost xgb.attr xgb.DMatrix xgb.importance xgb.train
## usethis namespace: end
NULL

# Column names used in ggplot2::aes() in plot_model_comparison() and the
# foreach() loop variable in filter_cpg_sites()
utils::globalVariables(
  c("AUC_Mean", "AUC_SD", "Model", "R2_Mean", "R2_SD", "chunk")
)
