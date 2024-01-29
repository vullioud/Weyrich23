#'full_300_ws_tables_1000
#'
#' A subset of the base table to run the analysis. It is composed of 2 data-frame.
#' The first consist of the read count, the normalised read count (using TMM),
#' for each samples and each windows of 300bp. This subset consist of the first
#' 1000 windows. In addition we have the [chr], [start], amd [stop] columns
#' refering to the scaffold, the first and last position of the windows.
#'
#' A subset of data
#' @usage data(full_300_ws_tables_1000)
#' @format `full_300_ws_tables_1000`
#' A list of 2 data-frames
#' \describe{
#'   \item{count_table}{data-frame with read count per samples per window}
#'   \item{covariate_table}{data-frame with information about the samples}
#' }
"full_300_ws_tables_1000"
NULL

#' model_1_one_fold
#'
#' The full set with its DMRs. It is organised as a nested tibble with 4 variables.
#' Used for the fitting of the models and the Random Forest.
#' Contains the result of the negative binomial model for the full dataset. It contains
#' one row only containing all the samples in both the train set and test set
#' variable. It have the same structure than the [model_1_six_fold] table.
#' @usage data(model_1_one_fold)
#' @format `model_1_one_fold`
#' A nested tibble with 1 row and 4 variables
#'  \describe{
#'   \item{set_ID}{the identifier of the set}
#'   \item{data_test}{A list column of  data-frame with the covariate table of the samples composing the set (all samples), The random forrest is trained on this set}
#'   \item{data_train}{The same list column as data_test}
#'   \item{DMR_train}{A list column of data-frame with the result of the neg-bin fit for the windows with significant difference between HR and LR}
#' }
"model_1_one_fold"
NULL


#' model_1_six_fold
#'
#' The six subsets with theirs DMRs. It is organised as a nested tibble with 4 variables,
#' and one row per subsets. Used for the robustness analysis of the DMRs.
#' Contains the result of the negative binomial model for the six subsets.
#' It contains six row each containing the samples used for the fitting in the
#' train set and the one removed from the analysis in the test set. The DMRs
#' are columns contains the position passing the cutof and the identified in
#' the full set.
#' It have the same structure than the [model_1_one_fold] table.
#'
#' @usage data(model_1_six_fold)
#' @format `model_1_six_fold`
#' A nested tibble with 1 row per subsets (6) and 4 variables
#' \describe{
#'   \item{set_ID}{the identifier of the set}
#'   \item{data_test}{A list column of  data-frame with the covariate table of the samples composing the subset}
#'   \item{data_train}{A list column of data-frame with the covariate table of the samples removed from the subset}
#'   \item{DMR_train}{A list column of data-frame with the result of the neg-bin fit for the windows with significant difference between HR and LR}
#' }
"model_1_six_fold"
NULL

#' DMR_group
#'
#' DMRs found in the full dataset with information about variance in the subsets
#' @usage data(DMR_group)
#' @format `DMR_group`
#' A nested tibble with 1 row per subset.
"DMR_group"
NULL

#' annotated_DMR_group
#'
#' DMRs found in the full dataset with information about variance in the subsets
#' @usage data(annotated_DMR_group)
#' @format `annotated_DMR_group`
#' A nested tibble with 1 row per subset.
"annotated_DMR_group"
NULL

#' loc_table
#'
#' table with loc names and corresponding annotation.
#' @usage data(loc_table)
#' @format `loc_table`
#' A tibble with loc names and gene names
"loc_table"
NULL

#' mito_table_DNA
#'
#' table with the mean cpm for mito and non-mito genes
#' fitted on a server with the same code as in the vignette (big files)
#' @usage data(mito_table_DNA)
#' @format `mito_table_DNA`
#' A tibble with mean cpm for mito and non mito genes per ID
"mito_table_DNA"
NULL

#' mito_table_RNA
#'
#' table with the mean cpm for mito and non-mito genes RNA
#'
#' @usage data(mito_table_RNA)
#' @format `mito_table_RNA`
#' A tibble with mean cpm for mito and non mito genes per ID
"mito_table_RNA"
NULL

#' DMR_adult_young_comp
#'
#' table with the effect size in the main DMR for young and adult fitted separately
#'
#' @usage data(DMR_full_comp)
#' @format `DMR_adult_young_comp`
#' DMRs found in the full dataset with information about variance in the subsets in adult and young sepparetely 
"DMR_adult_young_comp"
NULL
