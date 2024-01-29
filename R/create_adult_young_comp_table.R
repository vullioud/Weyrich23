#' Create fig3
#'
#' @param model_1 table with mean cpm for mito and non mito genes DNA
#' @param path_to_annotation_table table with mean cpm for mito and non mito genes RNA
#' @param full_tables otption to save the file, nee
#' @return A data frame with effect size for DMR position in both adult and young 
#' @examples
#' \dontrun{
#' data("model_1_one_fold")
#' load("../../Downloads/Weyrich_EpiRank/full_300_ws_tables.rda")
#' path_to_annotation_table <- "../../Downloads/Weyrich_EpiRank/crocuta_liftoff_Hhy_ASM300989v1_addPromoter.gtf"
#' DMR_adult_young_comp  <- refit_DMR_adult_young_separately(model_1_one_fold, full_300_ws_tables, path_to_annotation_table)
#' save(DMR_adult_young_comp, file = "data/DMR_adult_young_comp.rda")
#' }
#'
refit_DMR_adult_young_separately <- function(model_1, full_table, path_to_annotation_table){
  
  DMR_full <- reshape_DMR_table(set = model_1, coef = "groupLR", p_cut = 1, q_cut = 0.05)
  
  DMR_full <- add_annotation(path_to_annotation_table, table = DMR_full, keep = "all") %>%
    dplyr::mutate(annot_null = purrr::map_lgl(.data$annotation, ~ is.null(.x))) %>%
    dplyr::rename(logFC_full_set = .data$logFC,
                  Pvalue_full_set = .data$PValue,
                  FDR_full_set = .data$FDR) %>%
    dplyr::select(-.data$presence)
  
  ## first get all the DMR for rank
  DMR_full <- .add_merged_DMR(DMR_full) ## add block info
  indexes <- DMR_full$index
  
  
  # extract the read counts for the DMR
  full_table$count_table <- full_table$count_table[indexes, ]
  
  ## fit for adult
  full_300_ws_tables_adult <- full_table
  full_300_ws_tables_adult$covariate_table <- full_300_ws_tables_adult$covariate_table[full_300_ws_tables_adult$covariate_table$age_group == "adult", ]
  
  
  model_1_one_fold_adult <- prepare_balanced_k_set(k = 1, variable_to_balance = "group", tables =full_300_ws_tables_adult)
  
  ## fit negative binomial on the windows
  model_1_one_fold_adult <- add_DMRs_to_set(model_1_one_fold_adult,
                                            full_300_ws_tables_adult, formula = "~ group", #
                                            p_cut = 1)
  DMR_adult <- model_1_one_fold_adult$DMR_train[[1]] %>%
    dplyr::rename(logFC_adult = .data$logFC,
                  PValue_adult = .data$PValue) %>% 
    dplyr::select(-.data$FDR, -.data$coef, -.data$index)
  
  
  ################################## same for yougster
  
  ## fit for yougnster
  
  full_300_ws_tables_young <- full_table
  full_300_ws_tables_young$covariate_table <- full_300_ws_tables_young$covariate_table[full_300_ws_tables_young$covariate_table$age_group != "adult", ]
  
  
  model_1_one_fold_young <- prepare_balanced_k_set(k = 1, variable_to_balance = "group", tables =full_300_ws_tables_young)
  
  ## fit negative binomial on the windows
  model_1_one_fold_young <- add_DMRs_to_set(model_1_one_fold_young,
                                            full_300_ws_tables_young, formula = "~ group", #
                                            p_cut = 1)
  DMR_young <- model_1_one_fold_young$DMR_train[[1]] %>%
    dplyr::rename(logFC_young = .data$logFC,
                  PValue_young = .data$PValue) %>% 
    dplyr::select(-.data$FDR, -.data$coef, -.data$index)
  
  DMR_full %>% dplyr::left_join(DMR_adult) %>% 
    dplyr::left_join(DMR_young) %>%
    dplyr::mutate(logFC_diff = .data$logFC_adult - .data$logFC_young)
}
