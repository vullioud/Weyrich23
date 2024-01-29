#' Create bootstrad of Young / adult  comparison 
#'
#' @param model_1 fitted model 1 on the full datasets 
#' @param path_to_annotation_table annotation table 
#' @param full_tables full 300 ws tables 
#' @param n_boot number of bootstrap 
#' @return A data frame with effect size for DMR position in both adult and young 
#' @examples
#' \dontrun{
#' data("model_1_one_fold")
#' load("../../Downloads/Weyrich_EpiRank/full_300_ws_tables.rda")
#' path_to_annotation_table <- "../../Downloads/Weyrich_EpiRank/crocuta_liftoff_Hhy_ASM300989v1_addPromoter.gtf"
#' DMR_bootstrap_adult_young  <- refit_DMR_adult_young_random(model_1_one_fold, full_300_ws_tables, path_to_annotation_table, n_boot = 1000)
#' save(DMR_bootstrap_adult_young, file = "data/DMR_bootstrap_adult_young.rda")
#' }
#'
refit_DMR_adult_young_random <- function(model_1, full_tables, path_to_annotation_table, n_boot = 1000){
  
  
  DMR_full <- reshape_DMR_table(set = model_1, coef = "groupLR", p_cut = 1, q_cut = 0.05)
  
  full_300 <- full_tables
  out <- list()
  
  for (i in 1:n_boot){
    full_tables$count_table <- dplyr::sample_n(full_300$count_table, size = 321)
    #
    # ## fit for adult
    full_300_ws_tables_adult <- full_tables
    full_300_ws_tables_adult$covariate_table <- full_300_ws_tables_adult$covariate_table[full_300_ws_tables_adult$covariate_table$age_group == "adult", ]
    #
    #
    model_1_one_fold_adult <- prepare_balanced_k_set(k = 1, variable_to_balance = "group", tables =full_300_ws_tables_adult)
    model_1_one_fold_adult <- add_DMRs_to_set(model_1_one_fold_adult,
                                              full_300_ws_tables_adult, formula = "~ group", #
                                              p_cut = 1)
    DMR_adult_random <- model_1_one_fold_adult$DMR_train[[1]] %>%
      dplyr::rename(logFC_adult = .data$logFC,
                    PValue_adult = .data$PValue) %>% 
      dplyr::select(-.data$FDR, -.data$coef, -.data$index)
    #
    #
    # ################################## same for yougster
    #
    # ## fit for adult
    #
    full_300_ws_tables_young <- full_tables
    full_300_ws_tables_young$covariate_table <- full_300_ws_tables_young$covariate_table[full_300_ws_tables_young$covariate_table$age_group != "adult", ]
    #
    #
    model_1_one_fold_young <- prepare_balanced_k_set(k = 1, variable_to_balance = "group", tables =full_300_ws_tables_young)
    #
    # ## fit negative binomial on the windows
    model_1_one_fold_young <- add_DMRs_to_set(model_1_one_fold_young,
                                              full_300_ws_tables_young, formula = "~ group", #
                                              p_cut = 1)
    DMR_young_random <- model_1_one_fold_young$DMR_train[[1]] %>%
      dplyr::rename(logFC_young = .data$logFC,
                    PValue_young = .data$PValue) %>% 
      dplyr::select(-.data$FDR, -.data$coef, -.data$index)
    
    
    DMR_full_comp_random <- DMR_adult_random %>% 
      dplyr::left_join(DMR_young_random) %>%
      dplyr::mutate(logFC_diff = .data$logFC_adult - .data$logFC_young)
    
    
    x <- stats::cor.test(DMR_full_comp_random$logFC_young, DMR_full_comp_random$logFC_adult, method = "pearson")
    out1 <- data.frame(cor = x$estimate,
                       p = x$p.value)
    out[[i]] <- out1
  }
  dplyr::bind_rows(out)
}
