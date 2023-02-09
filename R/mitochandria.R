#' Compute mean cpm for mito and non-mito genes
#'
#' @param full_300_ws_tables full tables with count and cov tables
#' @param path_to_annotation_table path to annotated table stored in exdata
#' @return A tibble with covariate and mean cpm for mito genes and non-mito genes
#' @examples
#' \dontrun{
#' path_to_annotation_table <- "inst/extdata/crocuta_liftoff_Hhy_ASM300989v1_addPromoter.gtf"
#' load("data-raw/data/full_300_ws_tables.rda")
#'
#' out <- compute_mean_cpm_mito(full_300_ws_tables, path_to_annotation_table) ## TO run on the server as the files is around 30 GB
#'}
#' @export
compute_mean_cpm_mito <- function(full_300_ws_tables, path_to_annotation_table){

  mito <- c("ND1", "ND2", "ND3", "ND4", "ND5", "COX1", "COX2", "COX3")

  count_DNA <- add_annotation(path_to_annotation_table = path_to_annotation_table,
                              table = full_300_ws_tables$count_table, keep = "all") %>%
    dplyr::mutate(is_mito = purrr::map_lgl(.data$annotation, ~ any(.x$annot.gene_id %in% mito)))

  full_300_ws_tables$covariate_table %>%
    dplyr::mutate(mean_mito_TMM  = purrr::map_dbl(.data$count_names_TMM, ~ sum(count_DNA[count_DNA$is_mito, paste0(.x)], na.rm = T)/sum(count_DNA$is_mito)),  ## compute the mean CPM for  mito windows
                  mean_non_mito_TMM = purrr::map_dbl(.data$count_names_TMM, ~ sum(count_DNA[!count_DNA$is_mito, paste0(.x)], na.rm = T)/sum(!count_DNA$is_mito))) %>% ## compute the mean CPM for non mito windows
    dplyr::select(.data$ID, .data$group,  .data$mean_mito_TMM, .data$mean_non_mito_TMM) %>%
    tidyr::pivot_longer(dplyr::contains("mean"),names_to = "type_of_genes",values_to = "mean_count_TMM") %>%
    dplyr::mutate(type_of_genes = ifelse(.data$type_of_genes == "mean_mito_TMM", "mito", "non_mito"),
                  type_of_data  = "DNA")
}

#' Compute mean cpm for mito and non-mito genes RNA
#'
#' @param RNA_table full tables with RNA count and cov tables
#'
#' @return A tibble with covariate and mean cpm for mito genes and non-mito genes expression
#' @examples
#' \dontrun{
#' load("data-raw/data/RNA_table.rda")
#' mito_table_RNA <-compute_mean_cpm_mito_rna(RNA_table)
#'}
#' @export
compute_mean_cpm_mito_rna <- function(RNA_table){

  mito <- c("ND1", "ND2", "ND3", "ND4", "ND5", "COX1", "COX2", "COX3")
  count_RNA  <- RNA_table$count_table %>%
    dplyr::mutate(is_mito = .data$gene_id %in% mito)

cov_RNA <- RNA_table$covariate_table %>%
  dplyr::mutate(mean_mito_TMM = purrr::map_dbl(.data$count_names_TMM, ~ sum(count_RNA[count_RNA$is_mito, paste0(.x)], na.rm = T)/sum(count_RNA$is_mito)),
         mean_non_mito_TMM = purrr::map_dbl(.data$count_names_TMM, ~ sum(count_RNA[!count_RNA$is_mito, paste0(.x)], na.rm = T)/sum(!count_RNA$is_mito))) %>%
  dplyr::select(.data$ID, .data$group, .data$mean_mito_TMM, .data$mean_non_mito_TMM) %>%
  tidyr::pivot_longer(dplyr::contains("mean"),names_to = "type_of_genes",values_to = "mean_count_TMM") %>%
  dplyr::mutate(type_of_genes = ifelse(.data$type_of_genes == "mean_mito_TMM", "mito", "non_mito"),
         type_of_data = "RNA")
}
