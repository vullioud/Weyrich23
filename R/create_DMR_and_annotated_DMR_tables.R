#' Provide summary statistic on the performance of the classifier.
#' return the accuracy, the f1-score, the younden J-statistic, the recall and precision rates.
#'
#' @param DMR_full internal DMR table with indexes to merge
#'
#' @return A tibble 2 added col. block and block_num

.add_merged_DMR <- function(DMR_full) {

  indexes <- DMR_full$index

  ### just check how many of the indexes are merge add block !!! TO DO  be done in a clean way
  index_tables <- data.frame(index_base = DMR_full$index) %>%
    dplyr::mutate(index2 = dplyr::lag(.data$index_base),
           follow = .data$index_base == (.data$index2 +1)) ## check that the index + 1 is the same as the next DMR index (i.e. they are adjacent)

  ## add block of DMR ids
  out <- vector(mode = "numeric", length = length(index_tables$follow))
  out[1] <- 1 ## set the first block ID
  for(i in 2:length(out)){

    if(!index_tables$follow[i]) out[i] <- out[(i-1)] +1 else out[i] <- out[i-1]  ## get the same ID as previous DMRs if it follows it oitherwise + 1  (ugly but it works)
  }
  index_tables$block <- paste0("block_", out)
  DMR_full$block <- index_tables$block ## add block ID to the table
  DMR_full$block_num <- as.factor(as.numeric(stringr::str_remove(DMR_full$block,"block_"))) ## add a numerical ID to the table.
  DMR_full
}

#################################################################################

#' Provide summary statistic on the performance of the classifier.
#' return the accuracy, the f1-score, the younden J-statistic, the recall and precision rates.
#'
#' @param DMR_CV_table merged table with DMR from subsets and full dataset.
#'
#' @return A tibble 2 added col. block and block_num

.identify_incoherent_DMR <- function(DMR_CV_table){
  DMR_CV_table %>%
    dplyr::group_by(.data$block) %>%
    dplyr::mutate(effect_both_dir = any(.data$logFC <= 0) & any(.data$logFC >= 0))  %>%
    dplyr::group_by(.data$block) %>%
    dplyr::mutate(coherent_block = !any(.data$effect_both_dir | .data$max_p >= 0.05),
                  block_length = dplyr::n()/6)
}

#################################################################################

#' Create a clean table with DMRs from the full set with information about CV subsets.
#' return the accuracy, the f1-score, the younden J-statistic, the recall and precision rates.
#'
#' @param model_1_one_fold output of the fitting of model 1 fold
#' @param model_1_six_fold output of the fitting of model 6 fold
#' @param path_to_annotation_table path to annotation table stored in int/exdata
#' @param loc_table internal DMR table with indexes to merge
#'
#' @return A tibble with DMR with full DMRs information
#' @examples
#' \dontrun{
#' data(model_1_one_fold)
#' data(model_1_six_fold)
#' path_to_annotation_table <- "~/Documents/Weyrich23/inst/extdata/crocuta_liftoff_Hhy_ASM300989v1_addPromoter.gtf"
#' data(loc_table)
#' create_DMR_group_table(model_1_one_fold, model_1_six_fold, path_to_annotation_table, loc_table)
#' }
#' @export
create_DMR_group_table <- function(model_1_one_fold, model_1_six_fold, path_to_annotation_table, loc_table){

DMR_full <- reshape_DMR_table(set = model_1_one_fold, coef = "groupLR", p_cut = 1, q_cut = 0.05)

DMR_full <- add_annotation(path_to_annotation_table, table = DMR_full, keep = "all") %>%
  dplyr::mutate(annot_null = purrr::map_lgl(.data$annotation, ~ is.null(.x))) %>%
  dplyr::rename(logFC_full_set = .data$logFC,
         Pvalue_full_set = .data$PValue,
         FDR_full_set = .data$FDR) %>%
  dplyr::select(-.data$presence)
## first get all the DMR for rank

DMR_full <- .add_merged_DMR(DMR_full) ## add block info
indexes <- DMR_full$index
# merge now with the information about the CV
model_1_six_fold$DMR_train <- purrr::map2(model_1_six_fold$set_ID, model_1_six_fold$DMR_train, ~.y %>%
                                            dplyr::mutate(set_ID = .x)) ## add the set ID to the DMR table col.

DMR_CV2 <- dplyr::bind_rows(model_1_six_fold$DMR_train) %>%
  dplyr::filter(.data$index %in% indexes & .data$coef == "groupLR")%>%
  dplyr::group_by(.data$index) %>%
  dplyr::mutate(presence = dplyr::n(),
         presence_q = sum(.data$FDR <= 0.05),
         sd_effect = stats::sd(.data$logFC),
         sd_p = stats::sd(.data$PValue),
         sd_FDR = stats::sd(.data$FDR),
         mean_FDR = mean(.data$FDR),
         mean_P = mean(.data$PValue),
         mean_effect = mean(.data$logFC),
         range_effect = range(.data$logFC)[2] - range(.data$logFC)[1],
         range_p = range(.data$PValue)[2]- range(.data$PValue)[1],
         range_FDR = range(.data$FDR)[2] - range(.data$FDR)[1],
         min_effect = range(.data$logFC)[1],
         max_effect = range(.data$logFC)[2],
         min_p = range(.data$PValue)[1],
         max_p = range(.data$PValue)[2],
         min_q = range(.data$FDR)[1],
         max_q = range(.data$FDR)[2]) %>%
  dplyr::left_join(DMR_full) %>%
  dplyr::group_by() %>%
  dplyr::mutate(annotation = purrr::map(.data$annotation, ~ .rename_loc(.x, loc_table = loc_table)),
                gene_name = purrr::map_chr(.data$annotation, ~ if(is.null(.x)) "intergenic" else paste(unique(.x$annot.gene_id), collapse = "/")))

.identify_incoherent_DMR(DMR_CV2) %>% ##   flag incohenrent block
  dplyr::group_by()

}

#' Make a long version of the annotated DMR_group_table
#'
#'remove the intergenic region and keep 1 row per annotation so 1 DMR can represent
#'many rows.
#'
#' @param DMR_table The DMR table
#'
#' @return A long format tibble with annotated DMR
#' @export
create_annot_DMRs <- function(DMR_table){
  DMR_table %>%
    dplyr::filter(!.data$annot_null & .data$coherent_block) %>%
    tidyr::unnest_longer(.data$annotation) %>%
    tidyr::unnest(.data$annotation)
}


#' Get the annotation type
#'
#' Add intron information. Region are considered intron if they overlap with gene and no exon.
#'
#' @param annotated_DMR The annotated DMR table
#'
#' @return A tibble with gene annotation including intron
#' @export
get_annotation_type <- function(annotated_DMR){
annotated_DMR %>%
    dplyr::select(.data$block, .data$block_length, .data$coherent_block, .data$annot.type, .data$annot.gene_id) %>%
    dplyr::distinct() %>%
    dplyr::group_by(.data$annot.gene_id) %>%
    dplyr::mutate(is_intron = ifelse(any(c("gene") %in% .data$annot.type & !(c("exon") %in% .data$annot.type)), TRUE, FALSE)) %>% ## the rule is if in gene but no overlap with exon = intron
    dplyr::mutate(annot.type = ifelse(.data$annot.type == "gene" & .data$is_intron, "intron", paste0(.data$annot.type))) %>%
    dplyr::filter(!.data$annot.type  %in% c("transcript", "CDS", "gene"))
}


