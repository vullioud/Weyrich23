#' reshape the data for fig2.
#'
#' @param DMR_group table with DMR
#'
#' @return A tibble with data for fig.1
#'
.reshape_data_fig2 <- function(DMR_group){
  DMR_group %>%
    dplyr::group_by(.data$gene_name) %>%
    dplyr::summarise(mean_fold_change = mean(.data$logFC_full_set),
                     min = min(.data$min_effect),
                     max = max(.data$max_effect),
                     coherent_block = !any(!.data$coherent_block)) %>%
    dplyr::filter(.data$gene_name != "intergenic") %>%
    dplyr::arrange(.data$mean_fold_change) %>%
    dplyr::mutate(plot_pos = as.factor(1:dplyr::n()))
}

#' reshape the data for fig2.
#'
#' @param DMR_group table with DMR
#'
#' @return A tibble with data for fig.S3
#'
.reshape_data_fig_s3 <- function(DMR_group){
  
  DMR_group %>%
    dplyr::filter(.data$gene_name == "intergenic") %>%
    dplyr::group_by(.data$block_num) %>%
    dplyr::summarise(mean_fold_change = mean(.data$logFC_full_set),
                     min = min(.data$min_effect),
                     max = max(.data$max_effect),
                     coherent_block = !any(!.data$coherent_block)) %>%
    dplyr::arrange(.data$mean_fold_change) %>%
    dplyr::mutate(plot_pos = as.factor(1:dplyr::n()))
}

#' Create fig2
#'
#' @param DMR_group table with DMR
#' @param save otption to save the file, nee
#' @param sup logical to create the suplementary table
#' @param ... argument to pass to ggsave()
#' @return A plot
#' @examples
#' \dontrun{
#' data(DMR_group)
#' fig2(DMR_group, save = TRUE, sup = FALSE, filename = "figures/fig2_final.png", width = 8, height = 8, units = "in", dpi = 300, bg = "white")
#' fig2(DMR_group, save = TRUE, sup = TRUE, filename = "figures/figS3_final.png", width = 8, height = 14, units = "in", dpi = 300, bg = "white")
#' }
#' @export
fig2 <- function(DMR_group, save =FALSE, sup= FALSE, ...){
  
  if(sup) {
    data <- .reshape_data_fig_s3(DMR_group)
  } else {
    data <- .reshape_data_fig2(DMR_group)
  }
  plot <- ggplot2::ggplot(data) +
    ggplot2::aes(x =.data$mean_fold_change, y = .data$plot_pos) +
    ggplot2::geom_point(shape = 3) +
    ggplot2::geom_linerange(ggplot2::aes(xmin= .data$min, xmax= .data$max, col = .data$coherent_block), linetype=3) +
    ggplot2::geom_point(ggplot2::aes(x= .data$min, col = .data$coherent_block),size=1, shape = 4) +
    ggplot2::geom_point(ggplot2::aes(x= .data$max, col = .data$coherent_block),size=1, shape = 4) +
    ggplot2::scale_y_discrete("Genes with rankDMRs", labels = data$gene_name) +
    ggplot2::geom_vline(xintercept = 0, linetype = 3) +
    ggplot2::scale_x_continuous("Effect size [logFC]") +
    ggplot2::scale_color_manual("", labels = c("Excluded", "Selected"), values = c("#E69F00", "#56B4E9")) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = c(0.90, 0.1), 
                   legend.text =  ggplot2::element_text(size = 8))
  
  if(sup){
    plot <- plot + ggplot2::scale_y_discrete("rankDMRs", labels = data$block_num)
  }
  if (save) {
    ggplot2::ggsave(plot, ...) ##filename = "fig/manuscript_fig2.png", width = 130, height = 130, units = "mm")  // sup: filename = "fig/manuscript_sup_figS3.png", width = 120, height = 200, units = "mm"
  }
  plot
}
###############