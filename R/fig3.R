#' Create fig3
#'
#' @param mito_table_DNA table with mean cpm for mito and non mito genes DNA
#' @param mito_table_RNA table with mean cpm for mito and non mito genes RNA
#' @param save otption to save the file, nee
#' @param ... argument to pass to ggsave()
#' @return A plot
#' @export
#'
fig3 <- function(mito_table_DNA, mito_table_RNA, save = FALSE, ...){

plot_data <- dplyr::bind_rows(mito_table_DNA, mito_table_RNA)

plot <- ggplot2::ggplot(plot_data) +
  ggplot2::aes(x = .data$type_of_genes, y = .data$mean_count_TMM, col = .data$group) +
  ggplot2::geom_boxplot(fill = NA)+
  ggplot2::scale_y_log10("Mean Methylation") +
  ggplot2::scale_x_discrete("Genes", labels = c("Mitochondrial", "Nuclear")) +
  ggplot2::scale_color_manual("", values = c("#E69F00", "#56B4E9"), labels = c("High-ranking","Low-ranking")) +
  ggplot2::theme_classic() +
  ggplot2::facet_wrap( ~ .data$type_of_data, scales = "free")

if (save) {
ggplot2::ggsave(plot, ...) ## filename = "fig/manuscript_fig4.png", width = 150, height = 75, units = "mm"
  }
plot
}
