#' Create fig3
#'
#' @param mito_table_DNA table with mean cpm for mito and non mito genes DNA
#' @param mito_table_RNA table with mean cpm for mito and non mito genes RNA
#' @param save otption to save the file, nee
#' @param ... argument to pass to ggsave()
#' @return A plot
#' @export
#' @examples
#' \dontrun{
#' data("mito_table_DNA")
#' data("mito_table_RNA")
#' fig3(mito_table_DNA, mito_table_RNA, save = TRUE, filename = "figures/fig4_final.png", width = 8, height = 4, units = "in", dpi = 300)
#' fig3(mito_table_DNA, mito_table_RNA, save = TRUE, filename = "figures/fig4_final.pdf", width = 6.5, height = 3.2, units = "in", dpi = 300)
#' }
#'
fig3 <- function(mito_table_DNA, mito_table_RNA, save = FALSE, ...){
  
  plot_DNA <- ggplot2::ggplot(mito_table_DNA) +
    ggplot2::aes(x = .data$type_of_genes, y = .data$mean_count_TMM, col = .data$group) +
    ggplot2::geom_boxplot(fill = NA)+
    ggplot2::scale_y_log10("Mean Methylation [CpM]") +
    ggplot2::scale_x_discrete("Genes", labels = c("Mitochondrial", "Nuclear")) +
    ggplot2::scale_color_manual("", values = c("#E69F00", "#56B4E9")) +
    ggplot2::theme_classic() +
    ggplot2::theme(text =  ggplot2::element_text(size = 8), 
                   legend.position =  "none", 
                   axis.text =  ggplot2::element_text(size = 8), 
                   axis.title =  ggplot2::element_text(size = 8))
  
  ## plot RNA
  plot_RNA <- ggplot2::ggplot(mito_table_RNA) +
    ggplot2::aes(x = .data$type_of_genes, y = .data$mean_count_TMM, col = .data$group) +
    ggplot2::geom_boxplot(fill = NA)+
    ggplot2::scale_y_log10("Mean Expression [CpM]") +
    ggplot2::scale_x_discrete("Genes", labels = c("Mitochondrial", "Nuclear")) +
    ggplot2::scale_color_manual("", values = c("#E69F00", "#56B4E9"), labels = c("High-ranking","Low-ranking")) +
    ggplot2::theme_classic()+
    ggplot2::theme(text =  ggplot2::element_text(size = 8),
                   legend.position = c(0.82, 0.8), 
                   legend.text =  ggplot2::element_text(size = 8),
                   axis.text =  ggplot2::element_text(size = 8), 
                   axis.title =  ggplot2::element_text(size = 8))
  plot = cowplot::plot_grid(plot_DNA, plot_RNA, labels = c("a)", "b)"))
  if (save) {
    ggplot2::ggsave(plot, ...) ## filename = "fig/manuscript_fig4.png", width = 150, height = 75, units = "mm"
  }
  plot
  
}

