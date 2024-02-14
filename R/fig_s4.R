#' Create fig_S4
#'
#' @param DMR_group table with DMR on the full set
#' @param DMR_adult_young_comp table with Adult/youg fitted separately on DMR
#' @param DMR_bootstrap_adult_young table with adult/young fitted separately on random windows 
#' @param save otption to save the file, nee
#' @param ... argument to pass to ggsave()
#' @examples
#' \dontrun{
#'data("DMR_group")
#'data("DMR_adult_young_comp")
#'data("DMR_bootstrap_adult_young")
#'fig_s4(DMR_group, DMR_adult_young_comp, DMR_bootstrap_adult_young, save = TRUE, filename = "figures/figS4_final.png", device = "png", width = 12, units = "in", height = 8, dpi = 300, bg = "white")
#'fig_s4(DMR_group, DMR_adult_young_comp, DMR_bootstrap_adult_young, save = TRUE, filename = "figures/figS4_final.pdf", device = "pdf", width = 9, units = "in", height = 6, dpi = 300, bg = "white")

#' }
#' @return A plot
#' @export
fig_s4 <- function(DMR_group, DMR_adult_young_comp, DMR_bootstrap_adult_young, save = FALSE, ...){
  ## create sets 
  DMR_group <- DMR_group %>%
    dplyr::select(1:3, "block", "coherent_block", "gene_name") %>% 
    dplyr::distinct()
  
  DMR_comp <- DMR_adult_young_comp %>% 
    dplyr::left_join(DMR_group)  %>%
    dplyr::mutate(dif_adult = .data$logFC_full_set - .data$logFC_adult,
                  dif_young = .data$logFC_full_set - .data$logFC_young)
  
  
  DMR_comp_annot <- DMR_comp %>% 
    dplyr::filter(!.data$annot_null)
  
  length(unique(DMR_comp_annot$gene_name))
  
  #stats::cor.test(DMR_comp$logFC_young, DMR_comp$logFC_adult, method = "pearson")
  
  
  adult_full <- DMR_comp_annot %>%  dplyr::filter((.data$logFC_adult >= 0 & .data$logFC_young <= 0 & .data$logFC_full_set >= 0) | (.data$logFC_adult <= 0 & .data$logFC_young >= 0 & .data$logFC_full_set <= 0))
  young_full <- DMR_comp_annot %>%  dplyr::filter((.data$logFC_adult >= 0 & .data$logFC_young <= 0 & .data$logFC_full_set <= 0) | (.data$logFC_adult <= 0 & .data$logFC_young >= 0 & .data$logFC_full_set >= 0))
  ## flag the DMR going in the direction of adult // There is none going in the direction of youg
  
  adult_full$same_dir <- "Adult"
  x <- DMR_comp_annot %>% 
    dplyr::left_join(adult_full) 
  
  x$same_dir[is.na(x$same_dir)] <- "none"
  
  ### for annotated
  plot1 <- ggplot2::ggplot(x) +
    ggplot2::geom_text(ggplot2::aes(y = .data$logFC_adult, x = .data$logFC_young, col = .data$same_dir, label = .data$gene_name), angle = 0, size = 2.5) +
    ggplot2::scale_x_continuous("Effect size in cubs [logFC]", limits = c(-5, 6)) +
    ggplot2::scale_y_continuous("Effect size in adults [logFC]", limits = c(-5, 6)) +
    ggplot2::scale_color_manual("", labels = c("Incoherent effect direction", "Coherent effect direction"), values = c("#E69F00", "#56B4E9")) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2)+
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = c(0.7, 0.2)) +
    ggplot2::theme(text = ggplot2::element_text(size =8),
                   legend.text = ggplot2::element_text(size =8),
                   legend.title = ggplot2::element_text(size = 8),
                   axis.text = ggplot2::element_text(size = 8),
                   axis.title = ggplot2::element_text(size = 8))
  
  
  plot2 <- ggplot2::ggplot(DMR_comp) +
    ggplot2::geom_point(ggplot2::aes(x = .data$logFC_young, y = .data$logFC_adult, col = .data$coherent_block), size = 1) +
    ggplot2::scale_x_continuous("Effect size in cubs [logFC]", limits = c(-5, 6)) +
    ggplot2::scale_y_continuous("Effect size in adults [logFC]", limits = c(-5, 6)) +
    ggplot2::scale_color_manual("", labels = c("Excluded rankDMR", "Selected rankDMR"), values = c("#E69F00", "#56B4E9")) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2)+
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = c(0.7, 0.2))+
    ggplot2::theme(text = ggplot2::element_text(size =8), 
                   legend.text = ggplot2::element_text(size =8),
                   legend.title = ggplot2::element_text(size = 8),
                   axis.text = ggplot2::element_text(size = 8),
                   axis.title = ggplot2::element_text(size = 8))
  
  
  ########################### run som bootstrap to get the correlation 
  
  
  plot3 <- ggplot2::ggplot(DMR_bootstrap_adult_young) +
    ggplot2::aes(x = .data$cor) +
    ggplot2::geom_density() +
    ggplot2::geom_vline(xintercept = 0.4799619, col = "#E69F00") + ## 0.47.. represents the observed cor. obtained by the commented line above. 
    ggplot2::scale_x_continuous("Coef. of correlation", limits = c(-0.2, 0.7)) +
    ggplot2::scale_y_continuous("Density") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size =8),
                   legend.text = ggplot2::element_text(size =8),
                   legend.title = ggplot2::element_text(size = 8),
                   axis.text = ggplot2::element_text(size = 8),
                   axis.title = ggplot2::element_text(size = 8))
  
  plot <- cowplot::plot_grid(plot2, plot3, plot1,  labels = c("a)", "b)", "c)"))
  
  if (save) {
    ggplot2::ggsave(plot, ...)
  }
  plot
}
