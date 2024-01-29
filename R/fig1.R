#' reshape the data for fig1.
#'
#' @param DMR_group table with DMR
#'
#' @return A tibble with data for fig.1
#'
.reshape_data_fig1 <- function(DMR_group){
  DMR_group %>%
    dplyr::group_by(.data$block) %>%
    dplyr::mutate(block_length = length(unique(.data$index))) %>%
    dplyr::group_by() %>%
    dplyr::mutate(effect = ifelse(.data$logFC_full_set >= 0, "hyper", "hypo")) %>%
    dplyr::distinct(.data$annot_null, .data$effect, .data$block, .data$block_length) %>%
    dplyr::group_by(.data$block_length) %>%
    dplyr::summarise(n_hypo_annot = -sum(!.data$annot_null & .data$effect == "hypo"),
                     n_hyper_annot = sum(!.data$annot_null & .data$effect == "hyper"),
                     n_hypo = -sum(.data$effect == "hypo"),
                     n_hyper = sum(.data$effect == "hyper")) %>%
    tidyr::pivot_longer(-.data$block_length) %>%
    dplyr::mutate(color = ifelse(.data$name %in% c("n_hypo_annot", "n_hyper_annot"), "annotated", "intergenic"))
  
}


#' Create fig1
#'
#' @param DMR_group table with DMR
#' @param save otption to save the file, nee
#' @param ... argument to pass to ggsave()
#' @return A plot
#' @export
#' @examples
#' \dontrun{
#' data("DMR_group")
#'  fig1(DMR_group, save = TRUE, filename = "figures/fig1_final.png", width = 12, height = 8, units = "in", dpi = 300, bg = "white")
#' }
#' 
#'
fig1 <- function(DMR_group, save = FALSE, ...){
  data <- .reshape_data_fig1(DMR_group)
  
  plot <- ggplot2::ggplot(data) +
    ggplot2::aes(x = .data$block_length, y = .data$value,  fill = .data$color) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
    ggplot2::scale_x_continuous("rankDMRs length (bp)", breaks = c(1,2,3,4,5,6,7,8,9,10,11,12),
                                labels = c("300", "600", "900", "1200", "1500", "1800", "2100", "2400", "2700", "3000", "3300", "3600")) +
    ggplot2::scale_y_continuous(expression(paste("Count of rankDMRs")), breaks = c(-25, 0, 25, 50, 75, 100),
                                labels = c(25,0,25,50, 75, 100), limits = c(-40, 150)) +
    ggplot2::scale_fill_manual("", labels = c("Intragenic", "Intergenic"), values = c("#E69F00", "#56B4E9"))+
    ggplot2::annotate(geom = "text", y = 50, x = 10, label= "Hypermethylated", size = 5) +
    ggplot2::annotate(geom = "text", y = -20, x = 10, label= "Hypomethylated", size = 5) +
    ggplot2::theme_classic()  +
    ggplot2::theme(axis.text.x= ggplot2::element_text(angle = 60, hjust = 1)) +
    ggplot2::theme(legend.position = c(0.80, 0.80), 
                   legend.text = ggplot2::element_text(size =15),
                   axis.text = ggplot2::element_text(size = 12),
                   axis.title = ggplot2::element_text(size = 15))
  if(save){
    ggplot2::ggsave(plot, ...) 
  }
  plot
}