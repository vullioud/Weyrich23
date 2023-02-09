#' reshape the data for fig2.
#'
#' @param tables ttable with read_counts and covariates
#' @return A tibble with data for fig.1
#'
.reshape_data_fig_s2b <- function(tables){
  tables$covariate_table %>%
  dplyr::select(.data$group, .data$age_group, .data$lib.size,
                .data$norm.factors, .data$social_status, .data$age_at_sampling) %>%
  dplyr::mutate(adjusted_count = .data$norm.factors * .data$lib.size)

}

#' Create fig1
#'
#' @param tables table with read_counts and covariates
#' @param save otption to save the file, nee
#' @param ... argument to pass to ggsave()
#' @examples
#' \dontrun{
#' data(full_300_ws_tables_1000)
#' fig_s2b(full_300_ws_tables_1000,save= FALSE)
#' }
#' @return A plot
#' @export
#'
fig_s2b <- function(tables, save =FALSE, ...){
  cov_full <- .reshape_data_fig_s2b(tables)
plot <- ggplot2::ggplot(cov_full) +
  ggplot2::aes(y = .data$social_status, x = .data$age_at_sampling, col = .data$age_group, shape = .data$group, size = .data$adjusted_count) +
  ggplot2::geom_point() +
  ggplot2::geom_hline(yintercept = 0, col = "black") +
  ggplot2::geom_vline(xintercept = 2*365, col = "black") +
  ggplot2::geom_hline(yintercept = 0.33, col = "black", linetype = 2) +
  ggplot2::geom_hline(yintercept = -0.33, col = "black", linetype = 2) +
  ggplot2::scale_x_continuous("Age at sampling (in days)") +
  ggplot2::scale_y_continuous("Social status at sampling") +
  ggplot2::scale_shape_manual(name = "Rank group", values=c(3, 1)) +
  ggplot2::scale_color_manual(name = "Age group", values = c("#E69F00", "#56B4E9")) +
  ggplot2::scale_size(name = "Adjusted total read count") +
  ggplot2::theme_classic()

  if(save){
  ggplot2::ggsave(plot,...)
  }
 plot
}
