
#' describe the samples
#'
#' Give the mean, n, and sd for group and age_group
#'
#'
#' @param tables table with read_counts and covariates
#' @return a list of tibble with samples descriptiomns
#' @export
#'
describes_samples <- function(tables){

cov_full <- tables$covariate_table %>%
  dplyr::select(.data$group, .data$age_group, .data$lib.size, .data$norm.factors, .data$social_status, .data$age_at_sampling) %>%
  dplyr::mutate(adjusted_count = .data$norm.factors * .data$lib.size)

out <- list()

out$descriptive_stats <- cov_full %>%
  dplyr::group_by(.data$age_group, .data$group) %>%
  dplyr::summarise(n = dplyr::n(),
            mean_age = mean(.data$age_at_sampling),
            mean_rank = mean(.data$social_status),
            sd_age = stats::sd(.data$age_at_sampling),
            sd_rank = stats::sd(.data$social_status))

out$descriptive_stats_rank <- cov_full %>%
  dplyr::group_by(.data$group) %>%
  dplyr::summarise(n = dplyr::n(),
            mean_age = mean(.data$age_at_sampling),
            mean_rank = mean(.data$social_status),
            sd_age = stats::sd(.data$age_at_sampling),
            sd_rank = stats::sd(.data$social_status))

out$descriptive_stats_age <- cov_full %>%
  dplyr::group_by(.data$age_group) %>%
  dplyr::summarise(n = dplyr::n(),
            mean_age = mean(.data$age_at_sampling),
            mean_rank = mean(.data$social_status),
            sd_age = stats::sd(.data$age_at_sampling),
            sd_rank = stats::sd(.data$social_status))

out

}
