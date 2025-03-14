% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_adult_young_bootstrap_table.R
\name{refit_DMR_adult_young_random}
\alias{refit_DMR_adult_young_random}
\title{Create bootstrad of Young / adult  comparison}
\usage{
refit_DMR_adult_young_random(
  model_1,
  full_tables,
  path_to_annotation_table,
  n_boot = 1000
)
}
\arguments{
\item{model_1}{fitted model 1 on the full datasets}

\item{full_tables}{full 300 ws tables}

\item{path_to_annotation_table}{annotation table}

\item{n_boot}{number of bootstrap}
}
\value{
A data frame with effect size for DMR position in both adult and young
}
\description{
Create bootstrad of Young / adult  comparison
}
\examples{
\dontrun{
data("model_1_one_fold")
load("../../Downloads/Weyrich_EpiRank/full_300_ws_tables.rda")
path_to_annotation_table <- "../../Downloads/Weyrich_EpiRank/crocuta_liftoff_Hhy_ASM300989v1_addPromoter.gtf"
DMR_bootstrap_adult_young  <- refit_DMR_adult_young_random(model_1_one_fold, full_300_ws_tables, path_to_annotation_table, n_boot = 1000)
save(DMR_bootstrap_adult_young, file = "data/DMR_bootstrap_adult_young.rda")
}

}
