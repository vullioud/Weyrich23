% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/random_forrest.R
\name{summarize_cv}
\alias{summarize_cv}
\title{Provide summary statistic on the performance of the classifier.
return the accuracy, the f1-score, the younden J-statistic, the recall and precision rates.}
\usage{
summarize_cv(rf_outcome_table)
}
\arguments{
\item{rf_outcome_table}{the output table of the rf fitted on the cv set with run_rf_on_set()}
}
\value{
A tibble with the summary statistics for the predictive power of the RF
}
\description{
Provide summary statistic on the performance of the classifier.
return the accuracy, the f1-score, the younden J-statistic, the recall and precision rates.
}
