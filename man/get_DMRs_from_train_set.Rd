% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_DMRs.R
\name{get_DMRs_from_train_set}
\alias{get_DMRs_from_train_set}
\title{Wrapper around the get_DMRs() function. It first cut the count_table according to the data_train.
and then extract the DMR with get_DMRs()}
\usage{
get_DMRs_from_train_set(
  train_set,
  tables,
  formula,
  p_cut,
  coef,
  index = NULL,
  method = "GLM"
)
}
\arguments{
\item{train_set}{the talbe with the IDs in the train set.}

\item{tables}{the covariate and count tables created with prepare_ws_datasets()}

\item{formula}{the formula for the model}

\item{p_cut}{the p_value threshold on which to keep the DMR. (mostly used to keep object size decent)}

\item{coef}{the coeficient to extract the DMR for (can be multiple).}

\item{index}{a vector of intex to retain the value of the fit. THought for the CV to keep the value of the DMR in the full set}

\item{method}{GLM or QL}
}
\value{
A DF with DMRs
}
\description{
Wrapper around the get_DMRs() function. It first cut the count_table according to the data_train.
and then extract the DMR with get_DMRs()
}
