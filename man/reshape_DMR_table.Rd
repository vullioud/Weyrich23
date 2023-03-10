% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_DMRs.R
\name{reshape_DMR_table}
\alias{reshape_DMR_table}
\title{extract the DMRs for the given covariate and add the presence of CV.}
\usage{
reshape_DMR_table(set, coef, p_cut, q_cut)
}
\arguments{
\item{set}{the table created with prepare_balanced_k_set().}

\item{coef}{the coeficient on which to extract the DMRs.}

\item{p_cut}{the p_value threshold on which to keep the DMR. (if filter on q_value set to 1)}

\item{q_cut}{the FDR threshold on which to keep the DMR. (if filter on pvalue set to 1)}
}
\value{
A nested DF with DMRs for each subset
}
\description{
extract the DMRs for the given covariate and add the presence of CV.
}
\examples{
data(model_1_one_fold)
reshape_DMR_table(model_1_one_fold, coef = "groupLR", p_cut = 1, q_cut = 0.05)
}
