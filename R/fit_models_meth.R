#' Extract the count of a count table
#'
#' @param count_table A data frame of counts. Idealy created with the
#' Epiena::prepare_ROI_datasets()]or Epiena::prepare_ws_datasets() from the Epiena package.  Each count column must ends with "count".
#' The first three col of the table shoup be the "chr", "start" and "stop".
#' @param which either "count" or "count_TMM". It should match the end of the count columns.
#' @return A named matrix countaing the read counts.
#' @export
extract_count_matrix <- function(count_table, which = "count") {

  fcountDF <- count_table %>%
    dplyr::select(1:3, dplyr::ends_with(which)) %>%
    dplyr::mutate(name = paste(as.character(.data$chr), as.character(.data$start), as.character(.data$stop), sep = "_")) %>%
    dplyr::select(.data$name, dplyr::everything(), -.data$chr, -.data$start, -.data$stop)

  count_matrix <- as.matrix(fcountDF[, 2:ncol(fcountDF)])
  row.names(count_matrix) <- fcountDF$name
  count_matrix
}


#'  Extract the count of a count table
#'
#' @param covariate_table A data frame with the ID col.
#' @return A data.frame with the same column as the input plus a col with the
#'  corresponding names in the count_table.
#' @export
add_count_colnames_to_covariates_tables <- function(covariate_table){
  covariate_table %>%
    dplyr::mutate(count_names = paste0(.data$ID,"_count"),
                  count_names_TMM = paste0(.data$count_names, "_TMM"))
}


#'create a data.frame with the mean read count / mean TMM read count per pos per group.
#'
#'@param design A data frame with the design matrix of the fitted model
#'@param coef A character string specifying the coeficient of interest - need to be a col of the design df.
#'@param count_table A table of count, with colnames as "IDs_count".
#'@return A data.frame with the mean read count per group (coef == 1) for count and count_TMM.

.compute_mean_read_count <- function(design, coef, count_table){

  count_matrix <- extract_count_matrix(count_table, which = "count")
  count_matrix_TMM <- extract_count_matrix(count_table, which = "count_TMM")
  colnames(count_matrix_TMM) <- stringr::str_remove(colnames(count_matrix_TMM), "_TMM") ## need to change colnmaes to match rownames in design matrix

  IDs_coef <- rownames(design)[design[, coef] == 1]

  mean_count_group1 <- rowMeans(count_matrix[, (colnames(count_matrix) %in% IDs_coef)])
  mean_count_TMM_group1 <- rowMeans(count_matrix_TMM[, (colnames(count_matrix_TMM) %in% IDs_coef)])
  mean_count_group2 <- rowMeans(count_matrix[, !(colnames(count_matrix) %in% IDs_coef)])
  mean_count_TMM_group2 <- rowMeans(count_matrix_TMM[, !(colnames(count_matrix_TMM) %in% IDs_coef)])

  cbind(mean_count_group1, mean_count_group2, mean_count_TMM_group1, mean_count_TMM_group2)

}

#' Help function to find the possible coef given a certain formula.
#'
#'@param formula formula on which to fit the data. The variable used in the
#'formula should be present in the covariate_table.
#'@param covariate_table A data frame with the ID col. and other covariate.
#'@return A vector with coef names.
#'
#'@export
find_coef <- function(covariate_table, formula){

  design <- stats::model.matrix(object = stats::as.formula(formula), data = covariate_table)
  colnames(design)
}

#'Fit a negative Binomial model per position.
#'
#'@param count_table A table of count, with colnames as "IDs_count".
#' @param covariate_table A data frame with the ID col. and other covariate.
#'@param formula A formula on which to fit the data. The variable used in the
#'@param method  Specification of the method to use to fit the data quasi likelihood = "QL or GLM
#'@return A fitted edgeR model.
#'
#'@export
fit_model <- function(count_table, covariate_table, formula, method = "GLM"){

  count_matrix <- extract_count_matrix(count_table, which = "count")

  y <- edgeR::DGEList(counts = count_matrix)
  y <- edgeR::calcNormFactors(y, method="TMM") # calculate the effective library size

  design <- stats::model.matrix(object = stats::as.formula(formula), data = covariate_table)
  rownames(design) <- colnames(y)

  #  fitting the model pvalue from lrt
  y <- edgeR::estimateDisp(y, design)
  if(method == "GLM"){
    edgeR::glmFit(y, design)
  } else if (method == "QL") {
    edgeR::glmQLFit(y, design)
  } else {
    stop("method can only be QL or GLM")
  }
}

#' Reshape the outcome of the edgeR fit.
#'
#'code need to be
#'@param fit A edgeR-fit object.
#'@param coef A coeficient of interest on which the statistics will be computed.
#'@param method GLM or PQL
#'@return A table with the pvalue, fdr, and fitting results for the given coef.
#'
#'@export
reshape_edge_fit <- function(fit, coef, method ="GLM"){

  if(method == "GLM"){
    lrt <- edgeR::glmLRT(fit, coef = coef) ## can be change if needed especially if fitted with QL
  result_table <- edgeR::topTags(lrt, n= nrow(lrt$table), sort.by = "none")[["table"]] ## important to sort by none otherwise change the order

  result_table <-result_table %>%
    dplyr::mutate(pos = rownames(result_table)) %>%
    tidyr::separate(.data$pos, into = c("chr", "start", "stop")) %>%
    dplyr::select(.data$chr, .data$start, .data$stop, dplyr::everything(), -.data$logCPM, -.data$LR) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(coef = coef)

  } else if(method == "QL"){

    pql <- edgeR::glmQLFTest(fit, coef)
    result_table <- edgeR::topTags(pql, n= nrow(pql$table), sort.by = "none")[["table"]] ## important to sort by none otherwise change the order

    result_table <- result_table %>%
      dplyr::mutate(pos = rownames(result_table)) %>%
      tidyr::separate(.data$pos, into = c("chr", "start", "stop")) %>%
      dplyr::select(.data$chr, .data$start, .data$stop, dplyr::everything(), -.data$logCPM, -.data$F) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(coef = coef)
  } else {
    stop("method need to be QL or GLM")
  }
  result_table
}


#'Fit a negative Binomial model per position.
#'
#' @param count_table A table of count, with colnames as "IDs_count".
#' @param covariate_table A data frame with the ID col. and other covariate.
#' @param formula A formula on which to fit the data. The variable used in the
#' formula should be present in the covariate_table.
#' @param coef A coeficient of interest on which the statistics will be computed.
#' @param mean_read logical to specify if the mean read count per group should be computed. FALSE for num or fit with multiple covariate.
#' @param method Specification of the method to use to fit the data quasi likelihood = "QL or GLM.  Default GLM
#' @return A list with the fitted objects and a data.frame with the results.
#' @examples
#' data(full_300_ws_tables_1000)
#' count_table <- full_300_ws_tables_1000$count_table
#' covariate_table <- full_300_ws_tables_1000$covariate_table
#' fit_edgeR(count_table, covariate_table, formula = " ~ group + age_group",
#'           coef = "groupLR", mean_read = FALSE, method = "GLM")
#' @export
fit_edgeR <- function(count_table, covariate_table, formula = "~ group + age_group", coef = "groupLR", mean_read = FALSE, method = "GLM") {

  fit <- fit_model(count_table, covariate_table, formula, method = method)

  fit_table <- reshape_edge_fit(fit, coef, method = method)

  # add the mean read count for the coeficient of interest and the rest.
  if(mean_read){
    mean_count <- .compute_mean_read_count(design = fit$design, coef = coef, count_table =  count_table)

    out_table <- cbind(fit_table, mean_count) %>%
      tibble::as_tibble()

  } else {
    out_table <- fit_table %>%
      tibble::as_tibble()
  }


  list(result_table = out_table,
       fit = fit)

}


