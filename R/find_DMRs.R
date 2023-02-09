################################################################################
################################################################################
#### This script contains the functions used to create the sets used for the
#### DMR analysis. Both for theCross validation procedure and for the full set.
#### The differents tables are organised as list-col in a data frame. The df
#### has 4 column:
####
#### * set_ID: the ID of the set
#### * data_test: the covariate table including ONLY the IDs used for the train set
#### * data_train: the covariate table including the IDs NOT used for the train set
#### * DMR_train: the DMR observed in the train set. given a specific formula and p_cut.
####
#### The data_test and data_train for the full_set are equal and contains all IDs.
#### The table created with these functions represents the base table for the RF
#### input (in conjunction with the count_tables) and for subsequent analysis of the DMR
#### based on windows.
#### last update: 12.12.22
####
################################################################################
################################################################################


#' Create a data frame with list-col with data_train and data_test.
#'
#' Keeping the variable of interest (group) balanced.
#'
#' @param k the number of folds
#' @param tables the covariate and count tables created with prepare_ws_datasets()
#' @param variable_to_balance the variable name that need to be present in both sets.
#' @return A nested tibble with K row
#' @examples
#' data(full_300_ws_tables_1000)
#' prepare_balanced_k_set(6, full_300_ws_tables_1000, "group")
#'
#'@export
prepare_balanced_k_set <- function(k, tables, variable_to_balance){

  tables$covariate_table$set_ID <- groupdata2::fold(data = tables$covariate_table, cat_col = variable_to_balance, k = k)[[".folds"]]

  set_cov <- tables$covariate_table %>%
    dplyr::group_by(.data$set_ID) %>%
    tidyr::nest()

  set_cov$data_train <- purrr::map(set_cov$set_ID, ~ tibble::as_tibble(tables$covariate_table[tables$covariate_table$set_ID != .x, ]))
  set_cov <- set_cov %>% dplyr::rename(data_test = .data$data)

  if(k == 1){

    set_cov$data_train <- set_cov$data_test
    set_cov$data_train[[1]]$set_ID <- 1 ## just add the set_ID col to be sync with the 6 fold files

  }
  set_cov
}

################################################################################

#' Filter the count_table to keep ONLY the count and count_TMM for the IDs in the train set.
#'
#'
#' @param train_set the talbe with the IDs in the train set.
#' @param tables the covariate and count tables created with prepare_ws_datasets()
#' @return A nested tibble with K row
#'
#' @export

reduce_count_table <- function(train_set, tables){

  IDs <- train_set$count_names
  IDs2 <- train_set$count_names_TMM

  count_train <- tables$count_table %>%
    dplyr::select(1:4, !!IDs, !!IDs2)

  list(covariate_table = train_set,
       count_table = count_train)

}

################################################################################

#' Fit a negbin model using EdgeR to detect DMR.
#'
#' nb. It uses fit_edgeR that is part of the epiena package developped to work with this datasets.
#' and is a wrapper around the edgeR functions.
#'
#' @param tables the covariate and count tables created with prepare_ws_datasets()
#' @param formula the formula for the model
#' @param coef the coeficient to extract the DMR for (can be multiple).
#' @param p_cut the p_value threshold on which to keep the DMR. (mostly used to keep object size decent)
#' @param index a vector of intex to retain the value of the fit. THought for the CV to keep the value of the DMR in the full set
#' @param method GLM or QL
#'
#' @return A tibble with DMRs
#' @export
get_DMRs <- function(tables, formula, coef, p_cut, index = NULL,  method = "GLM"){

  fit <- fit_edgeR(count_table = tables$count_table,
                   covariate_table = tables$covariate_table,
                   formula = formula,
                   coef = coef[1], ## run first on the first coef
                   mean_read = F, ## hard coded as we might change the formula can be passed into an arg. if necessary
                   method = method)

  fit$result_table$index <- 1:nrow(fit$result_table)

  out <- list()
  if(is.null(index)){
    out[[1]] <- fit$result_table[fit$result_table$PValue <= p_cut, ]
  } else {
    out[[1]] <- fit$result_table[fit$result_table$PValue <= p_cut | fit$result_table$index %in% index, ] ## allow to retain value for specific position.
  }

  if(length(coef) > 1){ ## if more than 1 coef loop on the second forward. Add rows.

    for(i in 2:length(coef)){

      test <-reshape_edge_fit(fit$fit, coef = coef[i])
      test$index <- 1:nrow(test)
      if(is.null(index)){
        out[[i]] <- test[test$PValue <= p_cut , ]
      } else {
        out[[i]] <- test[test$PValue <= p_cut | test$index %in% index, ]
      }
    }
  }
  do.call(what = rbind, out)
}


################################################################################


#' Wrapper around the get_DMRs() function. It first cut the count_table according to the data_train.
#' and then extract the DMR with get_DMRs()
#' @param train_set the talbe with the IDs in the train set.
#' @param tables the covariate and count tables created with prepare_ws_datasets()
#' @param formula the formula for the model
#' @param coef the coeficient to extract the DMR for (can be multiple).
#' @param p_cut the p_value threshold on which to keep the DMR. (mostly used to keep object size decent)
#' @param index a vector of intex to retain the value of the fit. THought for the CV to keep the value of the DMR in the full set
#' @param method GLM or QL
#'
#' @return A DF with DMRs
#' @export
get_DMRs_from_train_set <- function(train_set, tables, formula,  p_cut, coef, index = NULL, method = "GLM"){

  tables <- reduce_count_table(train_set = train_set, tables)
  get_DMRs(tables = tables, formula = formula, coef = coef, p_cut = p_cut, index = index, method = method)


}


################################################################################

#'
#' Add the DMRs as a list_col to the set table.
#' nb. the coef are hard coded for the current usage but can easily be changed.
#'
#' @param set the table created with prepare_balanced_k_set().
#' @param full_tables the covariate and count tables created with prepare_ws_datasets()
#' @param formula the formula for the model
#' @param p_cut the p_value threshold on which to keep the DMR. (mostly used to keep object size decent)
#' @param index a vector of intex to retain the value of the fit. THought for the CV to keep the value of the DMR in the full set
#' @param method GLM or QL
#'
#' @return A nested DF with DMRs for each subset
#' @examples
#' data(full_300_ws_tables_1000)
#' set <- prepare_balanced_k_set(6, full_300_ws_tables_1000, "group")
#' set2 <- add_DMRs_to_set(set, full_300_ws_tables_1000, "~group + age_group", p_cut = 1)
#' @export
add_DMRs_to_set <- function(set, full_tables, formula, p_cut, index = NULL, method = "GLM"){

  set$DMR_train <- purrr::map(set$data_train,
                              get_DMRs_from_train_set,
                              tables = full_tables,
                              formula = formula,
                              coef = find_coef(full_tables$covariate_table, formula = formula)[2:length(find_coef(full_tables$covariate_table, formula = formula))], ### take all coef except intercept
                              p_cut = p_cut,
                              index = index,
                              method = method)
  set
}

################################################################################

#' extract the DMRs for the given covariate and add the presence of CV.
#'
#'
#' @param set the table created with prepare_balanced_k_set().
#' @param coef the coeficient on which to extract the DMRs.
#' @param p_cut the p_value threshold on which to keep the DMR. (if filter on q_value set to 1)
#' @param q_cut the FDR threshold on which to keep the DMR. (if filter on pvalue set to 1)
#'
#' @return A nested DF with DMRs for each subset
#' @examples
#' data(model_1_one_fold)
#' reshape_DMR_table(model_1_one_fold, coef = "groupLR", p_cut = 1, q_cut = 0.05)
#' @export
reshape_DMR_table <- function(set, coef, p_cut, q_cut){

  coef2 <- coef

  DMR_table <- dplyr::bind_rows(set$DMR_train) %>%
    dplyr::filter(.data$FDR <= q_cut & .data$PValue <= p_cut & .data$coef == coef2)

  DMR_table %>%
    dplyr::group_by(.data$index) %>%
    dplyr::mutate(presence = dplyr::n()) %>%
    dplyr::distinct(.data$index, .keep_all = TRUE) %>%
    dplyr::mutate(direction = ifelse(.data$logFC > 0, "hyper", "hypo")) ## Be aware will only take the value of the first. might differe between sets

}
