#'  Merge the information about the classes (age / rank) and the TMM count numbers
#'
#'  for the specific positions.
#'
#' @param cov_table covariate tables containing group and age_group.
#' @param count_table the count_tables for the specific ID and indexes
#' @param index the indexes of the DMR position to keep track in later steps.
#' @param covariate_from_formula the covariate used to fit the model to pass as covariate for the RF
#'
#' @return A tibble with CPM and covariate from DMRs postion
#' @export

pivote_and_trim_table <- function(cov_table, count_table, index,
                                  covariate_from_formula = c("group", "age_group")){

  x <- as.data.frame(t(as.matrix(count_table)))

  t <- dplyr::bind_cols(cov_table[, covariate_from_formula], x) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.character), as.factor))
  colnames(t) <- c(covariate_from_formula, index)
  t
}

################################################################################


#' Extract the indexes of the DMR positions. (needed to get the correct CPM in the count_table)
#'
#' @param table_with_DMR the DMR_train col of the set.
#' @param DMR_to_use The coef to use (either groupLR, or age_groupyoung)
#' @param  q_cut The cut-off value on the FDR value.
#'
#' @return A vector of indexes
#' @export

extract_indexes <- function(table_with_DMR, DMR_to_use, q_cut){
  table_with_DMR %>%
    dplyr::filter(.data$coef %in% DMR_to_use & .data$FDR <= q_cut) %>%
    dplyr::pull(.data$index)
}


################################################################################


#' Create the input table for the RF models. Combine information about the
#' samples and the CPM at the DMR positions.
#'
#' @param table_to_turn either the train_data or the test_data
#' @param count_table_full the full count table in which the CPM can be extracted.
#' @param index  the indexes of the DMR position to keep track in later steps.
#' @param covariate_from_formula the covariate used to fit the model to pass as covariate for the RF
#'
#' @return A tibble trimmed for the index position and for the ID to predict
#' @export


reshape_tables_for_RF <- function(table_to_turn, count_table_full,
                                  index,
                                  covariate_from_formula = c("group", "age_group")){
  new_count <- count_table_full[index, table_to_turn[["count_names_TMM"]]] ## indexes  select the correct row table_to_turn[count_names..] for the correct IDs.
  pivote_and_trim_table(table_to_turn,
                        new_count,
                        index,
                        covariate_from_formula)
}

###########################################################################


#'  Create the input table for the RF models. Combine information about the
#'  samples and the CPM at the DMR positions.
#'
#' @param set The full set including the col: data_test, data_train; DMR_train.
#' @param count_table_full the full count table in which the CPM can be extracted.
#' @param q_cut  The cut-off value on the FDR scale.
#' @param variable_to_test The variable to be classified. Either group or age_group.
#' @param DMR_to_use The coef of the DMR used to predict the variable_to_test.
#' @param index vector of index for targeted DMRs. Default to null, the indexes will be computed from the DMR_to_use and the q_cut
#' @param covariate_from_formula the covariate used to fit the model to pass as covariate for the RF
#'
#' @return A tibble with the summary results of the fitted RF.
#' @export

run_rf_on_set <- function(set,
                          count_table_full,
                          variable_to_test,
                          q_cut,
                          DMR_to_use,
                          index = NULL,
                          covariate_from_formula = c("group", "age_group")
){

  out <- list()

  for(i in seq_along(set$data_test)){

    ## extract indexes for the good threshold and the good variable
    if(is.null(index)){
      index <- extract_indexes(set$DMR_train[[i]], DMR_to_use, q_cut) ## will take the index of the DMR in each subsets
    } else {
      index <- index
    }

    if(length(index) < 1) stop("no DMR for the given covariate and q_cut, impossible to fit a RF")



    train_rf_data <- reshape_tables_for_RF(set$data_train[[i]], count_table_full, index, covariate_from_formula)
    test_rf_data <-  reshape_tables_for_RF(set$data_test[[i]], count_table_full, index, covariate_from_formula)


    fit <-   ranger::ranger(dependent.variable.name = variable_to_test, ## fit rf on the train set
                            data = train_rf_data,
                            importance = "impurity",
                            save.memory = T,
                            num.threads = 10, ## paralelll setting
                            num.trees = 10000)

    pred <- as.factor(stats::predict(fit, test_rf_data)$prediction) ## predict on the test set

    out[[i]] <- tibble::tibble(ID = set$data_test[[i]][["ID"]],
                               observed =  set$data_test[[i]][[variable_to_test]],
                               pred= pred,
                               set = set[["set_ID"]][i],
                               predicted_var = variable_to_test,
                               DMR_used = paste(DMR_to_use, collapse = "/"),
                               oob = fit$prediction.error, ## extract the oob error from the fitted model
                               importance = list(fit$variable.importance))

  }
  purrr::map_dfr(out, ~.x)
}

################################################################################


#' Provide summary statistic on the performance of the classifier.
#' return the accuracy, the f1-score, the younden J-statistic, the recall and precision rates.
#'
#' @param rf_outcome_table the output table of the rf fitted on the cv set with run_rf_on_set()
#'
#' @return A tibble with the summary statistics for the predictive power of the RF
#' @export
summarize_cv <- function(rf_outcome_table){


  group_summary <- rf_outcome_table %>%
    dplyr::mutate(observed = as.factor(.data$observed)) %>%
    dplyr::group_by(.data$set) %>%
    dplyr::summarise(n = dplyr::n(),
                     oob = .data$oob[1],
                     tp = sum(.data$pred == levels(.data$pred)[1]  & .data$observed == levels(.data$observed)[1]),
                     tn = sum(.data$pred == levels(.data$pred)[2]  & .data$observed == levels(.data$observed)[2]),
                     fp = sum(.data$pred == levels(.data$pred)[1]  & .data$observed == levels(.data$observed)[2]),
                     fn = sum(.data$pred == levels(.data$pred)[2]  & .data$observed == levels(.data$observed)[1]),
                     success = levels(.data$pred)[1])

  group_summary %>% dplyr::mutate(accuracy = as.numeric(cutpointr::accuracy(tp = .data$tp, fp = .data$fp, tn = .data$tn, fn = .data$fn)),
                           younden = as.numeric(cutpointr::youden(tp = .data$tp, fp = .data$fp, tn = .data$tn, fn = .data$fn)),
                           f1 = as.numeric(cutpointr::F1_score(tp = .data$tp, fp = .data$fp, tn = .data$tn, fn = .data$fn)),
                           recall = as.numeric(cutpointr::recall(tp = .data$tp, fp = .data$fp, tn = .data$tn, fn = .data$fn)),
                           precision = as.numeric(cutpointr::precision(tp = .data$tp, fp = .data$fp, tn = .data$tn, fn = .data$fn)))

}

################################################################################


#' create a list with indexes for DMRs -robust, non-robust, full, and randoms...
#'
#'
#' @param DMR_group tables with DMR and robustness info.
#' @param tables list with read count table and cov table
#'
#' @return A list with the indexes of coherent and non coherent DMR + libsize.
#' @export
.creat_indexes_list <- function(DMR_group, tables){

  index_full <- unique(DMR_group$index)
  index_robust <- unique(DMR_group$index[DMR_group$coherent_block])
  index_non_robust <- unique(DMR_group$index[!DMR_group$coherent_block])
  ## add the random
  index_random_full <- round(stats::runif(length(index_full), 1, nrow(tables$count_table)), 0)
  index_random_robust <- round(stats::runif(length(index_robust), 1, nrow(tables$count_table)), 0)
  index_random_non_robust <- round(stats::runif(length(index_non_robust), 1, nrow(tables$count_table)), 0)
  index_random_libsize <- round(stats::runif(1, 1, nrow(tables$count_table)), 0)## 1 positivon only for libsize.

  list(index_full = index_full,
       index_robust = index_robust,
       index_non_robust = index_non_robust,
       index_random_full = index_random_full,
       index_random_robust = index_random_robust,
       index_random_non_robust = index_random_non_robust,
       index_random_libsize = index_random_libsize)
}


#' create a list with indexes for DMRs -robust, non-robust, full, and randoms...
#'
#'
#' @param tables the full_300_ws_tables
#'
#' @return A list with the indexes of coherent and non coherent DMR + libsize.
#' @export
.extract_standardized_lib_size <-function(tables) {

  tables$covariate_table %>%
    dplyr::select(.data$group,  .data$lib.size, .data$norm.factors, .data$age_group) %>%
    dplyr::mutate(lib.size = .data$lib.size*.data$norm.factors) %>%
    dplyr::select(-.data$norm.factors) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.character), as.factor))

}

#' Fit a random forest model trained only on the ~ mean methylation.
#'
#'
#' @param tables the full_300_ws_tables
#'
#' @return A tibble with oob error and other variables
#' @examples
#' \dontrun{
#' data(full_300_ws_tables_1000)
#' rf_on_lib_size(full_300_ws_tables_1000)
#' }
#' @export
rf_on_lib_size <- function(tables){
  cov_group <- .extract_standardized_lib_size(tables)

  test <- ranger::ranger(dependent.variable.name = "group",
                         data = cov_group,
                         num.trees = 10000,
                         save.memory = T, num.threads = 10, importance = "impurity")## same param as for DMRs

  tibble::tibble(set = as.factor(1),
                 oob = test$prediction.error,
                 n = 42,
                 DMR = "lib_size", success = "HR") ## same output format to simplify integration

}

#' Fit random forest model trained on DMR and mean meth.
#'
#' run on the list of indexes and return a tibble
#'
#' @param tables the full_300_ws_tables
#' @param DMR_group table with DMRs for group variable
#' @param set the nested table containing the fitted model
#' @return A tibble with oob error and other variables
#' @examples
#' \dontrun{
#' load("data-raw/data/full_300_ws_tables.rda")
#' data(DMR_group)
#' data(model_1_one_fold)
#' fit_and_summarise_RF(DMR_group, model_1_one_fold, full_300_ws_tables)
#' }
#' @export
fit_and_summarise_RF <- function(DMR_group, set, tables){

  index_list <- .creat_indexes_list(DMR_group, tables)
  out <- list()
  for(i in seq_along(index_list)) {
    out[[i]] <-  run_rf_on_set(set = set,
                               count_table_full = tables$count_table,
                               variable_to_test = "group",
                               q_cut = 1,
                               index = index_list[[i]],
                               DMR_to_use = "groupLR") %>%
      summarize_cv() %>%
      dplyr::mutate(DMR = paste(names(index_list)[i])) %>%
      dplyr::select(.data$set, .data$oob, .data$n, .data$DMR, .data$success)
  }
  out_lib <- rf_on_lib_size(tables) ## add lib  size

  dplyr::bind_rows(out) %>% dplyr::bind_rows(out_lib)
}

#' Fit random forest model trained on DMR and mean meth.
#'
#' Bootstrap X times the fit_and_summarise_RF to have a less variable
#' estimate of the oob error. Especially for the randomise run.
#'
#' @inheritParams fit_and_summarise_RF
#' @param n number of bootstrap
#' @return A tibble with the mean oob error accross bootstrap

#' @export
estimate_oob <- function(DMR_group, set, tables, n = 1000){
  out <- list()
  for(i in 1:n){
    out[[i]] <- fit_and_summarise_RF(DMR_group, set, tables)
  }
  dplyr::bind_rows(out) %>%
    dplyr::group_by(.data$DMR) %>%
    dplyr::summarise(mean_oob = mean(.data$oob))
}

