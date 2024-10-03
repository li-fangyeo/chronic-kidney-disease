#'
#'
#' @param
#'
#' @return
#'
#' @export
nested_xgboost_cox_cv_outer_loop <- function(data,
                                             seed,
                                             iteration,
                                             event_incident,
                                             n_random = 6,
                                             n_outer_folds = 5,
                                             ...) {
  set.seed(seed + iteration)

    data_with_groups <- data %>%
    add_random_number_columns(n_random) %>% 
    dplyr::mutate(group = caret::createFolds(factor(data[[event_incident]]), k = n_outer_folds, list = FALSE))

  partial_nested_xgboost_cox_cv_inner_loop <- purrr::partial(nested_xgboost_cox_cv_inner_loop,
                                                             data = data_with_groups,
                                                             event_incident = event_incident,
                                                             ...)

  seq(n_outer_folds) %>%
    purrr::map(~partial_nested_xgboost_cox_cv_inner_loop(.x, seed = seed + iteration + .x))
}
