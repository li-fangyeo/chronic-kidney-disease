#'
#'
#' @param
#'
#' @return
#'
#' @export
nested_xgboost_cox_cv_load <- function(x) {
  readRDS(x) %>%
    purrr::map(~purrr::list_modify(.x, matrix = nested_xgboost_cox_cv_model_matrix(.x$data, .x$event_incident, .x$event_agediff, .x$vars)$dmatrix)) %>%
    purrr::map(~purrr::list_modify(.x, model = xgboost::xgb.load.raw(.x$model, as_booster = TRUE)))
}

