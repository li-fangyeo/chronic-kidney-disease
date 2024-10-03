#'
#'
#' @param
#'
#' @return
#'
#' @export
nested_xgboost_cox_cv_save <- function(x, saver = mysavefactory(), name = "nested_xgboost_cox_cv") {
  x %>%
    purrr::map(~purrr::list_modify(.x, matrix = purrr::zap())) %>%
    purrr::map(~purrr::list_modify(.x, model = xgboost::xgb.save.raw(.x$model, raw_format = "json"))) %>%
    saver(name = name)
}
