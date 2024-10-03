#'
#'
#' @param
#'
#' @return
#'
#' @export
nested_xgboost_cox_cv_harrell <- function(x) {
  pred <- predict(x$model, x$matrix)
  time <- xgboost::getinfo(x$matrix, "label")
  index <- ifelse(time > 0, 1, 0)
  intsurv::cIndex(abs(time), event = index, pred)
}
