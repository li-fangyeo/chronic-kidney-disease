#'
#'
#' @param
#'
#' @return
#'
#' @export
nested_xgboost_cox_cv <- function(data, index = 1, n_inner_folds = 5, n_outer_folds = 5, ...) {
  nested_xgboost_cox_cv_outer_loop(data = data,
                                   iteration = index,
                                   n_outer_folds = n_outer_folds,
                                   n_inner_folds = n_inner_folds,
                                   ...)
}
