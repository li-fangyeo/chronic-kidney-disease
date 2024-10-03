#' Round numeric columns in tibble
#'
#' @param
#'
#' @return
#'
#' @export
round_numeric_columns <- function(x, fmt = '%.4g') {
  dplyr::mutate(x, dplyr::across(where(is.numeric), sprintf, fmt = fmt))
}
