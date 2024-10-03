#'
#'
#' @param
#'
#' @return
#'
#' @export
my_adjust_p <- function(..., n = length(c(...))) {
  pmin(1, c(...)*n)
}
