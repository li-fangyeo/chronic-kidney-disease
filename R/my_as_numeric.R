#'
#'
#' @param
#'
#' @return
#'
#' @export
my_as_numeric <- function(x) {
  y <- rep(NA, length(x))
  ids_numeric <- stringr::str_detect(x, "^[0-9]+\\.?[0-9]*$")

  y[ids_numeric] <- as.numeric(x[ids_numeric])

  ids_small <- stringr::str_detect(x, "^<[0-9]+\\.?[0-9]*$")
  y[ids_small] <- 0

  y
}
