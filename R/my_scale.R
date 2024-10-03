#' Scale to zero mean and unit variance
#'
#' @param list of numeric values
#'
#' @return scaled list
#'
#' @export
  my_scale <- function(x) {
        (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
  }
