#' Set difference
#'
#' @param
#'
#' @return
#'
#' @export
`%difference%` <- function(a, b) {
    a[!a %in% b]
}

#' Set intersect
#'
#' @param
#'
#' @return
#'
#' @export
`%intersect%` <- function(a, b) {
    intersect(a, b)
}

#' Set union
#'
#' @param
#'
#' @return
#'
#' @export
`%union%` <- function(a, b) {
    c(a, b)
}

#' Set unique union
#'
#' @param
#'
#' @return
#'
#' @export
`%uniqunion%` <- function(a, b) {
    all <- c(a, b)
    all[!duplicated(all)]
}
