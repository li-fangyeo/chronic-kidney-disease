#' TSE filter
#'
#' @param
#'
#' @return
#'
#' @export
tse_filter <- function(.data, ...) {
  dots <- rlang::enquos(...)
  ids <-  tse_meta(.data) %>%
    dplyr::filter(!!!dots) %>%
    dplyr::pull(rownames)
  .data[, ids]
}
