#' TSE filter
#'
#' @param
#'
#' @return
#'
#' @export
tse_if_all <- function(.data, ...) {
  dots <- rlang::enquos(...)
  ids <-  tse_meta(.data) %>%
    dplyr::if_all(!!!dots) %>%
    dplyr::pull(rownames)
  .data[, ids]
}
