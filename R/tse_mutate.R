#' TSE mutate
#'
#' @param
#'
#' @return
#'
#' @export
tse_mutate <- function(.data, ...) {
  dots <- rlang::enquos(...)
  new_data <-  tse_meta(.data) %>%
    dplyr::mutate(!!!dots) %>%
    tibble::column_to_rownames(var = "rownames") %>%
    S4Vectors::DataFrame()
  SummarizedExperiment::colData(.data) <- new_data
  .data
}
