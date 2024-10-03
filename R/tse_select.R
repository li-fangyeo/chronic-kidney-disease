#' TSE select
#'
#' @param
#'
#' @return
#'
#' @export
tse_select <- function(.data, ...) {
  dots <- rlang::enquos(...)
  new_data <-  tse_meta(.data) %>%
    dplyr::select(rownames, !!!dots) %>%
    tibble::column_to_rownames(var = "rownames") %>%
    S4Vectors::DataFrame()
  SummarizedExperiment::colData(.data) <- new_data
  .data
}
