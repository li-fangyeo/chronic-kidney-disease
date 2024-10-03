#' TSE filter at
#'
#' @param
#'
#' @return
#'
#' @export
tse_filter_at <- function(x, vars, vars_predicate) {
  ids <- tse_meta(x) %>%
    dplyr::filter(dplyr::across({{vars}}, {{vars_predicate}})) %>%
    dplyr::pull(rownames)
  x[, ids]
}
