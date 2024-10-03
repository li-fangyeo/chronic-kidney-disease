#' TSE meta
#'
#' @param
#'
#' @return
#'
#' @export
tse_meta <- function(x, rownames = TRUE) {
  { if (inherits(x, "TreeSummarizedExperiment")) SummarizedExperiment::colData(x) else x } %>%
    { if(rownames) tibble::as_tibble(., rownames = "rownames") else tibble::as_tibble(.) }
}
