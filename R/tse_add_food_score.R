#' Add Koponen's Healthy Food Choices (HFC) score to the data frame. HFC is defined in DOI:10.1093/ajcn/nqab077.
#'
#' @param .data TSE object
#' @param col Column to store food score
#'
#' @return TSE with additional coldata column
#'
#' @export
tse_add_food_score <- function(.data, col) {
  stopifnot(!missing(.data), !missing(col))
  hfc <- function(x, high) { c(0.5, 1.5, 4.3, 8.6, 21.5, high)[x] }

  new_data <-  tse_meta(.data) %>%
    dplyr::mutate({{col}} := hfc(KY100_14, 60) +
                    hfc(KY100_17, 60) +
                    hfc(KY100_9, 60) +
                    hfc(KY100_1, 60) +
                    hfc(KY100_2, 60) +
                    hfc(KY100_15, 45) +
                    hfc(KY100_16, 45) +
                    hfc(KY100_18, 45) +
                    hfc(KY100_19, 45) +
                    hfc(KY100_20, 45) +
                    hfc(KY100_21, 45) +
                    hfc(FR02_101_3, 45) +
                    hfc(FR02_101_7, 45) +
                    hfc(FR02_101_8, 45)) %>%
    tibble::column_to_rownames(var = "rownames") %>%
    S4Vectors::DataFrame()
  SummarizedExperiment::colData(.data) <- new_data
  .data
}
