#' Add random number columns in data frame with specific seed to validate machine learning algorightms
#'
#' @param
#'
#' @return
#'
#' @export
add_random_number_columns <- function(df, n) {
  stringr::str_c("rnorm_", seq_len(n)) %>%
    purrr::reduce(~tibble::add_column(.x, '{.y}' := rnorm(n = nrow(df), mean = 0, sd = 1)), .init = df)
}
