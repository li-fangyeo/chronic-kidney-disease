#'
#'
#' @param
#'
#' @return
#'
#' @export
compute_or_load_result <- function(computation_function, file_path) {
  if (file.exists(file_path)) {
    result <- readRDS(file_path)
  } else {
    result <- computation_function()
    saveRDS(result, file_path)
  }
  result
}
