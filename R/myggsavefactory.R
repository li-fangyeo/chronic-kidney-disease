#' my ggsave factory
#'
#' @param
#'
#' @return
#'
#' @export
myggsavefactory <- function(filepath = "cache", time = Sys.time()) {
  datetime <- format(time, "%Y%m%d_%H%M%S")
  dir.create(filepath, showWarnings = FALSE)

  function(plot, name = NULL, draw = TRUE, extension = "png", ...) {
    if (is.null(name)) name <- deparse(substitute(plot))
    filename_with_time <- stringr::str_glue("{filepath}/{name}-{datetime}.{extension}")
    ggplot2::ggsave(filename = filename_with_time, plot = plot, ..., limitsize = FALSE)
    if (knitr::is_html_output() && draw)
      knitr::include_graphics(normalizePath(filename_with_time), rel_path = FALSE)
    else
      invisible(filename_with_time)
  }
}
