#' Function factory for saving objects to specific directory with constant time code
#'
#' @param dir: target directory
#'
#' @return save function
#'
#' @export
mysavefactory <- function(objectdir = "rds", plotdir = "cache", time = Sys.time()) {
  datetime <- format(time, '%Y%m%d_%H%M%S')
  dir.create(objectdir, showWarnings = FALSE)
  dir.create(plotdir, showWarnings = FALSE)

  objectsaver <- function(object, name, extension = "rds") {
    filename_with_time <- stringr::str_glue("{objectdir}/{name}-{datetime}.{extension}")
    saveRDS(object, filename_with_time)
    invisible(filename_with_time)
  }

  plotsaver <- function(object, name, draw = TRUE, extension = "png", ...) {
    filename_with_time <- stringr::str_glue("{plotdir}/{name}-{datetime}.{extension}")
    ggplot2::ggsave(filename = filename_with_time, plot = object, ..., limitsize = FALSE)
    if (knitr::is_html_output() && draw)
      knitr::include_graphics(normalizePath(filename_with_time), rel_path = FALSE)
    else
      invisible(filename_with_time)
  }

  function(object, name = NULL, ...) {
    if (is.null(name)) name <- deparse(substitute(object))
    if (inherits(object, "ggplot")) plotsaver(object, name, ...)
    else objectsaver(object, name, ...)
  }
}
