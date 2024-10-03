#' My file logger based on log4r
#'
#' @param x First used for log file and after that for log message
#' @param y Log level
#'
#' @return No return
#'
#' @export
mylogger <- function(x, y = info) {
  stopifnot(!missing(x))
  if (!rlang::env_has(rlang::global_env(), "logger")) {
    my_file_appender = log4r::file_appender(x, append = TRUE, layout = log4r::default_log_layout())
    my_logger <- log4r::logger(threshold = "INFO", appenders = list(my_file_appender))
    rlang::env_poke(rlang::global_env(), "logger", my_logger)
    return()
  }

  my_logger <- rlang::global_env()$logger
  arg <- deparse(substitute(y))

  if (arg == "info")
    log4r::info(my_logger, x)
  else if (arg == "warn")
    log4r::warn(my_logger, x)
  else if (arg == "debug")
    log4r::debug(my_logger, x)
  else if (arg == "error")
    log4r::error(my_logger, x)
  else
    stop("Invalid argument")
}
