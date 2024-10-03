#'
#'
#' @param
#'
#' @return
#'
#' @export
my_surv <- function(x, y) {
  if (is.factor(y)) {
    y <- y %>% as.character() %>% as.numeric()
  }
  survival::Surv(x, y)
}
