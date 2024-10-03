#'
#'
#' @param
#'
#' @return
#'
#' @export
nested_xgboost_cox_cv_model_matrix <- function(x, event_incident, event_agediff, ...) {
    dset <- x %>%
      dplyr::filter(!is.na(!!rlang::sym(event_incident))) %>%
      dplyr::mutate(event = ifelse(!!rlang::sym(event_incident) == 1,
                                   !!rlang::sym(event_agediff),
                                   -!!rlang::sym(event_agediff)))
    labels <- dset %>% dplyr::pull(event)
    matrix <- dset %>% dplyr::select(dplyr::one_of(c(...))) %>% data.matrix()
    dmatrix <- xgboost::xgb.DMatrix(matrix, label = labels)
    list(dmatrix = dmatrix, matrix = matrix, labels = labels)
  }
