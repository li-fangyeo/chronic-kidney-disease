#' Function performs single inner fold nested cross-validation step.
#'
#' 1. Split data frame to training and test sets.
#' 2. Perform model parameter tuning in training data using n-fold inner cross-validation
#' 3. Test tuned model using independent test data
#'
#' @param
#'
#' @return
#'
#' @export
nested_xgboost_cox_cv_inner_loop <- function(data,
                                             group,
                                             event_incident,
                                             event_agediff,
                                             vars,
                                             taxa,
                                             randoms,
                                             idcol,
                                             seed = 2024,
                                             booster = "gbtree",
                                             objective = "survival:cox",
                                             eval_metric = "cox-nloglik",
                                             training_fraction = 0.70,
                                             n_inner_folds = 5,
                                             nthread = 10,
                                             n_design = 50,
                                             opt_steps = 150,
                                             early_stopping_rounds = 100) {
  set.seed(seed)

  df_training <- data %>% dplyr::filter(group != {{group}})

  mtx_training <- nested_xgboost_cox_cv_model_matrix(df_training, event_incident, event_agediff, vars, taxa, randoms)

  .pb <- knitrProgressBar::progress_estimated(n_design + opt_steps)

  smoof_fun <- smoof::makeSingleObjectiveFunction(
                        name = "xgb_cv_bayes",
                        fn = function(x) {
                          xgbcv <- xgboost::xgb.cv(params = list(
                                                     booster = booster,
                                                     objective = objective,
                                                     eval_metric = eval_metric,
                                                     eta = x["eta"],
                                                     max_depth = x["max_depth"],
                                                     min_child_weight = x["min_child_weight"],
                                                     gamma = x["gamma"],
                                                     subsample = x["subsample"],
                                                     colsample_bytree = x["colsample_bytree"]),
                                                   data = mtx_training$dmatrix,
                                                   nround = x["nrounds"],
                                                   nfold = n_inner_folds,
                                                   stratified = TRUE,
                                                   early_stopping_rounds = early_stopping_rounds,
                                                   nthread = nthread,
                                                   maximize = FALSE,
                                                   verbose = 0)
                          knitrProgressBar::update_progress(.pb)
                          xgbcv$evaluation_log$test_cox_nloglik_mean %>% min()
                        },
                        par.set = ParamHelpers::makeParamSet(
                                                  ParamHelpers::makeIntegerParam("nrounds", lower = 50, upper = 10000),
                                                  ParamHelpers::makeNumericParam("eta", lower = 0.001, upper = 0.5),
                                                  ParamHelpers::makeNumericParam("gamma", lower = 0.1, upper = 100),
                                                  ParamHelpers::makeIntegerParam("max_depth", lower = 2, upper = 20),
                                                  ParamHelpers::makeIntegerParam("min_child_weight", lower = 1, upper = 200),
                                                  ParamHelpers::makeNumericParam("subsample", lower = 0.2, upper = 0.8),
                                                  ParamHelpers::makeNumericParam("colsample_bytree", lower = 0.2, upper = 0.8)),
                        minimize = TRUE)

  mlr_run <- mlrMBO::mbo(fun = smoof_fun,
                         design = ParamHelpers::generateDesign(n = n_design,
                                                               par.set = ParamHelpers::getParamSet(smoof_fun),
                                                               fun = lhs::randomLHS),
                         control = mlrMBO::setMBOControlTermination(mlrMBO::makeMBOControl(), iters = opt_steps),
                         show.info = TRUE)

  df_testing <- data %>% dplyr::filter(group == {{group}})
  mtx_testing <- nested_xgboost_cox_cv_model_matrix(df_testing, event_incident, event_agediff, vars, taxa, randoms)

  xgb_model <- xgboost::xgboost(params = mlr_run[["x"]] %>%
                                  purrr::list_modify(nrounds = purrr::zap(),
                                                     booster = booster,
                                                     objective = objective,
                                                     eval_metric = eval_metric),
                                data = mtx_testing$dmatrix,
                                nrounds = mlr_run[["x"]]$nrounds,
                                nthread = nthread,
                                early_stopping_rounds = early_stopping_rounds,
                                verbose = 0,
                                maximize = FALSE,
                                print_every_n = 500)

  list(data = df_testing,
       event_incident = event_incident,
       event_agediff = event_agediff,
       ids = df_testing[[idcol]],
       matrix = mtx_testing$dmatrix,
       vars = c(vars, taxa, randoms),
       params = mlr_run[["x"]],
       model = xgb_model)
}
