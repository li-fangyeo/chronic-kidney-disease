

#31052024
#Nephro-metagenomics in FINRISK 2002
library(magrittr)
devtools::load_all()

#Arguments
#makes it easier to run scripts (.sh) applicable when using emacs
args <- list(
  optparse::make_option("--east", action="store_true", default=FALSE, help="Exclude Eastern Finland subpopulation [default \"%default\"]"),
  optparse::make_option("--west", action="store_true", default=FALSE, help="Exclude Eastern Finland subpopulation [default \"%default\"]"),
  optparse::make_option("--detection", type = "numeric", default = 0.1/100, help = "Detection limit [default %default]"),
  optparse::make_option("--prevalence", type = "numeric", default = 1/100, help = "Prevalence limit [default %default]")) %>% 
  optparse::OptionParser(option_list = .) %>%
  optparse::parse_args()

args %>% tibble::enframe(name = "Option", value = "Argument") %>% DT::datatable()

#Formatting how outputs are saved
myggsave <- myggsavefactory()

#Formatting image outputs
{ ggthemes::theme_tufte(base_family = "sans", base_size = 12) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "black", fill = NA),
      axis.text = ggplot2::element_text(colour = "black", size = 10),
      aspect.ratio = 1) } %>%
  ggplot2::theme_set()

#Data definition
vars <- list(BL_AGE = "Age",
             MEN = "Men",
             BMI = "BMI",
             PREVAL_DIAB.col_from_endpoints = "Diabetes",
             SYSTM = "Systolic blood pressure",
             BP_TREAT = "Antihypertensive medication",
             PREVAL_HFAIL_STRICT.col_from_endpoints = "Heart failure",
             PREVAL_CKD = "Prevalent CKD",
             INCIDENT_CKD = "Incident CKD",
             CKD_AGEDIFF = "CKD Agediff",
             CURR_SMOKE = "Smoking",
             PREVAL_AUTOIMMUN.col_from_endpoints = "Autoimmune disease",
             KREA_ENTS = "Creatinine",
             GFR = "Glomerulal filtration rate",
             UAlbKrea = "Urine Albumin-Creatinine Ratio",
             EAST = "Eastern Finland",
             shannon = "Shannon diversity") 

#Read in data, calculate alpha, tidy data for selected variables
tse <- readRDS("data/tse_gg2_MGS_FR02.rds") %>%
  mia::transformAssay(assay.type = "counts", method = "relabundance") %>% 
  mia::estimateDiversity(assay.type = "counts", index = "shannon", name = "shannon") %>%
  tse_mutate(dplyr::across(c(MEN,
                             EAST,
                             BP_TREAT,
                             CURR_SMOKE,
                             dplyr::contains("INCIDENT"),
                             dplyr::contains("PREVAL")), as.factor)) %>% 
  tse_mutate(GFR = 0.993^round(BL_AGE) *
               dplyr::case_when(MEN == 0 & KREA_ENTS <= 62 ~ 144*(KREA_ENTS/61.9)^-0.329,
                                MEN == 0 & KREA_ENTS >  62 ~ 144*(KREA_ENTS/61.9)^-1.209,
                                MEN == 1 & KREA_ENTS <= 80 ~ 141*(KREA_ENTS/79.6)^-0.411,
                                MEN == 1 & KREA_ENTS >  80 ~ 141*(KREA_ENTS/79.6)^-1.209)) %>%
  tse_mutate(UAlbKrea = U_ALB/U_KREA) %>% 
  tse_filter(GRAVID %in% c(1, NA), BL_USE_RX_J01_1mo %in% c(0, NA)) %>% 
  tse_filter(GFR >= 60, PREVAL_CKD == 0, UAlbKrea <=3 | is.na(UAlbKrea)) %>%
  { if (args$east) tse_filter(., EAST == 0) else . } %>% 
  { if (args$west) tse_filter(., EAST == 1) else . } %>% 
  tse_filter(dplyr::if_all(dplyr::one_of(names(vars) %difference% "UAlbKrea"), not_na)) %>%
  tse_select(names(vars))

#Table one / characteristics of data
tse %>%
  tse_meta(rownames = FALSE) %>%
  mytableone(vars)

#Function use to define coxph model variables
#function for model with full list of variables
coxph_full_partial <-  purrr::partial(survival::coxph,
                                      formula = my_surv(CKD_AGEDIFF, INCIDENT_CKD) ~
                                      term + BL_AGE + MEN + BMI +
                                      PREVAL_DIAB.col_from_endpoints + SYSTM + BP_TREAT +
                                      PREVAL_HFAIL_STRICT.col_from_endpoints + CURR_SMOKE +
                                      PREVAL_AUTOIMMUN.col_from_endpoints,
                                      ties = "breslow")

#Function for cox model with minimum variable (age + sex)                      
coxph_minimum_partial <-  purrr::partial(survival::coxph,
                                         formula = my_surv(CKD_AGEDIFF, INCIDENT_CKD) ~
                                           term + BL_AGE + MEN,
                                         ties = "breslow")                                      

#Function formating for missing variables
cox_model_for_taxon <- function(df, term, fun, rank = NULL) {
  stopifnot(!missing(df), !missing(term), !missing(fun))
  message(term)
  { if (is.null(rank)) df else df[[rank]] } %>%
    dplyr::rename(term := {{term}}) %>% 
    fun(data = .)
}

##ALPHA DIVERSITY
#Univariate adj for sex and age
#Using alpha diversity to predict incidence of CKD, adjusting only for sex and age
tse %>%
  tse_meta(rownames = FALSE) %>% 
  cox_model_for_taxon(term = "shannon", fun = coxph_minimum_partial) %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  round_numeric_columns() %>%
  dplyr::mutate_if(is.numeric, ~round(.x, digits = 4)) %>% 
  DT::datatable(caption = "Alpha diversity")

tse %>%
  tse_meta(rownames = FALSE) %>%
  cox_model_for_taxon(term = "shannon", fun = coxph_full_partial) %>% 
  broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  round_numeric_columns() %>%
  dplyr::mutate_if(is.numeric, ~round(.x, digits = 4)) %>% 
  DT::datatable(caption = "Alpha diversity")

##03062924
#Approximate microbial data dimensions
df_counts <- tse %>%
  mia::subsetByPrevalentFeatures(detection = args$detection, prevalence = args$prevalence, as_relative = TRUE) %>% 
  SummarizedExperiment::assay("relabundance") %>%
  t %>%
  tibble::as_tibble()

tse

pca <- prcomp(df_counts)
n_independent_axes <- broom::tidy(pca, matrix = "eigenvalues") %>%
  dplyr::filter(cumulative < 0.9) %>%
  nrow %>%
  add(1)

#PCA vectors
df_pca <- pca$x %>%
  tibble::as_tibble() %>%
  dplyr::bind_cols(tse %>% tse_meta(rownames = FALSE))

tibble::tibble(model = glue::glue("PC{seq(n_independent_axes)}")) %>%
  dplyr::mutate(results = purrr::map(model, ~cox_model_for_taxon(df = df_pca, term = .x, fun = coxph_minimum_partial))) %>% 
  dplyr::mutate(results = purrr::map(results, ~broom::tidy(.x, exponentiate = TRUE, conf.int = TRUE))) %>%
  tidyr::unnest(results) %>% 
  dplyr::filter(stringr::str_detect(term, "term")) %>%
  dplyr::arrange(p.value) %>%
  dplyr::mutate(qval_my = my_adjust_p(p.value, n = n_independent_axes)) %>%
  dplyr::mutate(qval_fdr = p.adjust(p.value, method = "BH")) %>%
  dplyr::mutate_if(is.numeric, ~round(.x, digits = 4)) %>% 
  DT::datatable(caption = "PCA axes with GFR (minimum model)")

tibble::tibble(model = glue::glue("PC{seq(n_independent_axes)}")) %>%
  dplyr::mutate(results = purrr::map(model, ~cox_model_for_taxon(df = df_pca, term = .x, fun = coxph_full_partial))) %>% 
  dplyr::mutate(results = purrr::map(results, ~broom::tidy(.x, exponentiate = TRUE, conf.int = TRUE))) %>% 
  tidyr::unnest(results) %>% 
  dplyr::filter(stringr::str_detect(term, "term")) %>%
  dplyr::arrange(p.value) %>%
  dplyr::mutate(qval_my = my_adjust_p(p.value, n = n_independent_axes)) %>%
  dplyr::mutate(qval_fdr = p.adjust(p.value, method = "BH")) %>%
  dplyr::mutate_if(is.numeric, ~round(.x, digits = 4)) %>% 
  DT::datatable(caption = "PCA axes with GFR (full model)") 

#Common taxa
#Agglomerate prevalent features
subset_features_partial <-  purrr::partial(mia::subsetByPrevalentFeatures,
                                           x = tse,
                                           detection = args$detection,
                                           prevalence = args$prevalence,
                                           as_relative = TRUE) %>%
  purrr::possibly()

taxa_subsets <- c("Species") %>%
  rlang::set_names() %>%
  purrr::map(subset_features_partial)

taxa_subsets %>%
  purrr::map_df(~mia::getPrevalence(.x, detection = args$detection, as_relative = TRUE), .id = "rank") %>%
  tidyr::gather(Taxa, Abundance, -rank) %>%
  dplyr::arrange(Abundance) %>%
  dplyr::mutate_if(is.numeric, ~round(.x, digits = 4)) %>% 
  DT::datatable(caption = "Gut microbial abundances")

melt_tse <- function(x, method = "clr") {
  mia::transformAssay(x, method = method, pseudocount = 1) %>% 
    mia::meltAssay(add_row_data = TRUE, assay_name = method) %>%
    dplyr::mutate(FeatureID = glue::glue("GUT_{FeatureID}")) %>%
    dplyr::mutate(FeatureID = stringr::str_replace_all(FeatureID, c(" " = "_", "-" = "_"))) %>% 
    dplyr::select(SampleID, FeatureID, clr) %>%
    tidyr::spread(FeatureID, clr) %>% 
    dplyr::full_join(tse_meta(x), by = dplyr::join_by(SampleID == rownames))
}

dfs <- taxa_subsets %>%
  purrr::map(melt_tse, .progress = TRUE)

##COX's model
#Check available taxa
taxa_rank_list <- function(x, y) {
  list_names <- colnames(x) %>%
    stringr::str_subset("GUT_") %>%
    rlang::set_names()
}

taxa_rank_list <- dfs %>%
  purrr::imap(taxa_rank_list, .progress = TRUE)

#COX model adjusted for sex and age (minimum model)
taxa_rank_list %>%
  purrr::map_df(~tibble::tibble(taxa = .x), .id = "rank") %>% 
  tidyr::gather(rank, taxa) %>%
  dplyr::mutate(results = purrr::map2(rank, taxa, ~cox_model_for_taxon(df = dfs,
                                                                       term = .y,
                                                                       rank = .x,
                                                                       fun = coxph_minimum_partial), .progress = TRUE)) %>% 
  dplyr::mutate(results = purrr::map(results, ~broom::tidy(.x, exponentiate = TRUE, conf.int = TRUE)))  %>%
  tidyr::unnest(results) %>%
  dplyr::filter(stringr::str_detect(term, "term")) %>%
  dplyr::mutate(qval_my = my_adjust_p(p.value, n = n_independent_axes)) %>%
  dplyr::mutate(qval_fdr = p.adjust(p.value, method = "BH")) %>%
  dplyr::arrange(p.value) %>%
  dplyr::mutate_if(is.numeric, ~round(.x, digits = 4)) %>% 
  DT::datatable(caption = "Per taxa Cox's models (minimum)")

#COX model
#adjusted for all confounding 
taxa_rank_list %>%
  purrr::map_df(~tibble::tibble(taxa = .x), .id = "rank") %>% 
  tidyr::gather(rank, taxa) %>%
  dplyr::mutate(results = purrr::map2(rank, taxa, ~cox_model_for_taxon(df = dfs,
                                                                       term = .y,
                                                                       rank = .x,
                                                                       fun = coxph_full_partial), .progress = TRUE)) %>% 
  dplyr::mutate(results = purrr::map(results, ~broom::tidy(.x, exponentiate = TRUE, conf.int = TRUE)))  %>%
  tidyr::unnest(results) %>%
  dplyr::filter(stringr::str_detect(term, "term")) %>%
  dplyr::mutate(qval_my = my_adjust_p(p.value, n = n_independent_axes)) %>%
  dplyr::mutate(qval_fdr = p.adjust(p.value, method = "BH")) %>%
  dplyr::arrange(p.value) %>%
  dplyr::mutate_if(is.numeric, ~round(.x, digits = 4)) %>% 
  dplyr::filter(qval_fdr < 0.5)
  DT::datatable(caption = "Per taxa Cox's models (full)")


  
  colData(tse_gg2_MGS_FR02)  
  
