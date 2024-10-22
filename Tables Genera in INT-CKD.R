#Table for top 10 genus in Incident CKD
in.ckd<- taxa_rank_list %>%
  purrr::map_df(~tibble::tibble(taxa = .x), .id = "rank") %>% 
  tidyr::gather(rank, taxa) %>%
  dplyr::mutate(results = purrr::map2(rank, taxa, ~cox_model_for_taxon(df = dfs,
                                                                       term = .y,
                                                                       rank = .x,
                                                                       fun = coxph_full_partial), .progress = TRUE)) %>% 
  dplyr::mutate(results = purrr::map(results, ~broom::tidy(.x, exponentiate = TRUE, conf.int = TRUE)))  %>%
  tidyr::unnest(results) %>%
  dplyr::filter(stringr::str_detect(term, "term")) %>%
  #dplyr::mutate(qval_my = my_adjust_p(p.value, n = n_independent_axes)) %>%
  dplyr::mutate(qval_fdr = p.adjust(p.value, method = "BH")) %>%
  dplyr::arrange(p.value) %>%
  dplyr::mutate_if(is.numeric, ~round(.x, digits = 3)) #%>% 
  DT::datatable(caption = "Per taxa Cox's models (full)")
  
head(in.ckd)
in.ckd$taxa <- gsub("GUT_", "",x=in.ckd$taxa)
in.ckd <- as.data.frame(in.ckd)
k <- in.ckd[1:10,]
k <-k[, c("taxa", "estimate","conf.low","conf.high","p.value", "qval_fdr")]
k
k$"95% CI" <- paste(k$conf.low, "-", k$conf.high)
k
k<-k[,-(3:4)]

#Rearrange and rename column names
library(purrr)
library(dplyr)
k<- k[, c(1, 2, 5, 3,4)]
k<- k %>% 
  rename(
   Species = taxa,
   HR = estimate,
   "p-value" = p.value,
   "FDR q-value" = qval_fdr )
k
write.table(k, file = "Top10IncidentCKDsp.csv",
            sep = "\t", row.names = F)

#Species
in.ckd.sp<- taxa_rank_list %>%
  purrr::map_df(~tibble::tibble(taxa = .x), .id = "rank") %>% 
  tidyr::gather(rank, taxa) %>%
  dplyr::mutate(results = purrr::map2(rank, taxa, ~cox_model_for_taxon(df = dfs,
                                                                       term = .y,
                                                                       rank = .x,
                                                                       fun = coxph_full_partial), .progress = TRUE)) %>% 
  dplyr::mutate(results = purrr::map(results, ~broom::tidy(.x, exponentiate = TRUE, conf.int = TRUE)))  %>%
  tidyr::unnest(results) %>%
  dplyr::filter(stringr::str_detect(term, "term")) %>%
  #dplyr::mutate(qval_my = my_adjust_p(p.value, n = n_independent_axes)) %>%
  dplyr::mutate(qval_fdr = p.adjust(p.value, method = "BH")) %>%
  dplyr::arrange(p.value) %>%
  dplyr::mutate_if(is.numeric, ~round(.x, digits = 3)) #%>% 
DT::datatable(caption = "Per taxa Cox's models (full)")

head(in.ckd.sp)
in.ckd.sp$taxa <- gsub("GUT_", "",x=in.ckd.sp$taxa)
in.ckd.sp <- as.data.frame(in.ckd.sp)
k <- in.ckd.sp[1:10,]
k <-k[, c("taxa", "estimate","conf.low","conf.high","p.value", "qval_fdr")]
k
k$"95% CI" <- paste(k$conf.low, "-", k$conf.high)
k
k<-k[,-(3:4)]

#Rearrange and rename column names
library(purrr)
library(dplyr)
k<- k[, c(1, 2, 5, 3,4)]
k<- k %>% 
  rename(
    Genus = taxa,
    HR = estimate,
    "p-value" = p.value,
    "FDR q-value" = qval_fdr )
k
write.table(k, file = "Top10IncidentCKDspecies.csv",
            sep = "\t", row.names = F)
