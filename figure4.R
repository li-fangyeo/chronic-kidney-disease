#Figure 3A forest plot
library(grid)
library(forestploter)
library(dplyr)

t<- taxa_rank_list %>%
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
  dplyr::mutate_if(is.numeric, ~round(.x, digits = 4)) #%>% 
  DT::datatable(caption = "Per taxa Cox's models (minimum)")

f <- forest(t,
              est = t$estimate,
              lower = t$conf.low, 
              upper = t$conf.high,
              sizes = t$std.error,
              ci_column = 4,
              ref_line = 1,
             #arrow_lab = c("Placebo Better", "Treatment Better"),
              xlim = c(0, 4),
              ticks_at = c(0.5, 1, 2, 3),
              footnote = "This is the demo data. Please feel free to change\nanything you want.")
  
# Print plot
plot(f)
