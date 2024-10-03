#Beta diversity
#Error
outcome <- tse_subset_by_sample$PREVAL_CKD
covs_names <- names(vars)
covs_formula <- paste(covs_names, collapse = " + ")

#13062024
#make smaller dataset
library(mia)
library(dplyr)
library(magrittr)
library(stringr)
tse
#make a small dataset
tse_subset_by_sample <- tse[ , tse$BMI > 35]
tse_subset_by_sample

tse_subset_by_sample %>% mergeFeaturesByRank("Species") %>%
  transformSamples (method = "relabundance")
assay <- t(assay(tse_subset_by_sample,"relabundance"))

#Get abundance table
abund <- assay(tse_subset_by_sample, "relabundance")
abund_df <- as.data.frame(t(abund))

#BC dissimilarity
bc_matrix <- vegan::vegdist(abund_df, method = "bray")

set.seed(42)
vegan::adonis2((as.formula(str_glue("bc_matrix ~ INCIDENT_HTN"))), by = 'margin', 
               data = colData(tse_subset_by_sample), permutations =999)
vegan::adonis2((as.formula(str_glue("bc_matrix ~ INCIDENT_HTN  + BL_AGE + MEN + BMI + PREVAL_DIAB.col_from_pheno + HFC +
                    PREVAL_HFAIL_STRICT.col_from_pheno + CURR_SMOKE"))), by = 'margin', 
               data = colData(tse_subset_by_sample), permutations =999)

#Full dataset
tse_species <- tse %>% mergeFeaturesByRank("Species") %>%
  transformSamples (method = "relabundance")

#Get abundance table
abund <- assay(tse_species, "relabundance")
abund_df <- as.data.frame(t(abund))

#BC dissimilarity
bc_matrix <- vegan::vegdist(abund_df, method = "bray")

#order of the variable of interest
set.seed(423)
permanova_partial <- vegan::adonis2((as.formula(str_glue("bc_matrix ~ BL_AGE + MEN + INCIDENT_HTN"))), by = 'margin', 
               data = colData(tse_species), permutations =999)
permanova_full <- vegan::adonis2((as.formula(str_glue("bc_matrix ~ INCIDENT_HTN  + BL_AGE + MEN + BMI + PREVAL_DIAB.col_from_pheno + HFC +
                    PREVAL_HFAIL_STRICT.col_from_pheno + CURR_SMOKE"))), by = 'margin', 
               data = colData(tse_species), permutations =999)


#Learning codes
tse %>%
  tse_meta(rownames = FALSE) %>%
  dplyr::mutate(INCIDENT_HTN = factor(ifelse(INCIDENT_HTN == 1, "Incident HTN" , "Without HTN"))) %>%
  mytableone(vars,fo =  ~ .| INCIDENT_HTN )

#dplyr::case_when(MEN == 1 ~ ”Women”, MEN == 1 ~ ”Men”, )
#dplyr::mutate(MEN = factor(case_when(MEN == 1 ~ ”xxx”, MEN == 0 ~ ”yyy”, TRUE ~ ”zzz”))

tse %>%
  tse_meta(rownames = FALSE) %>%
  dplyr::mutate(INCIDENT_HTN = factor(ifelse(INCIDENT_HTN == 1, "Incident HTN" , "Without HTN"))) %>%
  mytableone(vars,fo =  ~ .| INCIDENT_HTN )

#counting UAC
pc_df %>% dplyr::filter(BL_AGE > 80) %>% dplyr::mutate(no_UAC = length(BL_AGE[BL_AGE >80]))
k %>% dplyr::filter(!is.na(UAC)) %>% dplyr::mutate(no_UAC = length(UAC[UAC >3]))
row.names(pc_df)