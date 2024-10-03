library(survival)
a <- as.data.frame(colData(tse))
a$INCIDENT_CKD <- as.numeric(a$INCIDENT_CKD)
cox_model <- coxph(Surv(BL_AGE, BL_AGE + CKD_AGEDIFF, INCIDENT_CKD) ~ shannon + MEN +
                     SYSTM + BP_TREAT + BMI + PREVAL_DIAB.col_from_endpoints +
                     PREVAL_HFAIL_STRICT.col_from_endpoints + CURR_SMOKE + 
                     PREVAL_AUTOIMMUN.col_from_endpoints, data = a)
summary(cox_model)


#visualizing outliers
boxplot(a$shannon)
a$shannon_log <- log10(a$shannon)
a$shannon_cube <-  a$shannon^(1/3)
par(mfrow=c(3,1))
i<- hist(a$shannon)
j<- hist(a$shannon_log)
k<- hist(a$shannon_cube)

#adding functional data to tse
library(mia)
joined_pathabundance <- read_delim("joined_pathabundance.tsv", 
  delim = "\t", escape_double = FALSE, trim_ws = TRUE)
#rename 1st column name
joined_pathabundance
a <- as.data.frame(colData(tse))
# Get the column names (excluding the first one, which is probably the path/feature column)
sample_names <- colnames(joined_pathabundance)[-1]

# Clean sample names by removing after .R1
cleaned_sample_names <- gsub("\\.R1.*","", sample_names, perl = TRUE)

# Update the column names with the cleaned sample names
colnames(joined_pathabundance)[-1] <- cleaned_sample_names

# Find the intersection of sample names between the cleaned joined_path_abundance sample names and FINRISK sample data
common_samples <- intersect(a$Barcode, cleaned_sample_names)

# Subset the joined_pathabundance file to include only the common samples
joined_path_abundance_subset <- joined_pathabundance %>%
  select(all_of(c("Pathway", common_samples)))  # "path_feature_column" should be replaced with actual first column name

# Write the subsetted data to a new file (if needed)
write_tsv(joined_path_abundance_subset, "joined_pathabundance_cleaned.tsv")

# Import the subsetted HUMAnN data (optional depending on your pipeline)
tse.hu <- mia::importHUMAnN("joined_pathabundance_cleaned.tsv")

# Add colData (e.g., sample metadata from FINRISK)
colData(tse.hu) <-DataFrame(a)

b <- as.data.frame(colData(tse.hu))
saveRDS(tse.hu, "tse_h3_functional.rds")

#27092024
#grep "s_" and filter to keep rows that are NOT zero
j <- joined_path_abundance_subset %>% filter(grepl('.s__', Pathway)) %>%
  filter(rowSums(across(where(is.numeric)))!=0)

#arrange by descending order 
mx <- j %>% arrange(desc(`820003769-1`)) %>%
  slice_head(n = 3) %>%
  t() %>%
  as.data.frame()
mi<- j %>% arrange(desc(`820003769-1`)) %>%
  slice_tail(n = 3)  %>%
  t()%>%
  as.data.frame()

colnames(mx) <- mx[1,]
mx <- mx[-1, ] 
# Column indices to be converted
i <- c(1 : 3)  
mx[, i] <- apply(mx[, i], 2, function(x) as.numeric(as.character(x)))
sapply(mx, class)

colnames(mi) <- mi[1,]
mi <- mi[-1, ] 
i <- c(1 : 3)  
mi[, i] <- apply(mi[, i], 2, function(x) as.numeric(as.character(x)))
sapply(mi, class)

hist(mx$`UNINTEGRATED|g__Bacteroides.s__Bacteroides_uniformis`)  

#Log (x+1) to all columns then plot histogram
mx_log <- mx %>%  mutate(across(c(1:3), function(x) log(x+1)))
hist(mx_log$`UNINTEGRATED|g__Bacteroides.s__Bacteroides_uniformis`)

