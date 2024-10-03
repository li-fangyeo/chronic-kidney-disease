#Alpha forest plot 
library(ggplot2)
library(survminer)
library(survival)
library(tibble)
library(mia)
library(survivalAnalysis)
#Change into dataframe and format Incident into integer, otherwise error
a <- as.data.frame(colData(tse))
a$INCIDENT_HTN <- as.integer(a$INCIDENT_HTN)

#Change the row names
a <- a %>% dplyr::mutate_at(vars(BMI, BL_AGE, SYSTM), ~./10) %>%
  dplyr::rename("Age" = BL_AGE ,
                "Sex"= MEN ,
                BMI = "BMI",
                "Diabetes" = PREVAL_DIAB.col_from_pheno,
                "Systolic blood pressure" = SYSTM,
                "Diastolic blood pressure" = DIASM,
                "Healthy Food Choices" = HFC,
                "Antihypertensive medication" = BP_TREAT ,
                "Heart failure" = PREVAL_HFAIL_STRICT.col_from_pheno,
                "Prevalent HTN" = PREVAL_HTN,
                "Incident HTN" = INCIDENT_HTN,
                "Hypertension time (yrs)" = HTN_AGEDIFF ,
                "Smoking" = CURR_SMOKE,
                "Eastern Finland" = EAST,
                "Shannon diversity" = shannon)

b<- coxph(Surv(`Hypertension time (yrs)`, `Incident HTN`) ~
            `Shannon diversity` + Age + Sex + BMI +
            Diabetes + `Healthy Food Choices` + `Heart failure` + Smoking,
            data = a)
ggforest(b) + theme_classic()
summary(b)

#PERMANOVA table
rda_info %>%
knitr::kable() %>%
write.csv("permanova-full.csv")

##calculate p-value from tableone
library(tableone)
CreateTableOne(data = a, strata = "INCIDENT_CKD")


##12072024
#forest plot
library(readr)
alpha <- read_csv("coxHTN.csv")
#alpha <- cox_incidentCKD
#alpha$'' <- cut(alpha$'P', breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 
alpha
#alpha$Test <- ifelse(is.na(alpha_nephro$Test), "", alpha_nephro$Test)
alpha$` ` <- paste(rep(" ", 20), collapse = " ")

# Create a confidence interval column to display
alpha$`HR (95% CI)` <- ifelse(is.na(alpha$se), "",
                                     sprintf("%.2f (%.2f to %.2f)",
                                             alpha$`exp(coef)`, alpha$`lower .95`, alpha$`upper .95`))
alpha

library(grid)
library(forestploter)
library(gridExtra)
#make the words center
tm <- forest_theme(base_size = 15,
                   rowhead=list(fg_params=list(hjust=0, x=0)),
                   colhead=list(fg_params=list(hjust=0.5, x=0.5)),
                  #core=list(fg_params=list(hjust=0.5, x=0.5)),
                   ci_Theight = 0.3)

alpha_ft <- forest(alpha[,c(1, 8, 7, 6)],
                est = alpha$`exp(coef)`,
                lower = alpha$`lower .95`, 
                upper = alpha$`upper .95`,
                #sizes = alpha$se,
                ci_column = 3,
                ref_line = 1,
                xlim = c(0, 5),
                x_trans = c("log10"),
                theme = tm)

# Print plot
plot(alpha_ft)

myggsave(alpha)

