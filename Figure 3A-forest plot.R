#Alpha forest plot 
library(ggplot2)
library(survminer)
library(survival)
library(tibble)
library(mia)
library(survivalAnalysis)
#Change into dataframe and format Incident into integer, otherwise error
a <- as.data.frame(colData(tse))
a$INCIDENT_CKD <- as.integer(a$INCIDENT_CKD)
#Change the row names
a <- a %>% dplyr::mutate_at(vars(BMI, BL_AGE, SYSTM), ~./10) %>%
  dplyr::rename("Age" = BL_AGE ,
              "Sex"= MEN ,
              BMI = "BMI",
              "Diabetes" = PREVAL_DIAB.col_from_endpoints,
              "Systolic blood pressure" = SYSTM,
              "Antihypertensive medication" = BP_TREAT ,
              "Heart failure" = PREVAL_HFAIL_STRICT.col_from_endpoints,
              "Prevalent CKD" = PREVAL_CKD,
              "Incident CKD" = INCIDENT_CKD,
              "CKD time (yrs)" = CKD_AGEDIFF ,
              "Smoker" = CURR_SMOKE,
              "Eastern Finland" = EAST,
              "Shannon diversity" = shannon,
              "Autoimmune Disease" = PREVAL_AUTOIMMUN.col_from_endpoints)
             

b<- coxph(Surv(`CKD time (yrs)`, `Incident CKD`) ~
            `Shannon diversity` + Age + Sex + BMI +
            Diabetes + `Systolic blood pressure` + `Antihypertensive medication` +
            `Heart failure` + Smoker +
            `Autoimmune Disease`, data = a)
ggforest(b) + theme_classic()
library(readr)
##manually change the column names and reordered based on ascending order
cox <- read_csv("cox-CKD.csv")

#cox$'' <- cut(cox$'P', breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 
cox

cox$'' <- paste(rep(" ", 28), collapse = " ")

# Create a confidence interval column to display
cox$`HR (95% CI)` <- ifelse(is.na(cox$std.error), "",
                                     sprintf("%.2f (%.2f to %.2f)",
                                             cox$estimate, cox$conf.low, cox$conf.high))
cox

library(grid)
library(forestploter)
library(gridExtra)
#make the words center
tm <- forest_theme(base_size = 15,
                   rowhead=list(fg_params=list(hjust=0, x=0)),
                   #core=list(fg_params=list(hjust = 0.5, x = 0.5)),
                   colhead=list(fg_params=list(hjust=0.5, x=0.5)),
                   ci_Theight = 0.2,
                   ci_lwd = 2)

cox.CKD <- forest(cox[,c(1,9,8,5)],
                est = cox$estimate,
                lower = cox$conf.low, 
                upper = cox$conf.high,
                ci_column = 3,
                ref_line = 1,
                xlim = c(0,10),
                ticks_at = c(0.25,0.5,1,2,4,6,8),
                theme = tm ,
                x_trans = c("log10"),
                xlab = "Hazard ratio")

# Print plot
plot(cox.CKD)
ggplot2::ggsave(filename = "coxCKD19092024.pdf", 
                plot = cox.CKD,
                #dpi = 300,
                width = 10,
                height = 4,
                units = "in" )

                