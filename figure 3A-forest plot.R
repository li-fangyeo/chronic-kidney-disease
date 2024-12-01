#Alpha forest plot 
library(ggplot2)
library(survminer)
library(survival)
library(tibble)
library(mia)
library(survivalAnalysis)
library(readr)
##manually change the column names and reordered based on ascending order and change p value to <0.001
cox <- read_csv("cox-CKD.csv")

#cox$'' <- cut(cox$'P', breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 
cox

cox$'' <- paste(rep(" ", 30), collapse = " ")

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

cox.CKD <- forest(cox[,c(2,10,9,6)],
                est = cox$estimate,
                lower = cox$conf.low, 
                upper = cox$conf.high,
                ci_column = 3,
                ref_line = 1,
                xlim = c(0,8),
                ticks_at = c(0.25,0.5,1,2,4,6,8),
                theme = tm ,
                x_trans = c("log10"),
                xlab = "Hazard ratio")

# Print plot
plot(cox.CKD)
ggplot2::ggsave(filename = "coxCKD28112024.pdf", 
                plot = cox.CKD,
                #dpi = 300,
                width = 10,
                height = 4,
                units = "in" )

                