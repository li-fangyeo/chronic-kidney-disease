#CKD figures
#Prevalence figures
#Figure 1 . Alpha diversity and individual variables in prevalent CKD
##12072024
#forest plot
library(readr)
alpha_nephro <- read_csv("alpha-nephro.csv")
#alpha_nephro$'' <- cut(alpha_nephro$'P', breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", "")) 
alpha_nephro
alpha_nephro$Outcome <- ifelse(is.na(alpha_nephro$Outcome), "", alpha_nephro$Outcome)
alpha_nephro$` ` <- paste(rep(" ", 20), collapse = " ")

# Create a confidence interval column to display
alpha_nephro$`β (95% CI)` <- ifelse(is.na(alpha_nephro$se), "",
                           sprintf("%.2f (%.2f to %.2f)",
                                   alpha_nephro$Shannon, alpha_nephro$conf.low, alpha_nephro$conf.high))
alpha_nephro

library(grid)
library(forestploter)
library(gridExtra)
#make the words center
tm <- forest_theme(base_size = 16,
                   core=list(fg_params=list(hjust=0.5,
                                            x= 0.5)),
                   colhead=list(fg_params=list(hjust= 0.5, x = 0.5)),
                   #rowhead=list(fg_params=list(hjust=0, x=0)),
                   ci_Theight = 0.2,
                   ci_lwd = 2)

forest <- forest(alpha_nephro[,c(1,2, 9, 8, 4)],
            est = alpha_nephro$Shannon,
            lower = alpha_nephro$conf.low, 
            upper = alpha_nephro$conf.high,
            ci_column = 4,
            #ref_line = 1,
            xlab = " β ",
            xlim = c(-2, 2),
            theme = tm)

# Print plot
plot(forest)
ggplot2::ggsave(filename = "forest.pdf", 
                plot = forest,
                #dpi = 300,
                width = 10,
                height = 4,
                units = "in" )

