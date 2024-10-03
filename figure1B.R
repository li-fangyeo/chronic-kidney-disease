#Figure 2 - beta diversity plot
#Figure 3b 
library(ggfortify) 
library(mia)
#convert tse to a dataframe.
pca
pc_df <- as.data.frame(colData(tse))

pca.plot <- autoplot(pca, 
                  data = pc_df, 
                  colour = 'GFR') 

dots<- pca.plot + 
      theme(title = element_text(size = 15)) +
  scale_color_gradient(high = "darkslategrey", low = "orange")

dots <- dots + labs(color = "eGFR") + theme_classic()+ 
        theme(legend.title = element_text(size = 16),
              legend.text = element_text(size = 16),
              axis.text=element_text(size=14),
              axis.title = element_text(size = 16))
dots
#17072024
##For two figures
library(ggplot2)
library("grid")
library("ggplotify")
library(patchwork)
forest
dots
forest<- as.ggplot(forest)

nested <- (forest / dots) +
  plot_annotation(tag_levels = 'A')
nested
ggplot2::ggsave(filename = "nested.pdf", 
                plot = nested,
                #dpi = 300,
                width = 10,
                height = 10,
                units = "in" )

####
library(patchwork)
alpha
taxa.heatmap
dots
#for three figures
nested <- ((alpha/dots)|taxa.heatmap)+
  plot_annotation(tag_levels = 'A')
nested


## PCA for Incident CKD - figure 3
pca.plot2 <- autoplot(pca, 
                      data = pc_df, 
                      colour = 'INCIDENT_CKD') 
dots2 <- pca.plot2 + 
  theme(title = element_text(size = 16)) +
  theme_classic() +
  scale_color_manual(values = c("darkslategrey", "orange"), 
                     labels = c("No","Yes")) +
  labs(color = "Incident CKD") 
  
  
dots2 + theme(legend.position = c(0.2, 0.8))+ 
        theme(legend.title = element_text(size = 16),
              legend.text = element_text(size = 16),
              axis.text=element_text(size=14),
              axis.title = element_text(size = 16))

##Table one
tse %>%
  tse_meta(rownames = FALSE) %>%
  dplyr::mutate(INCIDENT_CKD = factor(ifelse(INCIDENT_CKD == 1, "Incident CKD" , "No CKD"))) %>%
  mytableone(vars,fo =  ~ .| INCIDENT_CKD ) %>%
  write.csv("CKD-tableone.csv")
##calculate p-value from tableone
library(tableone)
a<- as.data.frame(colData(tse))
CreateTableOne(data = a, strata = "INCIDENT_CKD")
