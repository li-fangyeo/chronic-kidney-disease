#Figure 2A Heatmap of associated taxa with variables
library(miaViz)
library(scater)
library(ggplot2)
library(pheatmap)
library(dplyr)

gfr<- df_lm_gfr_results %>%
  dplyr::arrange(p.value) %>% 
  dplyr::filter(qval_fdr < 0.05)
  #DT::datatable(caption = "Linear model for GFR")

uac <- df_lm_uac_results %>%
  dplyr::arrange(p.value) %>% 
  dplyr::filter(qval_fdr < 0.05)
  #DT::datatable(caption = "Linear model for UAC")

crea<- df_lm_crea_results %>%
  dplyr::arrange(p.value) %>% 
  dplyr::filter(qval_fdr < 0.05)
  #DT::datatable(caption = "Linear model for Serum Creatinine")

h<- uac%>% select(taxa, estimate) %>%
          mutate(newcol = "UAC")
h
i<- gfr %>% select(taxa, estimate) %>%
            mutate(newcol = "GFR")
i

j<- crea %>% select(taxa, estimate) %>%
              mutate(newcol = "CREA")
j
k<- rbind(h,i,j)
k
#d$stars <- cut(d$'qval_fdr', breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

k$taxa <- gsub("GUT_Species:", "",x=k$taxa)
k$taxa<- factor(k$taxa, levels=rev(sort(unique(k$taxa))))
k<- k %>% arrange(estimate) 

p <- ggplot(k, aes(x = newcol, y = taxa, fill = estimate)) +
  geom_tile(color = "black") +
  #geom_text(aes(label = estimate), color = "black", size = 6) +
  #geom_text(aes(label=stars), color="black", size=5) +
  coord_fixed() +
  scale_fill_gradient2(low = 'cornflowerblue', high = 'orange', mid = "white") +
  theme_classic()
taxa.heatmap <- p + theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x=element_text(size=6, angle = 90),
          axis.text.y=element_text(size=6))
taxa.heatmap 

ggplot(data = k, aes(x = taxa, y = estimate, fill = newcol)) +
  geom_bar(stat = "identity", position=position_dodge()) + 
  scale_fill_manual(values = c("darkslategrey", "orange")) +
  theme(axis.text.x = element_text(angle = 90))

library(patchwork)
alpha
taxa.heatmap
dots

nested <- ((alpha/dots)|taxa.heatmap)+
  plot_annotation(tag_levels = 'A')
nested

##18072024
##Circular heatmap
library(ggtreeExtra)
library(ggtree)
library(dplyr)
library(ggnewscale)
library(tibble)
library(ComplexHeatmap)
library(circlize)

a<- uac%>% select(taxa, estimate)
b<- gfr %>% select(taxa, estimate)
c<- crea %>% select(taxa, estimate)
e<- full_join(a,b, by = "taxa",suffix = c(".uac",".gfr"))
e<- full_join(e,c, by = "taxa")
colnames(e)[4] <- "Creatinine"
colnames(e)[2] <- "UACR"
colnames(e)[3] <- "eGFR"
as.data.frame(e)
e$taxa <- gsub("GUT_Species:", "",x=e$taxa) 
e <- e%>% replace(is.na(.), 0) %>% arrange(taxa) %>% column_to_rownames(var = "taxa")
e
#e<- e %>% relocate(UACR, .after=Creatinine)
e <- subset(e, select = -c(UACR))
e

write.csv(e, "taxa-sp-DAA.csv")
#circular
circos.clear()
col_fun1 = colorRamp2(c(-2, 0, 2), c("orange", "white", "darkslategrey"))
circos.par(gap.after = c(10))

circos.heatmap(e, col = col_fun1, rownames.side = "outside"
               , cluster = TRUE, rownames.cex = 0.65) #track.height = 0.2)

#add legend
lgd = Legend(title = "β", col_fun = col_fun1)
circle_size = unit(1, "snpc")
pushViewport(viewport(x = 0.2, y = 0.5, width = circle_size, height = circle_size,
                      just = c("left", "center")))
draw(lgd, x = circle_size, just = "left")



#adding the group names
#circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
 # if(CELL_META$sector.numeric.index == 1) { # the last sector
  #  cn = colnames(e)
   # n = length(cn)
    #circos.text(rep(CELL_META$cell.xlim[1], n) + convert_x(1, "mm"), 
     #           1:n - 0.5, cn, 
      #          cex = 0.5, adj = c(0, 0.5), facing = "clockwise")
#  }
#}, bg.border = NA)



  


