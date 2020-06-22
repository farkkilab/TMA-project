### heatmap ----
## creating stuff 

library(RColorBrewer)
library(pheatmap)
library(circlize)
library(ComplexHeatmap)
tsi <- read.csv("C:/LocalData/nuplyyt/cyto/celltypes/round_april/tsi_for_plots.csv")

tsi <- tsi[,c(3:4,6,10:16,18:32,47:48)] #sample, ids, markers excl CD57 and CD45, HR and GCT
mc <- metaclusters[,c(4:5,7,51)]

snowflake <- function(x,y) {
  paste(x, y, sep = "_")
}

tsi$unique <- snowflake(tsi$cores, tsi$Cellid)
mc$unique <- snowflake(mc$cores, mc$Cellid)
tsi$metacluster <- tsi$GlobalCellType
tsi$metacluster <- as.character(tsi$metacluster)
mc$Metacluster <- as.character(mc$Metacluster)

#this will take a while
for (i in 1:nrow(mc)) {
  tsi$metacluster[which(tsi$unique == mc$unique[i])] <- mc$Metacluster[i]
}

tsi <- tsi[tsi$metacluster != "Other",]
tsih <- matrix(1:23, nrow = 1)

for (i in unique(tsi$metacluster)) {
  celltype <- tsi[tsi$metacluster == i,]
  ct <- 1
  
  for (x in names(tsi)[4:25]) {
    marker <- celltype[[x]]
    ct <- cbind(ct, median(marker))
  }
  
  ct <- cbind(ct, nrow(celltype))
  tsih <- rbind(tsih, ct[2:24])
  
}

tsih <- tsih[2:11,]
colnames(tsih)[23] <- "Total"
colnames(tsih)[1:22] <- names(tsi)[4:25]
rownames(tsih) <- unique(tsi$metacluster)
rownames(tsih)[8] <- "Functional epithelial"
rownames(tsih)[10] <- "Other"
colnames(tsih)[2] <- "53BP1"


## actual heatmaps

coul <- colorRampPalette(brewer.pal(8, "Blues"))(25)

# heatmap(tsih[c(6,9,8,4,7,5,2),c(18,19,13,9:10,14,17,21)], scale="column",
#         Colv = NA, main="Heatmap scaled by column", col=coul)
# heatmap(tsih[c(10,1,3,9,8,6,7,5,4,2),c(22,18:20,1,3:8,12,15,16)], scale="column", 
#         Colv = NA, Rowv = NA, col=coul)
# heatmap(tsih[c(9,8,6,7,5,4,2),c(18,19,13,2,9:11,14,17,21)], scale="column", 
#         Colv = NA, Rowv = NA, col=coul)

Number <- as.data.frame(tsih[c(2,4,5,7,6,8,9),23])
names(Number) <- "Number"

#trying with pheatmap
# I think this is the one I used for the actual heatmap and the legend....
pheatmap(tsih[c(4,2,5,7,6,9,8),c(18,19,13,9:10,14,17,21)], cluster_rows=F,
         cluster_cols=F, color=coul, border_color=NA, scale="column",
         fontsize=15, cellwidth=30, cellheight=30) #annotation_row=Number


#trying with complexheatmap
#....and then I would take the bar graph from here and combine them in PowerPoint

#Rowann <- HeatmapAnnotation(df=ann, which="row") #col=colours
ha <- rowAnnotation(cells=anno_barplot(Number$Number[order(Number$Number, decreasing = T)], width = unit(2, "cm")))
#colnames(percentages_of_cells) <- c("tumor", "stromal", "immune", "other")
#col_fun = colorRamp2(c(0,5,30, 60), c("lavenderblush","mistyrose","indianred1" ,"red2"))
# hmap <- Heatmap(tsih[c(2,4,5,7,6,8,9),c(18,19,13,9:10,14,17,21)], name="asd",
#                 row_names_gp = gpar(fontsize = 8),
#                 row_dend_width = unit(1.3, "cm"),row_title="patients",
#                 show_row_dend=T, show_column_dend=T,
#                 column_title = "Percentages of immune cells",
#                 height = unit(16, "cm") , right_annotation = ha, width = unit(4, "cm"),
#                 border="white", rect_gp = gpar(col = "white", lwd = 4), col=col_fun,
#               cell_fun = function(j, i, x, y, width, height, fill) {
# if(percentages_of_cells[i, j] > 0)
#   grid.text(sprintf("%.1f", percentages_of_cells[i, j]), x, y, gp = gpar(fontsize = 5))
# })
# draw(hmap)

coul <- colorRampPalette(brewer.pal(8, "Blues"))(25)
col_scale <- colorRamp2(c(min(tsih[c(4,2,5,7,6,9,8),c(18,19,13,9:10,14,17,21)]),
                          mean(tsih[c(4,2,5,7,6,9,8),c(18,19,13,9:10,14,17,21)]),
                          max(tsih[c(4,2,5,7,6,9,8),c(18,19,13,9:10,14,17,21)])),
                        c("lavenderblush","mistyrose","indianred1"))
col_scale = colorRamp2(c(0,5,30, 60), c("lavenderblush","mistyrose","indianred1" ,"red2"))
Heatmap(tsih[c(4,2,5,7,6,9,8),c(18,19,13,9:10,14,17,21)], right_annotation = ha,
        col=col_scale, cluster_rows=F, cluster_columns=F)
# so I don't think I used this heatmap above at all, just the annotation


### CD4 high and low -----
CD4hilo <- cbind(prop.i$Sample[order(prop.i$CD4)],
                 prop.i$CD4[order(prop.i$CD4)],
                 prop.i$Type[order(prop.i$CD4)],
                 c(rep(1, 44)))
CD4hilo <- as.data.frame(CD4hilo)

names(CD4hilo) <- c("Sample", "CD4", "HR", "cores")

for (s in CD4hilo$Sample){
  cores <- unique(tmasample$cores[tmasample$Sample==s])
  n <- length(cores)
  corelist <- ""
  for (l in 1:n){
    corelist <- paste(corelist, cores[l], sep = ", ")
  }
  
  CD4hilo$cores[CD4hilo$Sample==s] <- substr(corelist, 3, 20)
}

write.csv(CD4hilo, "C:/LocalData/nuplyyt/cyto/celltypes/images/maybe_paper/CD4highlow.csv",
          row.names = F)

### barplot ----

library(RColorBrewer)
library(scales)
library(cowplot)

tumor <- read.csv("C:/LocalData/nuplyyt/cyto/celltypes/round_april/tsi_for_plots_tumor_new.csv")
immune <- read.csv("C:/LocalData/nuplyyt/cyto/celltypes/round_april/tsi_for_plots_immune_new.csv")
stroma <- read.csv("C:/LocalData/nuplyyt/cyto/celltypes/round_april/stroma_for_plots_nonendothelial.csv")
endothelial <- read.csv("C:/LocalData/nuplyyt/cyto/celltypes/round_april/stroma_for_plots_endothelial.csv")

tumor$Subtype <- as.character(tumor$Subtype)
tumor$Subtype <- c(rep("Tumor",nrow(tumor)))
immune$Subtype <- as.character(immune$Subtype)
immune$Subtype <- c(rep("Immune",nrow(immune)))

names(stroma)[50] <- "Subtype"
stroma$Subtype <- as.character(stroma$Subtype)
stroma$Subtype <- c(rep("Stroma",nrow(stroma)))
names(endothelial)[50] <- "Subtype"
endothelial$Subtype <- as.character(endothelial$Subtype)
endothelial$Subtype <- c(rep("Endothelial", nrow(endothelial)))

barp <- rbind(tumor[,c(3,48)], immune[,c(3,49)], stroma[,c(4,50)], endothelial[,c(4,50)])

forbarplot <- data.frame("Sample"=c(rep(unique(barp$Sample),4)),
                         "Subtype"=c(rep("Tumor",44),rep("Immune",44),rep("Stroma",44),rep("Endothelial",44)),
                         "Number"=rep(0,176),
                         "Tumorcont"=rep(0,176))
for (i in 1:176){
  sample <- forbarplot[i,1]
  patient <- barp[barp$Sample==sample,]
  
  for (celltype in c("Tumor","Immune","Stroma","Endothelial")){
    nro <- nrow(patient[patient$Subtype == celltype,])
    forbarplot$Number[forbarplot$Sample==sample & forbarplot$Subtype==celltype] <- nro
  }
  
  forbarplot[i,4] <- nrow(patient[patient$Subtype == "Tumor",])/nrow(patient)
}

forbarplot$Sample <- as.factor(forbarplot$Sample)
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(4)

#lvls <- names(sort(tapply(forbarplot$Subtype == "Tumor", forbarplot$Sample, mean)))
#lvls <- names(sort(table(forbarplot[forbarplot$Subtype == "Tumor", "Sample"])))
#factor(x, levels = lvls)

#forbarplot2 <- forbarplot[order(forbarplot$Tumorcont, decreasing = TRUE),]
#forbarplot2$TCF <- factor(forbarplot2$Tumorcont, levels= unique(forbarplot2$Tumorcont))

forbarplot2 <- forbarplot
forbarplot2$Sample = with(forbarplot2[forbarplot2$Subtype == "Tumor", ], reorder(Sample, -Tumorcont, mean))

#actual barplot
h1 <- ggplot(data=forbarplot2, aes(x=Sample, fill=Subtype, y=Number))+
  geom_bar(stat="identity", position="fill")+
  theme_minimal()+
  scale_fill_manual(values=mycolors)+
  #theme(axis.text.x = element_text(angle=45))+
  #scale_x_discrete(name="Patient ID")+
  scale_y_continuous(name="Celltype proportion",labels=percent)+
  guides(fill=guide_legend(title="Celltype"))+
  theme(axis.text.x = element_blank(),
    axis.ticks = element_blank())+
  scale_x_discrete(breaks=NULL, name=NULL)

#BRCA status
clinical$Sample <- as.factor(clinical$Sample)
forbarplot2 <- merge(forbarplot2, clinical, by="Sample")

h2 <- ggplot(forbarplot2)+
  geom_bar(mapping = aes(x = Sample, y = 1, fill = Type), 
           stat = "identity", 
           width = 1)+
  scale_fill_manual(values=c("blue","red"))+
  theme(panel.spacing.x = unit(1, "mm"))+
  guides(fill=guide_legend(title="HR status"))+
  theme_void()

#PFI
h3 <- ggplot(forbarplot2)+
  geom_bar(mapping = aes(x = Sample, y = 1, fill = PFI), 
           stat = "identity", 
           width = 1)+
  scale_fill_manual(values=c("#3399CC","#003366"))+
  theme(panel.spacing.x = unit(1, "mm"))+
  guides(fill=guide_legend(title="PFI"))+
  theme_void()

#h3 <- h3 + theme(legend.position = "none")
#long "#3399CC" and "short "#003366"

plot_grid(h1, h2, h3, align = "v", ncol = 1, axis = "tb", rel_heights = c(15, 0.5, 0.5))



### misc ----
is.brca <- function(core){
  patient <- unique(tmasample$Sample[tmasample$cores==core])
  brca <- clinical$Type[clinical$Sample==patient]
  print(paste0("Patient ID ", patient, ", BRCA-status ", brca))
}
