# I really hope no one's gonna use this because it is a total mess
# I also didn't think it would be used so I didn't clean it up very well; message me if you need help

### functions ----

#install.packages("RColorBrewer")
#install.packages("ggplot2")
#install.packages("caret")
library(RColorBrewer)
library(ggplot2)
library(caret)

#I always changed the function when doing different kinds of umaps so it might be a bit messy
plot.umap <- function(df, colorvar, save=FALSE, fname = "umap_plot.pdf", xlim = 100, ylim=NULL, ratio = 1) {
  ntypes <- length(unique(df[,colorvar]))
  mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(ntypes)
  p <- ggplot(data.frame(df), aes_string(x="X1", y="X2", color=colorvar)) +
    geom_point(size=0.5, stroke=0, alpha=0.7) +
    scale_color_manual(values=mycolors) +
    theme_bw() +
    theme(legend.position="right", legend.text = element_text(size=5), legend.title = element_blank()) + #legend.title =element_text(size=7)
    #theme(legend.position="none") +
    coord_fixed(ratio=ratio) +
    guides(colour = guide_legend(override.aes = list(size=5, shape=15, alpha=1), direction="horizontal", ncol=2, label.position="bottom", byrow=F)) + #nrow=1
    xlab("umap1") + ylab("umap2")
  print(p)
  if(save) ggsave(fname)
}

plot.by.markers <- function(df, channels, save=FALSE, fname = paste("umap_julia_", channels, ".pdf", sep=""), xlim = NULL, ylim=NULL, ratio = 1) {
  if (save) pdf(fname)
  for( channel in channels) {
    p <- ggplot(data.frame(df), aes_string(x="X1", y="X2", color=channel)) +
      geom_point(size=1, stroke=0, alpha=0.3) +
      scale_color_gradientn(colours=rainbow(5)) +
      coord_fixed(ratio=ratio) +
      xlab("umap1") + ylab("umap2") +
      theme_bw()
    print(p)
  }
  if (save) dev.off()
}


plot.by.log <- function(df, channels, save=FALSE, fname = "umap_plot.pdf", xlim = NULL, ylim=NULL, ratio = 0.3) {
  if (save) pdf(fname)
  for( channel in channels) {
    p <- ggplot(data.frame(df), aes_string(x="X1", y="X2", color=channel)) +
      geom_point(size=0.6, stroke=0, alpha=0.1) +
      scale_color_gradientn(colours=rainbow(5), trans="pseudo_log", na.value = "#FF000011") +
      coord_fixed(ratio = ratio)  +
      xlab("umap1") + ylab("umap2") +
      theme_bw()
    print(p)
  }
  if (save) dev.off()
}


### creating for.umap objects ------
# tsi
all_mark <- read.delim("C:/LocalData/nuplyyt/cyto/celltypes/round_april/tsi_Cyto-new_clustering_0204-corrected_set_of_markers-results.csv")
tmasample <- read.csv("C:/LocalData/nuplyyt/cyto/celltypes/tma1_sample.csv")

all_mark$GlobalCellType <- all_mark$label

#old but don't remove (see below)
all_mark$GlobalCellType[which(all_mark$GlobalCellType %in% c(3,4,6,9,10,12,13,15,19))] <- "Immune"
all_mark$GlobalCellType[which(all_mark$GlobalCellType %in% c(1,2,5,8,11,14,17,18))] <- "Tumor"

#new and fine on their own
all_mark$GlobalCellType[which(all_mark$GlobalCellType %in% c(7,16))] <- "Stroma"
all_mark$GlobalCellType[which(all_mark$label == 16)] <- "Endothelial"
#all_mark$global2 <- all_mark$GlobalCellType
#all_mark$global2[which(all_mark$label == 16)] <- "Endothelial" #if needed
all_mark$GlobalCellType[which(all_mark$GlobalCellType %in% c(20))] <- "Other"

#note! the ones below require first running the old globalcelltype, the unique stuff and the subtype stuff!
#note! run both of these to remove monocytes and other.immune from immnue cells!
all_mark$GlobalCellType[which(all_mark$GlobalCellType %in% c(3,4,6,9,10,12,13,15,19))] <- "Immune"
all_mark$GlobalCellType[which(all_mark$subtype %in% c("Tumor", "Monocytes", "Other.immune"))] <- "Tumor"

for.umap2 <- all_mark[all_mark$GlobalCellType != "Other",]


## unique stuff pt1 - global
# so here I create unique cell ID's for each cell so I can remove the ones I want
all_mark$core <- tmasample$cores
all_mark$cellid <- tmasample$Cellid

snowflake <- function(x,y) {
  paste(x, y, sep = "_")
}

all_mark$unique <- snowflake(all_mark$core, all_mark$cellid)
all_mark$subtype <- all_mark$GlobalCellType


# immune (old)
# labels.2 <- read.delim("C:/LocalData/nuplyyt/cyto/celltypes/round2/cyto_round2.csv")
# all_mark4i <- read.csv("C:/LocalData/nuplyyt/cyto/celltypes/round1/immune1_1.csv")
# all_mark.2 <- cbind(all_mark4i[,1:25],labels.2[,13])

# immune (new)
all_mark.2 <- read.delim("C:/LocalData/nuplyyt/cyto/celltypes/round_april/immune_Cyto-april_clustering_immune-no_tumor-1304-results.csv")
actual_markers.2 <- read.csv("C:/LocalData/nuplyyt/cyto/celltypes/round_april/tsi_for_plots_immune.csv")

all_mark.2$subtype <- all_mark.2$label

all_mark.2$subtype[which(all_mark.2$subtype %in% c(1))] <- "Monocytes"
all_mark.2$subtype[which(all_mark.2$subtype %in% c(8))] <- "M1 macrophages"
all_mark.2$subtype[which(all_mark.2$subtype %in% c(3))] <- "M2 macrophages"
all_mark.2$subtype[which(all_mark.2$subtype %in% c(5))] <- "CD8 T-cells"
all_mark.2$subtype[which(all_mark.2$subtype %in% c(6,10))] <- "B cells"
#all_mark.2$subtype[which(all_mark.2$subtype %in% c(2))] <- "pDC"
#all_mark.2$subtype[which(all_mark.2$subtype %in% c(9))] <- "mDC"
#all_mark.2$subtype[which(all_mark.2$subtype %in% c("mDC", "pDC"))] <- "Dendritic.cells"
all_mark.2$subtype[which(all_mark.2$subtype %in% c(2,9))] <- "APC"
all_mark.2$subtype[which(all_mark.2$subtype %in% c(11))] <- "Tregs"
all_mark.2$subtype[which(all_mark.2$subtype %in% c(12,15))] <- "CD4 Effector T-cells"
all_mark.2$subtype[which(all_mark.2$subtype %in% c(14))] <- "Neutrophils"
#note: this now has m1 and pDC interchanged as should be
#note 2: monocytes and other immune aren't immune cells --> remove from for.umap and add to tumor with the code at the start
all_mark.2$subtype[which(all_mark.2$subtype %in% c(4,7,13,16))] <- "Other.immune"

#for.umap, with removing of the other immune and monocytes that aren't immune
for.umap1 <- all_mark.2[all_mark.2$subtype != "Other.immune",]
for.umap1 <- for.umap1[for.umap1$subtype != "Monocytes",]
for.umap1 <- for.umap1[,c(1:9,12:17)] #removing CD57, CD45 and unique stuff


## unique stuff pt2 - immune cells
all_mark.2$core <- actual_markers.2$cores
all_mark.2$cellid <- actual_markers.2$Cellid

snowflake <- function(x,y) {
  paste(x, y, sep = "_")
}

all_mark.2$unique <- snowflake(all_mark.2$core, all_mark.2$cellid)

for (i in 1:nrow(all_mark.2)) {
  all_mark$subtype[which(all_mark$unique == all_mark.2$unique[i])] <- all_mark.2$subtype[i]
}


## Using all markers for nearest neighbor object:
pre.umap1 <- uwot::tumap(for.umap1[,c(2:13,15)], ret_nn = TRUE) #not using sample or label
pre.umap2 <- uwot::tumap(for.umap2[,c(2:9,12:17,19)], ret_nn = TRUE) #only markers (no CD45 or CD57) and gct


### No target metric (don't use this) ----
umap.run <- uwot::umap(for.umap2, nn_method = pre.umap2$nn, n_threads = 19, n_epochs = 5000, a=1.58, b=0.75)
umap.df2 <- data.frame(cbind(for.umap2, umap.run))

umap.run.i <- uwot::umap(for.umap1, nn_method = pre.umap1$nn, n_threads = 19, n_epochs = 5000, a=1.58, b=0.75)
umap.df.i <- data.frame(cbind(for.umap1, umap.run.i))


# Create globalcelltype column
umap.df2$GlobalCellType <- umap.df2$label
umap.df.i$ImmuneCellType <- umap.df.i$label
#umap.df2$GlobalCellType[which(umap.df2$GlobalCellType %in% c("Antigen.presenting.cells", "B.cells", "CD4.T.cells", "CD8.T.cells", "Macrophages", "Neutrophils", "NK.cells"))] <- "Immune.cells"
for.umap2$GlobalCellType <- umap.df2$GlobalCellType


plot.umap(umap.df2, '.id')
plot.umap(umap.df2[,1:28], 'cell.type')
plot.umap(umap.df2, 'cell.subtype')
plot.umap(umap.df2, 'GlobalCellType')
plot.umap(umap.df.i, 'ImmuneCellType')

#ds and no other mine
umap.df.noother <- umap.df2[umap.df2$label != "Other.cells",]
noother.ds3 <- umap.df.noother[sample(nrow(umap.df.noother), 10000), ]
plot.umap(noother.ds3, 'GlobalCellType')

umap.df.noother$label[umap.df.noother$label == "Epithelial.cells"] <- 1
umap.df.noother$label[umap.df.noother$label == "Immune.cells"] <- 2
umap.df.noother$label[umap.df.noother$label == "Stroma.cells"] <- 3
umap.df.noother$label <- as.factor(umap.df.noother$label)
noother.ds <- downSample(umap.df.noother, umap.df.noother$label)
noother.ds2 <- downSample(noother.ds, noother.ds$label)
plot.umap(noother.ds, 'GlobalCellType')


### Target metric: Global cell types (use this) ------
# tsi
#set.seed(1)
#ct.cols <- c("CD11c","CD15","CD163","CD1c","C20","CD31","CD3d","CD4","CD45","CD57","CD8a","CK7","Ecadherin","FOXP3","IBA1","vimentin")

umap.run.global <- uwot::umap(for.umap2[,c(2:9,12:17,19)], n_epochs = 30, y=as.factor(for.umap2$GlobalCellType)) # spread=0.7, min_dist = 0.001,
umap.df2.global <- data.frame(cbind(for.umap2, umap.run.global))

umap.df2.global.ds <- umap.df2.global[sample(nrow(umap.df2.global), 30000), ]

# old stuff, ignore
# umap.df2.global.ds$focus <- umap.df2.global.ds$subtype
# for (i in 1:nrow(umap.df2.global.ds)) {
#   if(umap.df2.global.ds$focus[i] == "Other.immune") {
#     umap.df2.global.ds$focus[i] <- "Other.immune"
#   }
#   else if(umap.df2.global.ds$focus[i] == "Monocytes") {
#     umap.df2.global.ds$focus[i] <- "Monocytes"
#   }
#   else {
#     umap.df2.global.ds$focus[i] <- "Others"
#   }
# }

# I'm not sure if these are needed
umap.df2.global.ds$subtype[umap.df2.global.ds$subtype == 7] <-"Stroma"
umap.df2.global.ds$subtype[umap.df2.global.ds$subtype == 16] <-"Endothelial"

setwd("C:/LocalData/nuplyyt/cyto/celltypes/images/ssumap_april/tsei")
plot.umap(umap.df2.global.ds, 'sample', save = T, "umap_patient.pdf") #better to run with no legend (edit the function)
plot.umap(umap.df2.global.ds, 'GlobalCellType', save=T, "umap_globaltypes.pdf")
plot.umap(umap.df2.global.ds, 'subtype', save=T, "umap_subtypes_tsei.pdf")
plot.umap(umap.df2.global.ds, 'global2', save=T, "umap_global_tsei.pdf") #only stroma?
# plot.umap(umap.df2.global.ds, 'focus', save=T, "umap_monoc_other.pdf") #old
umap.df2.global.ds$label <- as.character(umap.df2.global.ds$label)
plot.umap(umap.df2.global.ds, 'label', save=T, "umap_clusters.pdf") #should be run with 4 cols (edit the function)
plot.by.markers(umap.df2.global.ds, names(umap.df2.global.ds)[2:17], save=T, "umap_markers.pdf")

#plot.by.log(umap.df2.global, ct.cols)
#write.csv(umap.df2.global, "globalCelltypesDownsampled.csv", sep="\t", row.names = FALSE)


# subtypes

umap.run.st <- uwot::umap(for.umap1[,c(2:13,15)], n_epochs = 30, y=as.factor(for.umap1$subtype)) # spread=0.7, min_dist = 0.001,
umap.df.st <- data.frame(cbind(for.umap1, umap.run.st))
#umap.df.st.ds <- umap.df.st[sample(nrow(umap.df.st), 30000), ]

setwd("C:/LocalData/nuplyyt/cyto/celltypes/images/ssumap_april/immune")
plot.umap(umap.df.st, 'sample', save = T, "umap_patient_subtype2.pdf") #run with no legend
plot.umap(umap.df.st, 'subtype', save=T, "umap_subtypes2.pdf")
umap.df.st$label <- as.character(umap.df.st$label)
plot.umap(umap.df.st, 'label', save=T, "umap_clusters_subtype.pdf")
plot.by.markers(umap.df.st, names(umap.df.st)[2:13], save=T, "umap_markers_subtype.pdf")


#tsi but on subtype umap run
umap.df.global$GlobalCellType <- umap.df.global$label

umap.df.global$GlobalCellType[which(umap.df.global$GlobalCellType %in% c(3,4,6,9,10,12,13,15,19))] <- "Immune"
umap.df.global$GlobalCellType[which(umap.df.global$GlobalCellType %in% c(7,16))] <- "Stroma"
umap.df.global$GlobalCellType[which(umap.df.global$GlobalCellType %in% c(1,2,5,8,11,14,17,18))] <- "Tumor"

umap.df.global.ds <- umap.df.global[sample(nrow(umap.df.global), 30000), ]
plot.umap(umap.df.global.ds, 'GlobalCellType', save=T, "umap_globalcelltype.pdf")

### miscellaneous -----
# to manuscript (old julia's stuff)
mycolors <- c("#56B4E9","#E69F00","#999999")
df <- umap.df2.global
p <- ggplot(data.frame(df), aes_string(x="V1", y="V2", color="GlobalCellType")) +
  geom_point(size=0.6, stroke=0, alpha=0.2) +
  scale_color_manual(values=mycolors) +
  theme_bw() +
  coord_fixed(ratio=0.5) +
  #guides(colour = guide_legend(override.aes = list(size=5, shape=15, alpha=1))) +
  xlab("umap1") + ylab("umap2")
print(p)
ggsave("globalCelltypes_test.pdf")


# to cyto
immune.subset <- all_mark[all_mark$GlobalCellType == "Immune.cells",]
immune.subset <- immune.subset[,1:17]
names(immune.subset)[1] <- "Sample"
write.csv(immune.subset, "C:/LocalData/nuplyyt/cyto/celltypes/round_april/immune_subset_for_cyto.csv", row.names = F)

## for inga-maria
# tsi
all_mark
actual_markers <- read.csv("C:/LocalData/nuplyyt/cyto/celltypes/tma1_sample.csv")
tsi.w.markers <- cbind(actual_markers, all_mark$GlobalCellType, all_mark$subtype)
tumor.w.markers <- cbind(actual_markers, all_mark$subtype)
names(tsi.w.markers)[47] <- "GlobalCellType"
names(tsi.w.markers)[48] <- "Subtype"
names(tumor.w.markers)[47] <- "Subtype"
tumor.w.markers <- tumor.w.markers[tumor.w.markers$Subtype %in% c("Tumor", "Monocytes", "Other.immune"),]

immune.w.markers <- tsi.w.markers[tsi.w.markers$GlobalCellType == "Immune",]
immune.w.markers <- immune.w.markers[immune.w.markers$Subtype != "Monocytes",]
immune.w.markers <- immune.w.markers[immune.w.markers$Subtype != "Other.immune",]

write.csv(tsi.w.markers, "C:/LocalData/nuplyyt/cyto/celltypes/round_april/tsi_for_plots.csv") #124,623 cells
#write.csv(tsi.w.markers[tsi.w.markers$GlobalCellType == "Tumor",], "C:/LocalData/nuplyyt/cyto/celltypes/round_april/tsi_for_plots_tumor.csv") #80,265
#write.csv(tsi.w.markers[tsi.w.markers$GlobalCellType == "Stroma",], "C:/LocalData/nuplyyt/cyto/celltypes/round_april/tsi_for_plots_stroma.csv") #27,197
#write.csv(tsi.w.markers[tsi.w.markers$GlobalCellType == "Immune",], "C:/LocalData/nuplyyt/cyto/celltypes/round_april/tsi_for_plots_immune.csv") #16,514
write.csv(tumor.w.markers, "C:/LocalData/nuplyyt/cyto/celltypes/round_april/tsi_for_plots_tumor_new.csv") #85,621 cells
write.csv(immune.w.markers, "C:/LocalData/nuplyyt/cyto/celltypes/round_april/tsi_for_plots_immune_new.csv") #11,158 cells
#other 647

# stroma
actual_markers.s <- read.csv("C:/LocalData/nuplyyt/cyto/celltypes/round_april/tsi_for_plots_stroma.csv")
stroma.l <- all_mark[which(all_mark$GlobalCellType %in% c(7,16)),]
stroma.im <- cbind(actual_markers.s, stroma.l$GlobalCellType)
names(stroma.im)[49] <- "Stroma.subtype"
stroma.im$Stroma.subtype[which(stroma.im$Stroma.subtype %in% c(16))] <- "Endothelial.cells"
stroma.im$Stroma.subtype[which(stroma.im$Stroma.subtype %in% c(7))] <- "Non-endothelial.cells"

write.csv(stroma.im[stroma.im$Stroma.subtype == "Endothelial.cells",],
          "C:/LocalData/nuplyyt/cyto/celltypes/round_april/stroma_for_plots_endothelial.csv") #3,790
write.csv(stroma.im[stroma.im$Stroma.subtype == "Non-endothelial.cells",],
          "C:/LocalData/nuplyyt/cyto/celltypes/round_april/stroma_for_plots_nonendothelial.csv") #23,407

# immune
actual_markers.2 <- read.csv("C:/LocalData/nuplyyt/cyto/celltypes/round_april/tsi_for_plots_immune.csv")
immune.w.markers <- cbind(actual_markers.2, all_mark.2$subtype)
names(immune.w.markers)[49] <- "Subtype"

#these are the old subtypes with mDC and pDC separately
write.csv(immune.w.markers, "C:/LocalData/nuplyyt/cyto/celltypes/round_april/immune_for_plots.csv") #16,514
write.csv(immune.w.markers[immune.w.markers$Subtype == "CD8",], "C:/LocalData/nuplyyt/cyto/celltypes/round_april/immune_for_plots_cd8.csv") #1,570
write.csv(immune.w.markers[immune.w.markers$Subtype == "mDC",], "C:/LocalData/nuplyyt/cyto/celltypes/round_april/immune_for_plots_mdc.csv") #884
write.csv(immune.w.markers[immune.w.markers$Subtype == "CD4",], "C:/LocalData/nuplyyt/cyto/celltypes/round_april/immune_for_plots_cd4.csv") #778
write.csv(immune.w.markers[immune.w.markers$Subtype == "pDC",], "C:/LocalData/nuplyyt/cyto/celltypes/round_april/immune_for_plots_pdc.csv") #1,259
write.csv(immune.w.markers[immune.w.markers$Subtype == "M1.macrop",], "C:/LocalData/nuplyyt/cyto/celltypes/round_april/immune_for_plots_m1.csv") #1,774
write.csv(immune.w.markers[immune.w.markers$Subtype == "M2.macrop",], "C:/LocalData/nuplyyt/cyto/celltypes/round_april/immune_for_plots_m2.csv") #1,641
write.csv(immune.w.markers[immune.w.markers$Subtype == "T.regs",], "C:/LocalData/nuplyyt/cyto/celltypes/round_april/immune_for_plots_tregs.csv") #640
write.csv(immune.w.markers[immune.w.markers$Subtype == "B.cells",], "C:/LocalData/nuplyyt/cyto/celltypes/round_april/immune_for_plots_Bcells.csv") #2,351
write.csv(immune.w.markers[immune.w.markers$Subtype == "Monocytes",], "C:/LocalData/nuplyyt/cyto/celltypes/round_april/immune_for_plots_monocytes.csv") #2,043
write.csv(immune.w.markers[immune.w.markers$Subtype == "Neutrophils",], "C:/LocalData/nuplyyt/cyto/celltypes/round_april/immune_for_plots_neutrophils.csv") #261
#other 3,313

# tumor
tumor_markers <- read.csv("C:/LocalData/nuplyyt/cyto/celltypes/round_april/tsi_for_plots_tumor_new.csv")
#tumor_clusters <- read.delim("C:/LocalData/nuplyyt/cyto/celltypes/round_april/tumor_Cyto_clustering_phenograph40_noCK7.csv")
tumor_clusters <- read.delim("C:/LocalData/nuplyyt/cyto/celltypes/round_april/tumor_Cyto_clustering_kmeans25_noCK7.csv")
tumor.final <- cbind(tumor_markers, tumor_clusters$label)
names(tumor.final)[49] <- "Cluster"
tumor.final$Metacluster <- tumor.final$Cluster
tumor.final$Metacluster2 <- tumor.final$Cluster

#phenograph 40
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(11))] <- "VIM+ECAD+Ki67+"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(5))]  <- "VIM+ECAD+Ki67-"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(16))] <- "VIM+ECAD-Ki67+"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(19))] <- "VIM+ECAD-Ki67-"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(1))]    <- "VIM=ECAD+Ki67+"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(8))]    <- "VIM=ECAD+Ki67="
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(9,12))] <- "VIM=ECAD+Ki67-"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(13))]   <- "VIM=ECAD+Ki67-cCasp+P21+PSTAT1+"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(7))]    <- "VIM=ECAD=Ki67-"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(20))]   <- "VIM=ECAD-Ki67-"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(4,10,14))]   <- "VIM-ECAD+Ki67+"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(2))]         <- "VIM-ECAD+Ki67="
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(18))]        <- "VIM-ECAD+Ki67-yH2AX+"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(3,6,15,17))] <- "VIM-ECAD+Ki67-"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(22))]        <- "VIM-ECAD=Ki67=PDL1+"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(21))]        <- "VIM-ECAD-Ki67-"
# 
# #k-means 25
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(13,14))] <- "VIM+Ki67+"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(4))]     <- "VIM+Ki67+ECAD-"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(10))]    <- "VIM+Ki67-"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(1,20))]  <- "VIM+Ki67-ECAD-"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(17))]   <- "VIM=Ki67+"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(3))]    <- "VIM=Ki67="
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(2,21))] <- "VIM=Ki67-"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(5,7,8,12,15,18))] <- "VIM-Ki67+"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(9,11,24))]        <- "VIM-Ki67="
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(16))]             <- "VIM-Ki67-ECAD-"
# tumor.final$Metacluster[which(tumor.final$Cluster %in% c(6,19,22,23,25))]  <- "VIM-Ki67-"

#k-means 25, MST
tumor.final$Metacluster[which(tumor.final$Cluster %in% c(3,13,14,17))] <- "Proliferating EMT"
tumor.final$Metacluster[which(tumor.final$Cluster %in% c(2,10,21))] <- "EMT"
tumor.final$Metacluster[which(tumor.final$Cluster %in% c(1,4,20))] <- "Mesenchymal"
tumor.final$Metacluster[which(tumor.final$Cluster %in% c(6,22,23,25))] <- "Epithelial"
tumor.final$Metacluster[which(tumor.final$Cluster %in% c(7,8,9,12,15,18,24))] <- "Proliferating epithelial"
tumor.final$Metacluster[which(tumor.final$Cluster %in% c(5,11))] <- "Apoptotic"
tumor.final$Metacluster[which(tumor.final$Cluster %in% c(19))] <- "Hyperfunctional epithelial" #New name is Functional epithelial
tumor.final$Metacluster[which(tumor.final$Cluster %in% c(16))] <- "Negative" #New name is other

#k-means 25, MST, only 3 metaclusters
tumor.final$Metacluster2[which(tumor.final$Cluster %in% c(2,3,10,13,14,17,21))] <- "EMT"
tumor.final$Metacluster2[which(tumor.final$Cluster %in% c(1,4,20))] <- "Mesenchymal"
tumor.final$Metacluster2[which(tumor.final$Cluster %in% c(5,6,7,8,9,11,12,15,16,18,19,22,23,24,25))] <- "Epithelial"


write.csv(tumor.final, "C:/LocalData/nuplyyt/cyto/celltypes/round_april/tumor_for_plots_metacluster_noCK7.csv")
write.csv(tumor.final, "C:/LocalData/nuplyyt/cyto/celltypes/round_april/tumor_for_plots_metacluster_noCK7_kmeans.csv")
write.csv(tumor.final, "C:/LocalData/nuplyyt/cyto/celltypes/round_april/tumor_for_plots_metacluster_noCK7_kmeans_MST.csv")
write.csv(tumor.final, "C:/LocalData/nuplyyt/cyto/celltypes/round_april/tumor_for_plots_metacluster_noCK7_kmeans_MST_splitEMT.csv")