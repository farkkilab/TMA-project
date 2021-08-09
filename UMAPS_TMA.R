
############### DATA ############################################################################################
library(viridis)
library(ggplot2)


all_celltypes <- read.csv("TMA_annotated_single_cell_data.csv")


#################################################################################################################
################   GLOBAL CELL TYPE UMAPS #######################################################################

df <- all_celltypes[, c("Sample","CD163","CD20","CD4","CD3d","CD8a", "CD45", "pSTAT1", "P21", "Ki67", "PDL1", "PD1", "cCasp3",
                                 "FOXP3","IBA1","CK7","CD11c","Ecadherin","vimentin", "CD31","Subtype", "GlobalCellType")]


df$Subtype <- as.character(df$Subtype)

df[which(df$GlobalCellType == "Endothelia"), "Subtype"] <- "Endothelia"
df[which(df$Subtype == "Stroma_Endothelia"), "Subtype"] <- "Stroma"

df$Subtype <- as.factor(df$Subtype)
#2:14
cell_type <- df$Subtype
umap_s = uwot::umap(df[,c(2:20)],n_neighbors = 80, scale=T ,spread=1, min_dist = 0.5,n_epochs = 50)
uplot = data.frame(x = umap_s[,1], y= umap_s[,2])
cell_type <- df$Subtype
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(4)
p = ggplot(uplot,aes(x,y, color=cell_type))+ geom_point(size=0.3, stroke=0, alpha=0.7)+ 
  scale_fill_viridis(discrete = TRUE) + 
  theme(legend.title = element_blank())+ guides(color = guide_legend(override.aes = list(size=5))) + 
  scale_color_manual(values=mycolors) + theme_bw() + coord_fixed(ratio=1)+
  guides(colour = guide_legend(override.aes = list(size=5, shape=15, alpha=1), direction="horizontal", ncol=1, label.position="bottom", byrow=F)) + 
  xlab("umap1") + ylab("umap2")+ theme(legend.position="none") + ylim(-10, 10) + theme(aspect.ratio=1)
print(p)



####################### COLORED BY MARKER ####################################################################

for (i in c(2:20)){
  markers <- df[, i]
  p = ggplot(uplot,aes(x,y, color=markers))+ geom_point(size=0.3, stroke=0, alpha=0.7)+ 
    guides(color = guide_legend(override.aes = list(size=5))) + 
    theme_bw() + coord_fixed(ratio=1)+ 
    xlab("umap1") + ylab("umap2") + ggtitle(colnames(df)[i])
  print(p)
  
}

####################### COLORED BY PATIENT OF ORIGIN #########################################################
Patient <- as.factor(df$Sample)
n <- 44
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
mycolors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
p1 = ggplot(uplot,aes(x,y, color=Patient))+ geom_point(size=0.3, stroke=0, alpha=0.7)+ 
  scale_fill_viridis(discrete = TRUE) + 
  theme(legend.title = element_blank())+ guides(color = guide_legend(override.aes = list(size=5))) + 
  theme_bw() + coord_fixed(ratio=1)+scale_color_manual(values=mycolors) +
  guides(colour = guide_legend(override.aes = list(size=5, shape=15, alpha=1), direction="horizontal", ncol=3, label.position="bottom", byrow=F)) + 
  xlab("umap1") + ylab("umap2") + theme(legend.position="none", aspect.ratio=1) + xlim(-10, 10) + ylim(-10, 10)
print(p1)


##############################################################################################################
##################### IMMUNE CELL UMAPS ######################################################################

all_celltypes_24092020 <- read.csv("TMA_annotated_single_cell_data.csv")

all_celltypes_24092020$GlobalCellType <- as.character(all_celltypes_24092020$GlobalCellType)

all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "CD11c+CD163+IBA1+Macrophages"), "GlobalCellType"] <- "IBA1+CD163+CD11c+ Macrophages"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "CD11c+IBA1+Macrophages"), "GlobalCellType"] <- "IBA1+CD11c+ Macrophages"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "CD163+IBA1+Macrophages"), "GlobalCellType"] <- "IBA1+CD163+ Macrophages"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "CD4"), "GlobalCellType"] <- "CD4+ Effector T-cells"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "CD8"), "GlobalCellType"] <- "CD8+ T-cells"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "CD163+Macrophages"), "GlobalCellType"] <- "CD163+ Macrophages"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "IBA1+Macrophages"), "GlobalCellType"] <- "IBA1+ Macrophages"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "CD11c+APC"), "GlobalCellType"] <- "CD11c+APC"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "FOXP3+CD4+Tregs"), "GlobalCellType"] <- "FOXP3+CD4+ T-regs"


all_cell_types_16092020_immune <- all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Immune"),]

df <- all_cell_types_16092020_immune[, c("Sample","CD163","CD20","CD4","CD3d","CD8a",
                                         "FOXP3","IBA1","CD11c", "CD1c", "GlobalCellType")]
df$GlobalCellType <- as.factor(df$GlobalCellType)


cell_type <- df$GlobalCellType

umap_s = uwot::umap(df[,c(2:10)],n_neighbors = 80, scale=T ,spread=2.5,n_epochs = 60, min_dist = 0.2)
uplot = data.frame(x = umap_s[,1], y= umap_s[,2])
mycolors <- colorRampPalette(brewer.pal(10, "Paired"))(10)
p2 = ggplot(uplot,aes(x,y, color=cell_type))+ geom_point(size=0.7, stroke=0, alpha=0.7)+ 
  scale_fill_viridis(discrete = TRUE) + 
  theme(legend.title = element_blank())+ guides(color = guide_legend(override.aes = list(size=5))) + 
  scale_color_manual(values=mycolors) + theme_bw() + coord_fixed(ratio=1)+
  guides(colour = guide_legend(override.aes = list(size=5, shape=15, alpha=1), direction="horizontal", ncol=1, label.position="right", byrow=F, reverse=TRUE)) + 
  xlab("umap1") + ylab("umap2")+ theme(legend.position="none", aspect.ratio=1)
print(p2)

############################ COLORED BY MARKER ################################################################

for (i in c(2:10)){
  Expression <- df[, i]
  p = ggplot(uplot,aes(x,y, color=Expression))+ geom_point(size=0.5, stroke=0, alpha=0.7)+  
    theme_bw() + coord_fixed(ratio=1)+ 
    xlab("umap1") + ylab("umap2") + ggtitle(colnames(df)[i])
  print(p)
}


################################################################################################################
####################### STROMAL CELL UMAPS #####################################################################

all_celltypes_24092020 <- read.csv("TMA_annotated_single_cell_data.csv")
all_celltypes_24092020$GlobalCellType <- as.character(all_celltypes_24092020$GlobalCellType)

all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "High-PDL1"), "GlobalCellType"] <- "Functional stroma"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "Non-proliferative_Stroma"), "GlobalCellType"] <- "Non-proliferative Stroma"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "Low_eccentricity_medium_vimentin"), "GlobalCellType"] <- "Low eccentricity"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "Proliferative_Stroma"), "GlobalCellType"] <- "Proliferative Stroma"
all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "High-proliferative_Stroma"), "GlobalCellType"] <- "High-proliferative Stroma"

all_cell_types_13102020_stromal <- all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Stroma_Endothelia"),]

df <- all_cell_types_13102020_stromal[, c("Sample","cCasp3", "pSTAT1", "Ki67",
                                          "PDL1","vimentin","CD31","P21","Eccentricity", "GlobalCellType")]
df$GlobalCellType <- as.factor(df$GlobalCellType)


cell_type <- df$GlobalCellType
umap_s = uwot::umap(df[,c(2:9)],n_neighbors = 100, scale=T ,spread=1.5,n_epochs = 50)
uplot = data.frame(x = umap_s[,1], y= umap_s[,2])
mycolors <- colorRampPalette(brewer.pal(9, "Set3"))(9)
p3 = ggplot(uplot,aes(x,y, color=cell_type))+ geom_point(size=0.8, stroke=0, alpha=0.7)+ 
  scale_fill_viridis(discrete = TRUE) + 
  theme(legend.title = element_blank())+ guides(color = guide_legend(override.aes = list(size=5))) + 
  scale_color_manual(values=mycolors) + theme_bw() + coord_fixed(ratio=1)+
  guides(colour = guide_legend(override.aes = list(size=5, shape=15, alpha=1), direction="horizontal", ncol=1, label.position="bottom", byrow=F, reverse=TRUE)) + 
  xlab("umap1") + ylab("umap2") + theme(legend.position="none") + xlim(-10, 10) + ylim(-10, 10)
print(p3)

################### COLORED BY MARKER ############################################################################


for (i in c(2:9)){
  expression <- df[, i]
  p = ggplot(uplot,aes(x,y, color=expression))+ geom_point(size=0.8, stroke=0, alpha=0.7)+
    theme_bw() + coord_fixed(ratio=1)+ 
    xlab("umap1") + ylab("umap2") + ggtitle(colnames(df)[i])
  print(p)
  
}


##################################################################################################################
########################## TUMOR CELL UMAPS ######################################################################

read.csv("TMA_annotated_single_cell_data.csv")
all_celltypes_24092020$GlobalCellType <- as.character(all_celltypes_24092020$GlobalCellType)

all_celltypes_24092020[which(all_celltypes_24092020$GlobalCellType == "Hyperfunctional epithelial"), "GlobalCellType"] <- "Functional epithelial"

all_cell_types_13102020_tumor <- all_celltypes_24092020[which(all_celltypes_24092020$Subtype == "Tumor"),]
all_cell_types_13102020_tumor <- all_cell_types_13102020_tumor[-which(all_cell_types_13102020_tumor$GlobalCellType == "Negative"),]

df <- all_cell_types_13102020_tumor[, c("Sample","cCasp3", "pSTAT1", "Ki67",
                                        "PDL1","vimentin","Ecadherin","P21", "GlobalCellType")]
df$GlobalCellType <- as.factor(df$GlobalCellType)


cell_type <- df$GlobalCellType

umap_s = uwot::umap(df[,c(2:8)],n_neighbors = 80, scale=T ,spread=0.7,n_epochs = 30)
uplot = data.frame(x = umap_s[,1], y= umap_s[,2])
mycolors <- rev(colorRampPalette(brewer.pal(7, "Accent"))(7))
p4 = ggplot(uplot,aes(x,y, color=cell_type))+ geom_point(size=0.3, stroke=0, alpha=0.7)+ 
  scale_fill_viridis(discrete = TRUE) + 
  theme(legend.title = element_blank())+ guides(color = guide_legend(override.aes = list(size=5))) + 
  scale_color_manual(values=mycolors) + theme_bw() + coord_fixed(ratio=1)+
  guides(colour = guide_legend(override.aes = list(size=5, shape=15, alpha=1), direction="horizontal", ncol=1, label.position="bottom", byrow=F, reverse=TRUE)) + 
  xlab("umap1") + ylab("umap2") + theme(aspect.ratio=1, legend.position="none") 
print(p4)


############################################### COLORED BY EXPRESSION ######################################

for (i in c(2:8)){
  expression <- df[, i]
  p = ggplot(uplot,aes(x,y, color=expression))+ geom_point(size=0.7, stroke=0, alpha=0.7)+
    theme_bw() + coord_fixed(ratio=1)+ 
    xlab("umap1") + ylab("umap2") + ggtitle(colnames(df)[i])
  print(p)
}





