library(survminer)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scater)
library(rmarkdown)
require("survival")

tma <- read.csv("C:/LocalData/nuplyyt/cyto/celltypes/images/maybe_paper/all_data_TMA_surv.csv")
tma <- tma[(tma$Identifier %in% counts$Sample),]
all_immune <- read.csv("C:/LocalData/nuplyyt/cyto/celltypes/round_april/tsi_for_plots_immune_new.csv")
all_immune$Subtype <- as.character(all_immune$Subtype)
all_immune$Subtype[all_immune$Subtype=="mDC"] <- "APC"
all_immune$Subtype[all_immune$Subtype=="pDC"] <- "APC"
tumor <- read.csv("C:/LocalData/nuplyyt/cyto/celltypes/round_april/tumor_for_plots_metacluster_noCK7_kmeans_MST_splitEMT.csv")
tumor$Metacluster <- as.character(tumor$Metacluster)
tumor$Metacluster[tumor$Metacluster=="Hyperfunctional epithelial"] <- "Functional epithelial"
tumor <- tumor[tumor$Metacluster!="Negative",]


### functional markers ----


# funct + immune + PFI
for (n in unique(all_immune$Subtype)){
  
  subtype <- all_immune[which(all_immune$Subtype == n),]
  meex <- subtype %>% group_by(Sample) %>% summarise(median_Ki67 = median(Ki67),
                                                     median_pSTAT1 = median(pSTAT1),
                                                     median_P21 = median(P21),
                                                     median_PD1=median(PD1),
                                                     median_PDL1=median(PDL1),
                                                     median_cCasp3=median(cCasp3))
  
  meex <- merge(meex, tma[,c(1,5,20,32)], by.x="Sample", by.y="Identifier")
  
  distlist <- list()
  plotlist <- list()
  
  for (i in seq_along(colnames(meex))[2:7]){
    #p_hist <- hist(median_expressions[, i], breaks = 50)
    #p_dist <- ggplot(meex, aes_string(x=colnames(meex)[i], fill="Type")) + 
    #  geom_density(alpha=0.2)+scale_fill_brewer(palette="Dark2")
    #distlist[[i-1]] <- p_dist
    
    tma.cat <- as.data.frame(cbind(meex$PFI_time, meex$status_PFI, meex[[i]], meex$Type))
    names(tma.cat) <- c("PFI_time","Status_PFI","Marker","HR")
    med <- median(tma.cat$Marker)
    
    tma.cat$Marker[which(tma.cat$Marker >= med)] <- "high"
    tma.cat$Marker[which(tma.cat$Marker < med)] <- "low"
    
    #BRCA survival curve
    tma.cat.brca <- tma.cat[which(tma.cat$HR == 1),]
    formula <- paste0("fit <- survfit(Surv(PFI_time,Status_PFI) ~ ",colnames(tma.cat.brca)[3],
                      ", data = tma.cat.brca)")
    eval(parse(text=formula))
    pb <- ggsurvplot(fit, data=tma.cat.brca, cencor=T, pval=T, pval.coord=c(2500, 0.75),
                     title=paste("BRCA1/2", n, colnames(meex)[i],"PFI"), risk.table=T,
                     risk.table.height = 0.3, xlim = c(0, 3000), break.time.by=500, palette=c("red","blue"))
    plotlist[[(i-1)*2-1]] <- pb
    #print(pb)
    
    #HR survival curve
    tma.cat.hr <- tma.cat[which(tma.cat$HR == 2),]
    formula <- paste0("fit <- survfit(Surv(PFI_time,Status_PFI) ~ ",colnames(tma.cat.hr)[3],
                      ", data = tma.cat.hr)")
    eval(parse(text=formula))
    ph <- ggsurvplot(fit, data=tma.cat.hr, cencor=T, pval=T, pval.coord=c(750, 0.75),
                     title=paste("HRwt", n, colnames(meex)[i],"PFI"), risk.table=T,
                     risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500, palette=c("red","blue"))
    plotlist[[(i-1)*2]] <- ph
    #print(ph)
    
    
    #all survival curve
    formula <- paste0("fit <- survfit(Surv(PFI_time,Status_PFI) ~ ",colnames(tma.cat)[3],
                      ", data = tma.cat)")
    eval(parse(text=formula))
    pa <- ggsurvplot(fit, data=tma.cat, cencor=T, pval=T, pval.coord=c(750, 0.75),
                     title=paste("All patients, Expression of", colnames(meex)[i], "in", n, "and PFI"), risk.table=T,
                     risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500, palette=c("red","blue"))
    #plotlist[[i*2]] <- pa
  }
  
  pdf(paste0("survival_funct_immune_PFI_", n, "_redblue_PDL1.pdf"))
  multiplot(plotlist = plotlist)
  dev.off()
  
  #arrange_ggsurvplots(plotlist)
  
}

# funct + immune + OS
for (n in unique(all_immune$Subtype)){
  
  subtype <- all_immune[which(all_immune$Subtype == n),]
  meex <- subtype %>% group_by(Sample) %>% summarise(median_Ki67 = median(Ki67),
                                                     median_pSTAT1 = median(pSTAT1),
                                                     median_P21 = median(P21),
                                                     median_PD1=median(PD1),
                                                     median_PDL1=median(PDL1),
                                                     median_cCasp3=median(cCasp3))
  
  meex <- merge(meex, tma[,c(1,5,27,33)], by.x="Sample", by.y="Identifier")
  
  #distlist <- list()
  #plotlist <- list()
  
  for (i in seq_along(colnames(meex))[2:7]){
    #p_hist <- hist(median_expressions[, i], breaks = 50)
    #p_dist <- ggplot(meex, aes_string(x=colnames(meex)[i], fill="Type")) + 
    #  geom_density(alpha=0.2)+scale_fill_brewer(palette="Dark2")
    #distlist[[i-1]] <- p_dist
    
    tma.cat <- as.data.frame(cbind(meex$OS_time, meex$status_OS, meex[[i]], meex$Type))
    names(tma.cat) <- c("OS_time","Status_OS","Marker","HR")
    med <- median(tma.cat$Marker)
    
    tma.cat$Marker[which(tma.cat$Marker >= med)] <- "high"
    tma.cat$Marker[which(tma.cat$Marker < med)] <- "low"
    
    #BRCA survival curve
    tma.cat.brca <- tma.cat[which(tma.cat$HR == 1),]
    formula <- paste0("fit <- survfit(Surv(OS_time,Status_OS) ~ ",colnames(tma.cat.brca)[3],
                      ", data = tma.cat.brca)")
    eval(parse(text=formula))
    pb <- ggsurvplot(fit, data=tma.cat.brca, cencor=T, pval=T, pval.coord=c(2500, 1),
                     title=paste("BRCA1/2", n, colnames(meex)[i],"OS"), risk.table=T,
                     risk.table.height = 0.3, xlim = c(0, 3000), break.time.by=500, palette=c("red","blue"))
    plotlist[[(i-1)*2-1]] <- pb
    #print(pb)
    
    #HR survival curve
    tma.cat.hr <- tma.cat[which(tma.cat$HR == 2),]
    formula <- paste0("fit <- survfit(Surv(OS_time,Status_OS) ~ ",colnames(tma.cat.hr)[3],
                      ", data = tma.cat.hr)")
    eval(parse(text=formula))
    ph <- ggsurvplot(fit, data=tma.cat.hr, cencor=T, pval=T, pval.coord=c(750, 0.75),
                     title=paste("HRwt", n, colnames(meex)[i],"OS"), risk.table=T,
                     risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500, palette=c("red","blue"))
    plotlist[[(i-1)*2]] <- ph
    #print(ph)
    
    
    #all survival curve
    formula <- paste0("fit <- survfit(Surv(OS_time,Status_OS) ~ ",colnames(tma.cat)[3],
                      ", data = tma.cat)")
    eval(parse(text=formula))
    pa <- ggsurvplot(fit, data=tma.cat, cencor=T, pval=T, pval.coord=c(750, 0.75),
                     title=paste("All patients, Expression of", colnames(meex)[i], "in", n, "and OS"), risk.table=T,
                     risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500, palette=c("red","blue"))
    #plotlist[[i*2]] <- pa
  }
  
  pdf(paste0("survival_funct_immune_OS_", n, "_redblue.pdf"))
  multiplot(plotlist = plotlist)
  dev.off()
  
  #arrange_ggsurvplots(plotlist)
  
}

# funct + tumor + PFI
for (n in unique(tumor$Metacluster)){
  
  subtype <- tumor[which(tumor$Metacluster == n),]
  meex <- subtype %>% group_by(Sample) %>% summarise(median_Ki67 = median(Ki67),
                                                     median_pSTAT1 = median(pSTAT1),
                                                     median_P21 = median(P21),
                                                     median_PD1=median(PD1),
                                                     median_PDL1=median(PDL1),
                                                     median_cCasp3=median(cCasp3))
  
  meex <- merge(meex, tma[,c(1,5,20,32)], by.x="Sample", by.y="Identifier")
  
  #distlist <- list()
  #plotlist <- list()
  
  for (i in seq_along(colnames(meex))[2:7]){
    
    #p_dist <- ggplot(meex, aes_string(x=colnames(meex)[i], fill="Type")) + 
    #  geom_density(alpha=0.2)+scale_fill_brewer(palette="Dark2")
    #distlist[[i-1]] <- p_dist
    
    tma.cat <- as.data.frame(cbind(meex$PFI_time, meex$status_PFI, meex[[i]], meex$Type))
    names(tma.cat) <- c("PFI_time","Status_PFI","Marker","HR")
    med <- median(tma.cat$Marker)
    
    tma.cat$Marker[which(tma.cat$Marker >= med)] <- "high"
    tma.cat$Marker[which(tma.cat$Marker < med)] <- "low"
    
    #BRCA survival curve
    tma.cat.brca <- tma.cat[which(tma.cat$HR == 1),]
    formula <- paste0("fit <- survfit(Surv(PFI_time,Status_PFI) ~ ",colnames(tma.cat.brca)[3],
                      ", data = tma.cat.brca)")
    eval(parse(text=formula))
    pb <- ggsurvplot(fit, data=tma.cat.brca, cencor=T, pval=T, pval.coord=c(2500, 0.75),
                     title=paste("BRCA1/2", n, colnames(meex)[i],"PFI"), risk.table=T,
                     risk.table.height = 0.3, xlim = c(0, 3000), break.time.by=500)
    plotlist[[(i-1)*2-1]] <- pb
    #print(pb)
    
    #HR survival curve
    tma.cat.hr <- tma.cat[which(tma.cat$HR == 2),]
    formula <- paste0("fit <- survfit(Surv(PFI_time,Status_PFI) ~ ",colnames(tma.cat.hr)[3],
                      ", data = tma.cat.hr)")
    eval(parse(text=formula))
    ph <- ggsurvplot(fit, data=tma.cat.hr, cencor=T, pval=T, pval.coord=c(750, 0.75),
                     title=paste("HRwt", n, colnames(meex)[i],"PFI"), risk.table=T,
                     risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500)
    plotlist[[(i-1)*2]] <- ph
    #print(ph)
    
    # all survival curve
    formula <- paste0("fit <- survfit(Surv(PFI_time,Status_PFI) ~ ",colnames(tma.cat)[3],
                      ", data = tma.cat)")
    eval(parse(text=formula))
    pa <- ggsurvplot(fit, data=tma.cat, cencor=T, pval=T, pval.coord=c(750, 0.75),
                     title=paste("All patients, Expression of", colnames(meex)[i], "in", n, "and PFI"), risk.table=T,
                     risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500, palette=c("red","blue"))
  }
  
  pdf(paste0("survival_funct_tumor_PFI_", n, ".pdf"))
  multiplot(plotlist = plotlist)
  dev.off()
  
  #arrange_ggsurvplots(plotlist)
  
}

# funct + tumor + OS
for (n in unique(tumor$Metacluster)){
  
  subtype <- tumor[which(tumor$Metacluster == n),]
  meex <- subtype %>% group_by(Sample) %>% summarise(median_Ki67 = median(Ki67),
                                                     median_pSTAT1 = median(pSTAT1),
                                                     median_P21 = median(P21),
                                                     median_PD1=median(PD1),
                                                     median_PDL1=median(PDL1),
                                                     median_cCasp3=median(cCasp3))
  
  meex <- merge(meex, tma[,c(1,5,27,33)], by.x="Sample", by.y="Identifier")
  
  #distlist <- list()
  #plotlist <- list()
  
  for (i in seq_along(colnames(meex))[2:7]){
    #p_dist <- ggplot(meex, aes_string(x=colnames(meex)[i], fill="Type")) + 
    #  geom_density(alpha=0.2)+scale_fill_brewer(palette="Dark2")
    #distlist[[i-1]] <- p_dist
    
    tma.cat <- as.data.frame(cbind(meex$OS_time, meex$status_OS, meex[[i]], meex$Type))
    names(tma.cat) <- c("OS_time","Status_OS","Marker","HR")
    med <- median(tma.cat$Marker)
    
    tma.cat$Marker[which(tma.cat$Marker >= med)] <- "high"
    tma.cat$Marker[which(tma.cat$Marker < med)] <- "low"
    
    #BRCA survival curve
    tma.cat.brca <- tma.cat[which(tma.cat$HR == 1),]
    formula <- paste0("fit <- survfit(Surv(OS_time,Status_OS) ~ ",colnames(tma.cat.brca)[3],
                      ", data = tma.cat.brca)")
    eval(parse(text=formula))
    pb <- ggsurvplot(fit, data=tma.cat.brca, cencor=T, pval=T, pval.coord=c(2500, 0.75),
                     title=paste("BRCA1/2", n, colnames(meex)[i],"OS"), risk.table=T,
                     risk.table.height = 0.3, xlim = c(0, 3000), break.time.by=500)
    plotlist[[(i-1)*2-1]] <- pb
    #print(pb)
    
    #HR survival curve
    tma.cat.hr <- tma.cat[which(tma.cat$HR == 2),]
    formula <- paste0("fit <- survfit(Surv(OS_time,Status_OS) ~ ",colnames(tma.cat.hr)[3],
                      ", data = tma.cat.hr)")
    eval(parse(text=formula))
    ph <- ggsurvplot(fit, data=tma.cat.hr, cencor=T, pval=T, pval.coord=c(750, 0.75),
                     title=paste("HRwt", n, colnames(meex)[i],"OS"), risk.table=T,
                     risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500)
    plotlist[[(i-1)*2]] <- ph
    #print(ph)
    
    
    # all survival curve
    formula <- paste0("fit <- survfit(Surv(OS_time,Status_OS) ~ ",colnames(tma.cat)[3],
                      ", data = tma.cat)")
    eval(parse(text=formula))
    pa <- ggsurvplot(fit, data=tma.cat, cencor=T, pval=T, pval.coord=c(750, 0.75),
                     title=paste("All patients, Expression of", colnames(meex)[i], "in", n, "and OS"), risk.table=T,
                     risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500, palette=c("red","blue"))
    
    
  }
  
  pdf(paste0("survival_funct_tumor_OS_", n, ".pdf"))
  multiplot(plotlist = plotlist)
  dev.off()
  
  #arrange_ggsurvplots(plotlist)
  
}


### proportions ----

# prop + immune + PFI
prop.i <- data.frame("Sample"=1:44,"CD8"=rep(0,44),"M2.macrop"=rep(0,44),"B.cells"=rep(0,44),
                     "APC"=rep(0,44),"M1.macrop"=rep(0,44),"T.regs"=rep(0,44),
                     "Neutrophils"=rep(0,44),"CD4"=rep(0,44))
prop.i$Sample <- unique(all_immune$Sample)
for (i in unique(all_immune$Sample)){
  patient <- all_immune[all_immune$Sample==i,]
  row <- which(prop.i$Sample==i)
  
  for (celltype in names(prop.i)[2:9]){
    prop.i[row, celltype] <- nrow(patient[patient$Subtype == celltype,])/nrow(patient)
  }
}
prop.i <- merge(prop.i, tma[,c(1,5,20,32)], by.x="Sample", by.y="Identifier")
plotlist <- list()
#here you need to check the plotlist[] if you want all proportions in the same image
for (i in c(1:8)){
  
  #  p_dist <- ggplot(meex, aes_string(x=colnames(meex)[i], fill="Type")) + 
  #    geom_density(alpha=0.2)+scale_fill_brewer(palette="Dark2")
  #  distlist[[i-1]] <- p_dist
    
    tma.cat <- as.data.frame(cbind(prop.i$PFI_time, prop.i$status_PFI, prop.i[[i+1]], prop.i$Type))
    names(tma.cat) <- c("PFI_time","Status_PFI","Subtype","HR")
    med <- median(tma.cat$Subtype)
    
    tma.cat$Subtype[which(tma.cat$Subtype >= med)] <- "high"
    tma.cat$Subtype[which(tma.cat$Subtype < med)] <- "low"
    
    #BRCA survival curve
    tma.cat.brca <- tma.cat[which(tma.cat$HR == 1),]
    formula <- paste0("fit <- survfit(Surv(PFI_time,Status_PFI) ~ ",colnames(tma.cat.brca)[3],
                     ", data = tma.cat.brca)")
    eval(parse(text=formula))
    pb <- ggsurvplot(fit, data=tma.cat.brca, cencor=T, pval=T, pval.coord=c(2500, 0.75),
                     title=paste("BRCA1/2 Proportion of", colnames(prop.i)[i+1],"PFI"), risk.table=T,
                     risk.table.height = 0.3, xlim = c(0, 3000), break.time.by=500, palette=c("red","blue"))
    plotlist[[1]] <- pb
    #plotlist[[i*2-1]] <- pb
    #print(pb)

    #HR survival curve
    tma.cat.hr <- tma.cat[which(tma.cat$HR == 2),]
    formula <- paste0("fit <- survfit(Surv(PFI_time,Status_PFI) ~ ",colnames(tma.cat.hr)[3],
                      ", data = tma.cat.hr)")
    eval(parse(text=formula))
    ph <- ggsurvplot(fit, data=tma.cat.hr, cencor=T, pval=T, pval.coord=c(750, 0.75),
                     title=paste("HRwt Proportion of", colnames(prop.i)[i+1],"PFI"), risk.table=T,
                     risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500, palette=c("red","blue"))
    plotlist[[2]] <- ph
    #plotlist[[i*2]] <- ph
    #print(ph)
  
    
    #all survival curve
    formula <- paste0("fit <- survfit(Surv(PFI_time,Status_PFI) ~ ",colnames(tma.cat)[3],
                      ", data = tma.cat)")
    eval(parse(text=formula))
    pa <- ggsurvplot(fit, data=tma.cat, cencor=T, pval=T, pval.coord=c(750, 0.75),
                     title=paste("All patients, Proportion of", colnames(prop.i)[i+1],"and PFI"), risk.table=T,
                     risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500, palette=c("red","blue"))
    plotlist[[3]] <- pa
    
    
  #arrange_ggsurvplots(plotlist)
  
}
#setwd("C:/LocalData/nuplyyt/cyto/celltypes/images/maybe_paper/Survival")
pdf("survival_prop_immune_PFI_redblue.pdf")
multiplot(plotlist = plotlist)
dev.off()


# prop + immune + OS
prop.i <- data.frame("Sample"=1:44,"CD8"=rep(0,44),"M2.macrop"=rep(0,44),"B.cells"=rep(0,44),
                     "APC"=rep(0,44),"M1.macrop"=rep(0,44),"T.regs"=rep(0,44),
                     "Neutrophils"=rep(0,44),"CD4"=rep(0,44))
prop.i$Sample <- unique(all_immune$Sample)
for (i in unique(all_immune$Sample)){
  patient <- all_immune[all_immune$Sample==i,]
  row <- which(prop.i$Sample==i)
  
  for (celltype in names(prop.i)[2:9]){
    prop.i[row, celltype] <- nrow(patient[patient$Subtype == celltype,])/nrow(patient)
  }
}
prop.i <- merge(prop.i, tma[,c(1,5,27,33)], by.x="Sample", by.y="Identifier")
plotlist <- list()
#here also you need to edit some if you want some other than CD8
for (i in c(1:8)){
  
  #  p_dist <- ggplot(meex, aes_string(x=colnames(meex)[i], fill="Type")) + 
  #    geom_density(alpha=0.2)+scale_fill_brewer(palette="Dark2")
  #  distlist[[i-1]] <- p_dist
  
  tma.cat <- as.data.frame(cbind(prop.i$OS_time, prop.i$status_OS, prop.i[[i+1]], prop.i$Type))
  names(tma.cat) <- c("OS_time","Status_OS","CD8","HR")
  med <- median(tma.cat$CD8)
  
  tma.cat$CD8[which(tma.cat$CD8 >= med)] <- "high"
  tma.cat$CD8[which(tma.cat$CD8 < med)] <- "low"
  
  #BRCA survival curve
  tma.cat.brca <- tma.cat[which(tma.cat$HR == 1),]
  formula <- paste0("fit <- survfit(Surv(OS_time,Status_OS) ~ ",colnames(tma.cat.brca)[3],
                    ", data = tma.cat.brca)")
  eval(parse(text=formula))
  pb <- ggsurvplot(fit, data=tma.cat.brca, cencor=T, pval=T, pval.coord=c(2500, 1),
                   title=paste("BRCA1/2 Proportion of", colnames(prop.i)[i+1],"OS"), risk.table=T,
                   risk.table.height = 0.3, xlim = c(0, 3000), break.time.by=500, palette=c("red","blue"))
  plotlist[[1]] <- pb
  #plotlist[[i*2-1]] <- pb
  #print(pb)

  #HR survival curve
  tma.cat.hr <- tma.cat[which(tma.cat$HR == 2),]
  formula <- paste0("fit <- survfit(Surv(OS_time,Status_OS) ~ ",colnames(tma.cat.hr)[3],
                    ", data = tma.cat.hr)")
  eval(parse(text=formula))
  ph <- ggsurvplot(fit, data=tma.cat.hr, cencor=T, pval=T, pval.coord=c(750, 1),
                   title=paste("HRwt Proportion of", colnames(prop.i)[i+1],"OS"), risk.table=T,
                   risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500, palette=c("red","blue"))
  plotlist[[2]] <- ph
  #plotlist[[i*2]] <- ph
  #print(ph)
  
  #all survival curve
  formula <- paste0("fit <- survfit(Surv(OS_time,Status_OS) ~ ",colnames(tma.cat)[3],
                    ", data = tma.cat)")
  eval(parse(text=formula))
  pa <- ggsurvplot(fit, data=tma.cat, cencor=T, pval=T, pval.coord=c(750, 0.75),
                   title=paste("All patients, Proportion of", colnames(prop.i)[i+1],"and OS"), risk.table=T,
                   risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500, palette=c("red","blue"))
  plotlist[[3]] <- pa
  
  #arrange_ggsurvplots(plotlist)
  
}
pdf("survival_prop_immune_CD8_OS_redblue.pdf")
multiplot(plotlist = plotlist)
dev.off()


# prop + tumor + PFI
prop <- data.frame("Sample"=1:44,"Epithelial"=rep(0,44),"Proliferating epithelial"=rep(0,44),
                     "EMT"=rep(0,44),"Mesenchymal"=rep(0,44),"Proliferating EMT"=rep(0,44),
                     "Functional epithelial"=rep(0,44),"Apoptotic"=rep(0,44))
prop$Sample <- unique(tumor$Sample)
names(prop)[3] <- "Proliferating epithelial"
names(prop)[6] <- "Proliferating EMT"
names(prop)[7] <- "Functional epithelial"
for (i in unique(tumor$Sample)){
  patient <- tumor[tumor$Sample==i,]
  row <- which(prop$Sample==i)
  
  for (celltype in names(prop)[2:8]){
    prop[row, celltype] <- nrow(patient[patient$Metacluster == celltype,])/nrow(patient)
  }
}
prop <- merge(prop, tma[,c(1,5,20,32)], by.x="Sample", by.y="Identifier")
distlist <- list()
plotlist <- list()
#again, this might need some editing
for (i in c(1:7)){
  
  #  p_dist <- ggplot(meex, aes_string(x=colnames(meex)[i], fill="Type")) + 
  #    geom_density(alpha=0.2)+scale_fill_brewer(palette="Dark2")
  #  distlist[[i-1]] <- p_dist
  
  tma.cat <- as.data.frame(cbind(prop$PFI_time, prop$status_PFI, prop[[i+1]], prop$Type))
  names(tma.cat) <- c("PFI_time","Status_PFI","Subtype","HR")
  med <- median(tma.cat$Subtype)
  
  tma.cat$Subtype[which(tma.cat$Subtype >= med)] <- "high"
  tma.cat$Subtype[which(tma.cat$Subtype < med)] <- "low"
  
  #BRCA survival curve
  tma.cat.brca <- tma.cat[which(tma.cat$HR == 1),]
  formula <- paste0("fit <- survfit(Surv(PFI_time,Status_PFI) ~ ",colnames(tma.cat.brca)[3],
                    ", data = tma.cat.brca)")
  eval(parse(text=formula))
  pb <- ggsurvplot(fit, data=tma.cat.brca, cencor=T, pval=T, pval.coord=c(2500, 0.75),
                   title=paste("BRCA1/2 Proportion of", colnames(prop)[i+1],"PFI"), risk.table=T,
                   risk.table.height = 0.3, xlim = c(0, 3000), break.time.by=500)
  plotlist[[i*2-1]] <- pb
  #print(pb)
  
  #HR survival curve
  tma.cat.hr <- tma.cat[which(tma.cat$HR == 2),]
  formula <- paste0("fit <- survfit(Surv(PFI_time,Status_PFI) ~ ",colnames(tma.cat.hr)[3],
                    ", data = tma.cat.hr)")
  eval(parse(text=formula))
  ph <- ggsurvplot(fit, data=tma.cat.hr, cencor=T, pval=T, pval.coord=c(750, 0.75),
                   title=paste("HRwt Proportion of", colnames(prop)[i+1],"PFI"), risk.table=T,
                   risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500)
  plotlist[[i*2]] <- ph
  #print(ph)
  
  #arrange_ggsurvplots(plotlist)
  
}
pdf("survival_prop_tumor_PFI.pdf")
multiplot(plotlist = plotlist)
dev.off()


# prop + tumor + OS
prop <- data.frame("Sample"=1:44,"Epithelial"=rep(0,44),"Proliferating epithelial"=rep(0,44),
                   "EMT"=rep(0,44),"Mesenchymal"=rep(0,44),"Proliferating EMT"=rep(0,44),
                   "Functional epithelial"=rep(0,44),"Apoptotic"=rep(0,44))
prop$Sample <- unique(tumor$Sample)
names(prop)[3] <- "Proliferating epithelial"
names(prop)[6] <- "Proliferating EMT"
names(prop)[7] <- "Functional epithelial"
for (i in unique(tumor$Sample)){
  patient <- tumor[tumor$Sample==i,]
  row <- which(prop$Sample==i)
  
  for (celltype in names(prop)[2:8]){
    prop[row, celltype] <- nrow(patient[patient$Metacluster == celltype,])/nrow(patient)
  }
}
prop <- merge(prop, tma[,c(1,5,27,33)], by.x="Sample", by.y="Identifier")
distlist <- list()
plotlist <- list()
#again, this might need some editing
for (i in c(1:7)){
  
  #  p_dist <- ggplot(meex, aes_string(x=colnames(meex)[i], fill="Type")) + 
  #    geom_density(alpha=0.2)+scale_fill_brewer(palette="Dark2")
  #  distlist[[i-1]] <- p_dist
  
  tma.cat <- as.data.frame(cbind(prop$OS_time, prop$status_OS, prop[[i+1]], prop$Type))
  names(tma.cat) <- c("OS_time","Status_OS","Subtype","HR")
  med <- median(tma.cat$Subtype)
  
  tma.cat$Subtype[which(tma.cat$Subtype >= med)] <- "high"
  tma.cat$Subtype[which(tma.cat$Subtype < med)] <- "low"
  
  #BRCA survival curve
  tma.cat.brca <- tma.cat[which(tma.cat$HR == 1),]
  formula <- paste0("fit <- survfit(Surv(OS_time,Status_OS) ~ ",colnames(tma.cat.brca)[3],
                    ", data = tma.cat.brca)")
  eval(parse(text=formula))
  pb <- ggsurvplot(fit, data=tma.cat.brca, cencor=T, pval=T, pval.coord=c(2500, 0.75),
                   title=paste("BRCA1/2 Proportion of", colnames(prop)[i+1],"OS"), risk.table=T,
                   risk.table.height = 0.3, xlim = c(0, 3000), break.time.by=500)
  plotlist[[i*2-1]] <- pb
  #print(pb)
  
  #HR survival curve
  tma.cat.hr <- tma.cat[which(tma.cat$HR == 2),]
  formula <- paste0("fit <- survfit(Surv(OS_time,Status_OS) ~ ",colnames(tma.cat.hr)[3],
                    ", data = tma.cat.hr)")
  eval(parse(text=formula))
  ph <- ggsurvplot(fit, data=tma.cat.hr, cencor=T, pval=T, pval.coord=c(750, 0.75),
                   title=paste("HRwt Proportion of", colnames(prop)[i+1],"OS"), risk.table=T,
                   risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500)
  plotlist[[i*2]] <- ph
  #print(ph)
  
  #arrange_ggsurvplots(plotlist)
  
}
pdf("survival_prop_tumor_OS.pdf")
multiplot(plotlist = plotlist)
dev.off()



### heterogeneity ----

plotlist <- list()

# immune/tumor diversity + PFI (check merge and column names)
div.surv <- merge(simp[,c(1,2)], tma[,c(1,5,20,27,32,33)], by.x="sample", by.y="Identifier")

tma.cat <- as.data.frame(cbind(div.surv$PFI_time, div.surv$status_PFI, div.surv$simpson, div.surv$Type))
names(tma.cat) <- c("PFI_time","Status_PFI","TumorHG","HR")
med <- median(tma.cat$TumorHG)

tma.cat$TumorHG[which(tma.cat$TumorHG >= med)] <- "high"
tma.cat$TumorHG[which(tma.cat$TumorHG < med)] <- "low"

#BRCA survival curve
tma.cat.brca <- tma.cat[which(tma.cat$HR == 1),]
formula <- paste0("fit <- survfit(Surv(PFI_time,Status_PFI) ~ ",colnames(tma.cat.brca)[3],
                  ", data = tma.cat.brca)")
eval(parse(text=formula))
pb <- ggsurvplot(fit, data=tma.cat.brca, cencor=T, pval=T, pval.coord=c(2500, 0.75),
                 title=paste("BRCA1/2 Tumor heterogeneity and PFI"), risk.table=T,
                 risk.table.height = 0.3, xlim = c(0, 3000), break.time.by=500, palette=c("red","blue"))
plotlist[[1]] <- pb
#print(pb)

#HR survival curve
tma.cat.hr <- tma.cat[which(tma.cat$HR == 2),]
formula <- paste0("fit <- survfit(Surv(PFI_time,Status_PFI) ~ ",colnames(tma.cat.hr)[3],
                  ", data = tma.cat.hr)")
eval(parse(text=formula))
ph <- ggsurvplot(fit, data=tma.cat.hr, cencor=T, pval=T, pval.coord=c(750, 0.75),
                 title=paste("HRwt Tumor heterogeneity and PFI"), risk.table=T,
                 risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500, palette=c("red","blue"))
plotlist[[2]] <- ph
#print(ph)

#all survival curve
formula <- paste0("fit <- survfit(Surv(PFI_time,Status_PFI) ~ ",colnames(tma.cat)[3],
                  ", data = tma.cat)")
eval(parse(text=formula))
pa <- ggsurvplot(fit, data=tma.cat, cencor=T, pval=T, pval.coord=c(750, 0.75),
                 title=paste("All patients, Tumor heterogeneity and PFI"), risk.table=T,
                 risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500, palette=c("red","blue"))
plotlist[[3]] <- pa



# immune/tumor diversity + OS (again check if changing something)
div.surv <- merge(simp[,c(1,2)], tma[,c(1,5,20,27,32,33)], by.x="sample", by.y="Identifier")

tma.cat <- as.data.frame(cbind(div.surv$OS_time, div.surv$status_OS, div.surv$simpson, div.surv$Type))
names(tma.cat) <- c("OS_time","Status_OS","TumorHG","HR")
med <- median(tma.cat$TumorHG)

tma.cat$TumorHG[which(tma.cat$TumorHG >= med)] <- "high"
tma.cat$TumorHG[which(tma.cat$TumorHG < med)] <- "low"

#BRCA survival curve
tma.cat.brca <- tma.cat[which(tma.cat$HR == 1),]
formula <- paste0("fit <- survfit(Surv(OS_time,Status_OS) ~ ",colnames(tma.cat.brca)[3],
                  ", data = tma.cat.brca)")
eval(parse(text=formula))
pb <- ggsurvplot(fit, data=tma.cat.brca, cencor=T, pval=T, pval.coord=c(2500, 0.75),
                 title=paste("BRCA1/2 Tumor heterogeneity and OS"), risk.table=T,
                 risk.table.height = 0.3, xlim = c(0, 3000), break.time.by=500, palette=c("red","blue"))
plotlist[[4]] <- pb
#print(pb)

#HR survival curve
tma.cat.hr <- tma.cat[which(tma.cat$HR == 2),]
formula <- paste0("fit <- survfit(Surv(OS_time,Status_OS) ~ ",colnames(tma.cat.hr)[3],
                  ", data = tma.cat.hr)")
eval(parse(text=formula))
ph <- ggsurvplot(fit, data=tma.cat.hr, cencor=T, pval=T, pval.coord=c(750, 0.75),
                 title=paste("HRwt Tumor heterogeneity and OS"), risk.table=T,
                 risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500, palette=c("red","blue"))
plotlist[[5]] <- ph
#print(ph)


#all survival curve
formula <- paste0("fit <- survfit(Surv(OS_time,Status_OS) ~ ",colnames(tma.cat)[3],
                  ", data = tma.cat)")
eval(parse(text=formula))
pa <- ggsurvplot(fit, data=tma.cat, cencor=T, pval=T, pval.coord=c(750, 0.75),
                 title=paste("All patients, Tumor heterogeneity and OS"), risk.table=T,
                 risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500, palette=c("red","blue"))
plotlist[[6]] <- pa


pdf("survival_hgimmune_PFI_OS_redblue.pdf")
pdf("survival_hgtumor_PFI_OS_redblue.pdf")
#arrange_ggsurvplots(plotlist)
multiplot(plotlist = plotlist)
dev.off()

#render("x.Rmd", output_format = "word_document")




### mega emt ----

# mega EMT + PFI

prop <- data.frame("Sample"=1:44,"Epithelial"=rep(0,44),"Proliferating epithelial"=rep(0,44),
                   "EMT"=rep(0,44),"Mesenchymal"=rep(0,44),"Proliferating EMT"=rep(0,44),
                   "Functional epithelial"=rep(0,44),"Apoptotic"=rep(0,44))
prop$Sample <- unique(tumor$Sample)
names(prop)[3] <- "Proliferating epithelial"
names(prop)[4] <- "Mega-EMT"
names(prop)[6] <- "Proliferating EMT"
names(prop)[7] <- "Functional epithelial"
for (i in unique(tumor$Sample)){
  patient <- tumor[tumor$Sample==i,]
  row <- which(prop$Sample==i)
  
  for (celltype in names(prop)[2:8]){
    prop[row, celltype] <- nrow(patient[patient$Metacluster == celltype,])/nrow(patient)
  }
}
for (i in unique(prop$Sample)){
  prop[which(prop$Sample==i),4] <- prop[which(prop$Sample==i),4]+prop[which(prop$Sample==i),6]
}

megaemt <- merge(prop[,c(1,4)], tma[,c(1,5,20,27,32,33)], by.x="Sample", by.y="Identifier")

tma.cat <- as.data.frame(cbind(megaemt$PFI_time, megaemt$status_PFI, megaemt$`Mega-EMT`, megaemt$Type))
names(tma.cat) <- c("PFI_time","Status_PFI","MegaEMT","HR")
med <- median(tma.cat$MegaEMT)

tma.cat$MegaEMT[which(tma.cat$MegaEMT >= med)] <- "high"
tma.cat$MegaEMT[which(tma.cat$MegaEMT < med)] <- "low"

#BRCA survival curve
tma.cat.brca <- tma.cat[which(tma.cat$HR == 1),]
formula <- paste0("fit <- survfit(Surv(PFI_time,Status_PFI) ~ ",colnames(tma.cat.brca)[3],
                  ", data = tma.cat.brca)")
eval(parse(text=formula))
pb <- ggsurvplot(fit, data=tma.cat.brca, cencor=T, pval=T, pval.coord=c(2500, 0.75),
                 title=paste("BRCA1/2 Proportion of mega-EMT and PFI"), risk.table=T,
                 risk.table.height = 0.3, xlim = c(0, 3000), break.time.by=500, palette=c("red","blue"))
plotlist[[1]] <- pb
#print(pb)

#HR survival curve
tma.cat.hr <- tma.cat[which(tma.cat$HR == 2),]
formula <- paste0("fit <- survfit(Surv(PFI_time,Status_PFI) ~ ",colnames(tma.cat.hr)[3],
                  ", data = tma.cat.hr)")
eval(parse(text=formula))
ph <- ggsurvplot(fit, data=tma.cat.hr, cencor=T, pval=T, pval.coord=c(750, 0.75),
                 title=paste("HRwt Proportion of mega-EMT and PFI"), risk.table=T,
                 risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500, palette=c("red","blue"))
plotlist[[2]] <- ph
#print(ph)

#all survival curve
formula <- paste0("fit <- survfit(Surv(PFI_time,Status_PFI) ~ ",colnames(tma.cat)[3],
                  ", data = tma.cat)")
eval(parse(text=formula))
pa <- ggsurvplot(fit, data=tma.cat, cencor=T, pval=T, pval.coord=c(750, 0.75),
                 title=paste("All patients, Proportion of mega-EMT and PFI"), risk.table=T,
                 risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500, palette=c("red","blue"))
plotlist[[3]] <- pa



# immune diversity + OS
tma.cat <- as.data.frame(cbind(megaemt$OS_time, megaemt$status_OS, megaemt$`Mega-EMT`, megaemt$Type))
names(tma.cat) <- c("OS_time","Status_OS","MegaEMT","HR")
med <- median(tma.cat$MegaEMT)

tma.cat$MegaEMT[which(tma.cat$MegaEMT >= med)] <- "high"
tma.cat$MegaEMT[which(tma.cat$MegaEMT < med)] <- "low"

#BRCA survival curve
tma.cat.brca <- tma.cat[which(tma.cat$HR == 1),]
formula <- paste0("fit <- survfit(Surv(OS_time,Status_OS) ~ ",colnames(tma.cat.brca)[3],
                  ", data = tma.cat.brca)")
eval(parse(text=formula))
pb <- ggsurvplot(fit, data=tma.cat.brca, cencor=T, pval=T, pval.coord=c(2500, 0.75),
                 title=paste("BRCA1/2 Proportion of mega-EMT and OS"), risk.table=T,
                 risk.table.height = 0.3, xlim = c(0, 3000), break.time.by=500, palette=c("red","blue"))
plotlist[[4]] <- pb
#print(pb)

#HR survival curve
tma.cat.hr <- tma.cat[which(tma.cat$HR == 2),]
formula <- paste0("fit <- survfit(Surv(OS_time,Status_OS) ~ ",colnames(tma.cat.hr)[3],
                  ", data = tma.cat.hr)")
eval(parse(text=formula))
ph <- ggsurvplot(fit, data=tma.cat.hr, cencor=T, pval=T, pval.coord=c(750, 0.75),
                 title=paste("HRwt Proportion of mega-EMT and OS"), risk.table=T,
                 risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500, palette=c("red","blue"))
plotlist[[5]] <- ph
#print(ph)


#all survival curve
formula <- paste0("fit <- survfit(Surv(OS_time,Status_OS) ~ ",colnames(tma.cat)[3],
                  ", data = tma.cat)")
eval(parse(text=formula))
pa <- ggsurvplot(fit, data=tma.cat, cencor=T, pval=T, pval.coord=c(750, 0.75),
                 title=paste("All patients, Proportion of mega-EMT and OS"), risk.table=T,
                 risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500, palette=c("red","blue"))
plotlist[[6]] <- pa


pdf("survival_megaemt_PFI_OS_redblue.pdf")
#arrange_ggsurvplots(plotlist)
multiplot(plotlist = plotlist)
dev.off()



### only brca vs hr ----

plotlist <- list()

#OS
tma.cat <- as.data.frame(cbind(megaemt$OS_time, megaemt$status_OS, megaemt$Type))
names(tma.cat) <- c("OS_time","Status_OS","HR")
formula <- paste0("fit <- survfit(Surv(OS_time,Status_OS) ~ ",colnames(tma.cat)[3],
                  ", data = tma.cat)")
eval(parse(text=formula))
po <- ggsurvplot(fit, data=tma.cat, cencor=T, pval=T, pval.coord=c(750, 0.65),
                 title=paste("All patients, BRCA status and OS"), risk.table=T,
                 risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500,
                 palette=c("blue","red"), legend.labs=c("BRCA1/2mut","HRwt"))
plotlist[[1]] <- po

#PFI
tma.cat <- as.data.frame(cbind(megaemt$PFI_time, megaemt$status_PFI, megaemt$Type))
names(tma.cat) <- c("PFI_time","Status_PFI","HR")
formula <- paste0("fit <- survfit(Surv(PFI_time,Status_PFI) ~ ",colnames(tma.cat)[3],
                  ", data = tma.cat)")
eval(parse(text=formula))
pp <- ggsurvplot(fit, data=tma.cat, cencor=T, pval=T, pval.coord=c(750, 0.75),
                 title=paste("All patients, BRCA status and PFI"), risk.table=T,
                 risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500,
                 palette=c("blue","red"), legend.labs=c("BRCA1/2mut","HRwt"))
plotlist[[2]] <- pp


setwd("C:/LocalData/nuplyyt/cyto/celltypes/images/maybe_paper/Survival")
pdf("survival_BRCAstatus_PFI_OS_redblue_strata.pdf")
#arrange_ggsurvplots(plotlist)
multiplot(plotlist = plotlist)
dev.off()








asd
### old stuff -------

CD4_surv <- all_immune[which(all_immune$Subtype == "CD4"),]
#median expressions in each patient
meex <- CD4_surv %>% group_by(Sample) %>% summarise(median_Ki67 = median(Ki67),
                                                    median_pSTAT1 = median(pSTAT1),
                                                    median_P21 = median(P21),
                                                    median_PD1=median(PD1),
                                                    median_PDL1=median(PDL1),
                                                    median_cCasp3=median(cCasp3))

meex <- merge(meex, tma[,c(1,5,20,32)], by.x="Sample", by.y="Identifier")

for (i in seq_along(colnames(meex))[2:7]){
  #p_hist <- hist(median_expressions[, i], breaks = 50)
  p_dist <- ggplot(meex, aes_string(x=colnames(meex)[i], fill="Type")) + 
    geom_density(alpha=0.2)+scale_fill_brewer(palette="Dark2")
  print(p_dist)
  
  #seuraava jakaa datan kategorioihin (all_data_TMA_cat)
  #all_data_TMA_cut <- surv_cutpoint(median_expressions, time="PFI_time", event="status_PFI",
                                    #variables=colnames(median_expressions)[i])
  #p_density <- plot(all_data_TMA_cut, colnames(median_expressions)[i], palette="npg")
  #print(p_density)
  
  tma.cat <- as.data.frame(cbind(meex$PFI_time, meex$status_PFI, meex[[i]], meex$Type))
  names(tma.cat) <- c("PFI_time","Status_PFI","Marker","HR")
  med <- median(tma.cat$Marker)
  
  tma.cat$Marker[which(tma.cat$Marker > med)] <- "high"
  tma.cat$Marker[which(tma.cat$Marker < med)] <- "low"
  
  #to categorize
  #all_data_TMA_cat <- surv_categorize(all_data_TMA_cut)
  #all_data_TMA_cat <- cbind(all_data_TMA_cat, median_expressions$Type)
  #colnames(all_data_TMA_cat)[4] <- "HR" 
  
  #BRCA survival curve
  tma.cat.brca <- tma.cat[which(tma.cat$HR == 1),]
  formula <- paste0("fit <- survfit(Surv(PFI_time,Status_PFI) ~ ",colnames(tma.cat.brca)[3],
                    ", data = tma.cat.brca)")
  eval(parse(text=formula))
  pb <- ggsurvplot(fit, data=tma.cat.brca, cencor=T, pval=T, pval.coord=c(2500, 0.75),
                  title=paste("BRCA1/2 CD4", colnames(meex)[i],"PFI"), risk.table=T,
                  risk.table.height = 0.3, xlim = c(0, 3000), break.time.by=500)
  print(pb)
  
  #HR survival curve
  tma.cat.hr <- tma.cat[which(tma.cat$HR == 2),]
  formula <- paste0("fit <- survfit(Surv(PFI_time,Status_PFI) ~ ",colnames(tma.cat.hr)[3],
                    ", data = tma.cat.hr)")
  eval(parse(text=formula))
  ph <- ggsurvplot(fit, data=tma.cat.hr, cencor=T, pval=T, pval.coord=c(750, 0.75),
                  title=paste("HRwt CD4", colnames(meex)[i],"PFI"), risk.table=T,
                  risk.table.height = 0.3, xlim = c(0, 1000), break.time.by=500)
  print(ph)
  
  arrange_ggsurvplots(plots)
}

res <- arrange_ggsurvplots(plots, print = FALSE)
ggsave("myfile.pdf", res)

multiplot(plotlist=plotlist,cols = 2)