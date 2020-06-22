### creation of simpson indexes ----

## tumor

#install.packages("vegan")
#install.packages("sm")
#install.packages("ggpubr")
#install.packages("Hmisc")
library(Hmisc)
library(vegan)
library(sm)
library(ISLR)
library(ggpubr)
library(grid)
library(gridExtra)

metaclusters <- read.csv("C:/LocalData/nuplyyt/cyto/celltypes/round_april/tumor_for_plots_metacluster_noCK7_kmeans_MST_splitEMT.csv")

counts <- data.frame("Sample"=1:44,"Epithelial"=rep(0,44),"Proliferating epithelial"=rep(0,44),
                     "EMT"=rep(0,44),"Mesenchymal"=rep(0,44),"Proliferating EMT"=rep(0,44),
                     "Hyperfunctional epithelial"=rep(0,44),"Apoptotic"=rep(0,44),
                     "Negative"=rep(0,44),"Total"=rep(0,44))

counts$Sample <- unique(metaclusters$Sample)
names(counts)[3] <- "Proliferating epithelial"
names(counts)[6] <- "Proliferating EMT"
names(counts)[7] <- "Functional epithelial"

for (i in unique(metaclusters$Sample)){
  patient <- metaclusters[metaclusters$Sample==i,]
  row <- which(counts$Sample==i)
  
  for (celltype in names(counts)[2:9]){
    counts[row, celltype] <- nrow(patient[patient$Metacluster == celltype,])
  }
  
  counts[row,10] <- nrow(patient)
}

simpson <- diversity(counts[,1:9], MARGIN = 1)
simps <- data.frame(sample=unique(metaclusters$Sample), simpson)


## immune

immuneclusters <- read.csv("C:/LocalData/nuplyyt/cyto/celltypes/round_april/tsi_for_plots_immune_new.csv")
immuneclusters$Subtype <- as.character(immuneclusters$Subtype)
immuneclusters$Subtype[which(immuneclusters$Subtype == "pDC")] <- "APC"
immuneclusters$Subtype[which(immuneclusters$Subtype == "mDC")] <- "APC"


counts.i <- data.frame("Sample"=1:44,"CD8"=rep(0,44),"M2.macrop"=rep(0,44),
                     "B.cells"=rep(0,44),"APC"=rep(0,44),"M1.macrop"=rep(0,44),
                     "T.regs"=rep(0,44),"Neutrophils"=rep(0,44),"CD4"=rep(0,44),"Total"=rep(0,44))


counts.i$Sample <- unique(immuneclusters$Sample)


for (i in unique(immuneclusters$Sample)){
  patient <- immuneclusters[immuneclusters$Sample==i,]
  row <- which(counts.i$Sample==i)
  
  for (celltype in names(counts.i)[2:9]){
    counts.i[row, celltype] <- nrow(patient[patient$Subtype == celltype,])
  }
  
  counts.i[row,10] <- nrow(patient)
}

simpson.i <- diversity(counts.i, MARGIN = 1)
simps.i <- data.frame(sample=unique(immuneclusters$Sample), simpson.i)


compositions <- read.csv("C:/LocalData/nuplyyt/heterogeneity/immune_annotations_m2_split.csv")
compositions$Sample <- sub("id","",compositions$Sample)
# you need to order Sample column if you merge this with others!


## merging simpson and clinical

clinical <- read.csv("C:/LocalData/nuplyyt/all_data_TMA_annotation.csv")
clinical$Sample <- sub("id","",clinical$Sample)
clinical$PFI[which(clinical$PFI=="medium")] <- "short"
clinical$PFI <- factor(clinical$PFI)
simp <- cbind(simps, simps.i[,2], clinical[,2:3])
names(simp)[3] <- "simpson.i"

brca  <- simp[simp$Type == "BRCA1/2 mutated",]
brca$Type <- rep("BRCAmut",31)
wt    <- simp[simp$Type == "HRwt",]
data  <- rbind(brca,wt)
data$Type <- as.factor(data$Type)
short <- data[data$PFI == "short",]
long  <- data[data$PFI == "long",]

### density plots ----

makeplot <- function(x, groups, title, group1, group2){
  sm.density.compare(x, groups, xlab="Simpson diversity index")
  title(main = title)
  colfill<-c(2:(2+length(levels(groups))))
  #legend("topright", inset=c(0,0), levels(groups), fill=colfill, cex=0.9)
  
  w <- wilcox.test(group1, group2)
  s <- signif(w$p.value,2)
  p <- paste("p=",as.character(s))
  mtext(p, 3, adj=1, cex=0.8)
}


#tumor hg in comparison to HR status and PFI - density plots

par(mfrow=c(2,3))

makeplot(data$simpson, data$Type, title="Heterogeneity in BRCAmut vs HRwt",
         group1 = brca$simpson, group2 = wt$simpson)
makeplot(brca$simpson, brca$PFI, title="Heterogeneity in BRCAmut",
         group1 = brca$simpson[brca$PFI=="short"], group2 = brca$simpson[brca$PFI=="long"])
makeplot(wt$simpson, wt$PFI, title="Heterogeneity in HRwt",
         group1 = wt$simpson[wt$PFI=="short"], group2 = wt$simpson[wt$PFI=="long"])
makeplot(data$simpson, data$PFI, title="Heterogeneity in short vs long PFI",
         group1 = short$simpson, group2 = long$simpson)
makeplot(short$simpson, short$Type, title="Heterogeneity in short PFI",
         group1 = short$simpson[short$Type=="BRCAmut"], group2 = short$simpson[short$Type=="HRwt"])
makeplot(long$simpson, long$Type, title="Heterogeneity in long PFI",
         group1 = long$simpson[long$Type=="BRCAmut"], group2 = long$simpson[long$Type=="HRwt"])


#same but with immune hg - density plots

par(mfrow=c(2,3))

makeplot(data$simpson.i, data$Type, title="Immune heterogeneity in BRCAmut vs HRwt",
         group1 = brca$simpson.i, group2 = wt$simpson.i)
makeplot(short$simpson.i, short$Type, title="Immune heterogeneity in short PFI",
         group1 = short$simpson.i[short$Type=="BRCAmut"], group2 = short$simpson.i[short$Type=="HRwt"])
makeplot(long$simpson.i, long$Type, title="Immune heterogeneity in long PFI",
         group1 = long$simpson.i[long$Type=="BRCAmut"], group2 = long$simpson.i[long$Type=="HRwt"])
makeplot(data$simpson.i, data$PFI, title="Immune heterogeneity in short vs long PFI",
         group1 = short$simpson.i, group2 = long$simpson.i)
makeplot(brca$simpson.i, brca$PFI, title="Immune heterogeneity in BRCAmut",
         group1 = brca$simpson.i[brca$PFI=="short"], group2 = brca$simpson.i[brca$PFI=="long"])
makeplot(wt$simpson.i, wt$PFI, title="Immune heterogeneity in HRwt",
         group1 = wt$simpson.i[wt$PFI=="short"], group2 = wt$simpson.i[wt$PFI=="long"])


### box plots - HR and PFI ----

#tumor hg in comparison to HR status and PFI - box plots

boxes <- function(a,b,names,main,level){
  boxplot(a, b, names = names, ylab="Simpson's diversity index",main=main,outline=FALSE)
  
  levelProportions <- summary(level)/(length(a)+length(b))

  myjitter <- jitter(rep(1, length(a)), amount=levelProportions[1]/2)
  points(myjitter, a, pch=20, col="blue") #"#003366"
    
  myjitter2 <- jitter(rep(2, length(b)), amount=levelProportions[2]/2)
  points(myjitter2, b, pch=20, col="red") #"#3399CC" 
  
  w <- wilcox.test(a, b)
  s <- signif(w$p.value,2)
  p <- paste("p=",as.character(s))
  mtext(p, 3, adj=1, cex=0.8)
  #if(s<0.05){
  #  mtext(p, 3, adj=1, cex=0.8, col="red")
  #} else{
  #  mtext(p, 3, adj=1, cex=0.8)
  #} 
}

par(mfrow=c(2,3))
par(mfrow=c(2,2))

boxes(brca$simpson, wt$simpson, names=c("BRCAmut","HRwt"), main="Tumor heterogeneity in BRCAmut vs HRwt",level=data$Type)
boxes(brca$simpson[brca$PFI=="short"], wt$simpson[wt$PFI=="short"], names = c("BRCAmut","HRwt"),
      main="Tumor heterogeneity in short PFI",level=data$Type[data$PFI=="short"])
boxes(brca$simpson[brca$PFI=="long"], wt$simpson[wt$PFI=="long"], names = c("BRCAmut","HRwt"),
      main="Tumor heterogeneity in long PFI",level=data$Type[data$PFI=="long"])
boxes(short$simpson, long$simpson, names=c("Short","Long"), main="Tumor heterogeneity in short vs long PFI",level=data$PFI)
boxes(short$simpson[short$Type=="BRCAmut"], long$simpson[long$Type=="BRCAmut"], names = c("short","long"),
      main="Tumor heterogeneity in BRCAmut",level=data$PFI[data$Type=="BRCAmut"])
boxes(short$simpson[short$Type=="HRwt"], long$simpson[long$Type=="HRwt"], names = c("short","long"),
      main="Tumor heterogeneity in HRwt",level=data$PFI[data$Type=="HRwt"])

#immune hg in comparison to HR status and PFI - box plots

par(mfrow=c(2,3))

boxes(brca$simpson.i, wt$simpson.i, names=c("BRCAmut","HRwt"), main="Immune heterogeneity in BRCAmut vs HRwt",level=data$Type)
boxes(brca$simpson.i[brca$PFI=="short"], wt$simpson.i[wt$PFI=="short"], names = c("BRCAmut","HRwt"),
      main="Immune heterogeneity in short PFI",level=data$Type[data$PFI=="short"])
boxes(brca$simpson.i[brca$PFI=="long"], wt$simpson.i[wt$PFI=="long"], names = c("BRCAmut","HRwt"),
      main="Immune heterogeneity in long PFI",level=data$Type[data$PFI=="long"])
boxes(short$simpson.i, long$simpson.i, names=c("Short","Long"), main="Immune heterogeneity in short vs long PFI",level=data$PFI)
boxes(short$simpson.i[short$Type=="BRCAmut"], long$simpson.i[long$Type=="BRCAmut"], names = c("short","long"),
      main="Immune heterogeneity in BRCAmut",level=data$PFI[data$Type=="BRCAmut"])
boxes(short$simpson.i[short$Type=="HRwt"], long$simpson.i[long$Type=="HRwt"], names = c("short","long"),
      main="Immune heterogeneity in HRwt",level=data$PFI[data$Type=="HRwt"])


### correlation between immune hg and tumor hg ----

par(mfrow=c(1,1))

plot(data$simpson, data$simpson.i, main="Heterogeneity",
     xlab="Tumor heterogeneity", ylab="Immune heterogeneity",
     col=c("red","orange")[data$Type], pch=19)
legend(x="topright", legend = levels(data$Type), col=c("red","orange"), pch=19)
abline(lm(data$simpson.i~data$simpson), col="blue") # regression line (y~x)

c1 <- cor.test(data$simpson, data$simpson.i, method=c("spearman"))
c1s <- signif(c1$p.value,2)
c1r <- signif(c1$estimate,2)
c1p <- paste("p=",as.character(c1s))
c1rho <- paste("S: rho=",as.character(c1r))
mtext(c1p, 3, adj=1, cex=0.8)
mtext(c1rho, 3, adj=0, cex=0.8)

#different version with ggplot
library(RColorBrewer)

logdata <- data
logdata[,2:3] <- log(logdata[,2:3])

#one line
reg<-lm(logdata$simpson.i ~ logdata$simpson)
co=coefficients(reg)

c1 <- cor.test(logdata$simpson, logdata$simpson.i, method=c("spearman"))
c1s <- signif(c1$p.value,2)
c1r <- signif(c1$estimate,2)
c <- paste0("p=",as.character(c1s),", R=",as.character(c1r))

mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(3)

ggplot(logdata, aes(x=simpson, y=simpson.i, color=Type))+
  geom_point(size=3)+
  scale_colour_manual(values=c("blue","red"), name="HR status")+
  geom_abline(intercept = co[1], slope = co[2], color="#666666")+
  labs(subtitle = c, x="Log of Simpson's diversity index", y="Log of Simpson's diversity index", tag="C")+
  theme(plot.subtitle=element_text(size=10), plot.title=element_text(size=13),
        axis.title=element_text(size=9), axis.text=element_text(size=8),
        legend.position=c(0.2,0.85))+
  coord_fixed()


#split line
logdata <- logdata[order(logdata$simpson),]

reg_brca<-lm(logdata$simpson.i[logdata$Type=="BRCAmut"] ~ logdata$simpson[logdata$Type=="BRCAmut"])
reg_wt<-lm(logdata$simpson.i[logdata$Type=="HRwt"] ~ logdata$simpson[logdata$Type=="HRwt"])
reg_mut<-lm(simpson.i ~ simpson, data = logdata)
summary(reg_mut)

cob=coefficients(reg_brca)
cow=coefficients(reg_wt)

c1 <- cor.test(logdata$simpson[logdata$Type=="BRCAmut"], logdata$simpson.i[logdata$Type=="BRCAmut"], method=c("spearman"))
c1s <- signif(c1$p.value,2)
c1r <- signif(c1$estimate,2)
c <- paste0("BRCAmut: p=",as.character(c1s),", R=",as.character(c1r))

c1.w <- cor.test(logdata$simpson[logdata$Type=="HRwt"], logdata$simpson.i[logdata$Type=="HRwt"], method=c("spearman"))
c1s.w <- signif(c1.w$p.value,2)
c1r.w <- signif(c1.w$estimate,2)
cw <- paste0(c,"  /  HRwt: p=",as.character(c1s.w),", R=",as.character(c1r.w))

mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(3)

ggplot(logdata, aes(x=simpson, y=simpson.i, color=Type))+
  geom_point(size=3)+
  scale_colour_manual(values=mycolors, name="HR status")+
  geom_abline(intercept = cob[1], slope = cob[2], color="#1B9E77")+
  #geom_abline(intercept = cow[1], slope = cow[2], color="#A66753")+
  #labs(subtitle = cw, x="Log of tumor heterogeneity", y="Log of immune heterogeneity", tag="A")+
  theme(plot.subtitle=element_text(size=10), plot.title=element_text(size=13),
        axis.title=element_text(size=9), axis.text=element_text(size=8),
        legend.position=c(0.2,0.85))+
        coord_fixed()

### multiple scatterplots - composition ----
#correlation between tumor hg and immune composition

comp <- compositions[order(compositions$Sample),]
data <- data[order(data$sample),]
comp <- cbind(comp, data[,2:5])

comp$immune_cluster <- as.character(comp$immune_cluster)
comp$immune_cluster[which(comp$immune_cluster == "APC rich")] <- "APC"
comp$immune_cluster[which(comp$immune_cluster == "APC+M1 macrophage rich")] <- "APC+M1"
comp$immune_cluster[which(comp$immune_cluster == "B-cell rich")] <- "Bcell"
comp$immune_cluster[which(comp$immune_cluster == "CD4+ effector T-cell rich")] <- "CD4"
comp$immune_cluster[which(comp$immune_cluster == "CD8+ T-cell rich")] <- "CD8"
comp$immune_cluster[which(comp$immune_cluster == "M2 macrophage + CD8+ T-cell rich")] <- "M2+CD8"
comp$immune_cluster[which(comp$immune_cluster == "M2 macrophage rich")] <- "M2"
comp$immune_cluster[which(comp$immune_cluster == "myeloid cell rich")] <- "myeloid"
comp$immune_cluster[which(comp$immune_cluster == "neutrophil rich")] <- "NP"
comp$immune_cluster[which(comp$immune_cluster == "Treg rich")] <- "Treg"

# this needs to be done separately for each comparison
my_comparisons <- list(c("M2","Bcell"))
ggboxplot(comp[,c(2:3)],
          x = "immune_cluster",
          y = "simpson",
          color = "immune_cluster",
          add = "jitter",
          outlier.shape = NA)+
  stat_compare_means(method = "kruskal.test")+
  stat_compare_means(comparisons = my_comparisons, label.y = 0.22)+
  theme(legend.position="none")+
  ggtitle("Heterogeneity and immune composition")



#correlation between immune hg and tumor composition

t.comp <- cbind(data, counts[,2:10])

par(mfrow=c(1,1))
require(stats)

#NOTE! this same function is used in many places below, so you need to check carefully
#1) whether you want to use i.comp or t.comp
#2) "immune" or "tumor" in the plot subtitles and axes
#3) the text in "d" if you want to recreate these plots
plot2 <- function(a,b,tumor){
  reg<-lm(b ~ a, data = i.comp)
  co=coefficients(reg)
  
  c1 <- cor.test(a, b, method=c("spearman"))
  c1s <- signif(c1$p.value,2)
  c1r <- signif(c1$estimate,2)
  
  c <- paste0("p=",as.character(c1s),", R=",as.character(c1r))
  d <- paste0("Proportion of ", tumor)
  
  ggplot(i.comp, aes(x=a, y=b))+
    geom_point()+
    geom_abline(intercept = co[1], slope = co[2], color="blue")+
    labs(subtitle = c, title=paste0("Tumor heterogeneity and ",tumor),
         x="Simpson diversity index for tumor cells", y=d)+
    theme(plot.subtitle=element_text(size=5), plot.title=element_text(size=7),
          axis.title=element_text(size=5), axis.text=element_text(size=4))
}

par(mfrow=c(2,1))

p1 <- plot2(t.comp$simpson.i, (t.comp$Epithelial/t.comp$Total),tumor="epithelial")
p2 <- plot2(t.comp$simpson.i, (t.comp$`Proliferating epithelial`/t.comp$Total),tumor="proliferating epithelial")
p3 <- plot2(t.comp$simpson.i, (t.comp$EMT/t.comp$Total),tumor="EMT")
p4 <- plot2(t.comp$simpson.i, (t.comp$Mesenchymal/t.comp$Total),tumor="mesenchymal")
p5 <- plot2(t.comp$simpson.i, (t.comp$`Proliferating EMT`/t.comp$Total),tumor="proliferating EMT")
p6 <- plot2(t.comp$simpson.i, (t.comp$`Hyperfunctional epithelial`/t.comp$Total),tumor="hyperfunctional epithelial")
p7 <- plot2(t.comp$simpson.i, (t.comp$Apoptotic/t.comp$Total),tumor="apoptotic")
grid.arrange(p1, p2, p3, p4, p5, p6, p7, nrow = 3)

# for (i in names(t.comp)[6:13]){
#   plot2(t.comp$simpson.i, (t.comp[[i]]/t.comp$Total),tumor=i)
# }


#correlation between tumor hg and immune composition

i.comp <- cbind(data, counts.i[,2:10])

#remember to change plot2 function
p1 <- plot2(i.comp$simpson, (i.comp$CD8/i.comp$Total),tumor="CD8 T cells")
p2 <- plot2(i.comp$simpson, (i.comp$M2.macrop/i.comp$Total),tumor="M2 macrophages")
p3 <- plot2(i.comp$simpson, (i.comp$B.cells/i.comp$Total),tumor="B cells")
p4 <- plot2(i.comp$simpson, (i.comp$APC/i.comp$Total),tumor="APC")
p5 <- plot2(i.comp$simpson, (i.comp$M1.macrop/i.comp$Total),tumor="M1 macrophages")
p6 <- plot2(i.comp$simpson, (i.comp$T.regs/i.comp$Total),tumor="Tregs")
p7 <- plot2(i.comp$simpson, (i.comp$Neutrophils/i.comp$Total),tumor="neutrophils")
p8 <- plot2(i.comp$simpson, (i.comp$CD4/i.comp$Total),tumor="CD4 effector T cells")
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 3)


#different version for splitting the line
#NOTE! check for immune/tumor and t/i.comp before use!
data <- data[order(data$sample),]
i.comp <- i.comp[order(i.comp$sample),]
logdata <- cbind(data,i.comp$CD8, i.comp$Total)
names(logdata)[6:7] <- c("CD8", "Total")
#logdata[,2:3] <- log(logdata[,2:3])
#logdata <- logdata[order(logdata$simpson),]

reg_brca<-lm((logdata$CD8[logdata$Type=="BRCAmut"]/logdata$Total[logdata$Type=="BRCAmut"]) ~ logdata$simpson[logdata$Type=="BRCAmut"])
reg_wt<-lm((logdata$CD8[logdata$Type=="HRwt"]/logdata$Total[logdata$Type=="HRwt"]) ~ logdata$simpson[logdata$Type=="HRwt"])
reg_all<-lm((logdata$CD8/logdata$Total) ~ logdata$simpson)

cob=coefficients(reg_brca)
cow=coefficients(reg_wt)
coa=coefficients(reg_all)

c1 <- cor.test(logdata$simpson[logdata$Type=="BRCAmut"],
               (logdata$CD8[logdata$Type=="BRCAmut"]/logdata$Total[logdata$Type=="BRCAmut"]),
               method=c("spearman"))
c1s <- signif(c1$p.value,2)
c1r <- signif(c1$estimate,2)
c <- paste0("BRCAmut: p=",as.character(c1s),", R=",as.character(c1r))

c1.w <- cor.test(logdata$simpson[logdata$Type=="HRwt"],
                 (logdata$CD8[logdata$Type=="HRwt"]/logdata$Total[logdata$Type=="HRwt"]),
                 method=c("spearman"))
c1s.w <- signif(c1.w$p.value,2)
c1r.w <- signif(c1.w$estimate,2)
cw <- paste0(c,"  /  HRwt: p=",as.character(c1s.w),", R=",as.character(c1r.w))

c1.a <- cor.test(logdata$simpson,(logdata$CD8/logdata$Total),method=c("spearman"))
c1s.a <- signif(c1.a$p.value,2)
c1r.a <- signif(c1.a$estimate,2)
#ca <- paste0(cw,"  /  All: p=",as.character(c1s.a),", R=",as.character(c1r.a))
ca <- paste0("p=", as.character(c1s.a),", R=", as.character(c1r.a))

mycolors <- c("blue","red")
#colorRampPalette(brewer.pal(8, "Dark2"))(3)

ggplot(logdata, aes(x=simpson, y=(CD8/Total), color=Type))+
  geom_point(size=3)+
  scale_colour_manual(values=mycolors, name="HR status")+
  geom_abline(intercept = cob[1], slope = cob[2], color="blue")+ ##1B9E77
  geom_abline(intercept = cow[1], slope = cow[2], color="red")+ ##A66753
  #geom_abline(intercept = coa[1], slope = coa[2], color="black")+
  labs(subtitle = cw, x="Simpson's diversity index in tumor", y="Proportion of CD8 cells")+
  theme(plot.subtitle=element_text(size=10), plot.title=element_text(size=13),
        axis.title=element_text(size=9), axis.text=element_text(size=8),
        legend.position=c(0.2,0.85))+
  coord_fixed()


### multiple scatterplots - markers ----

#correlation between hg's and marker activity

#tumor marker activity

medex <- metaclusters[1:44,c(4,12,20:22,24:25,28,32)]
medex$Sample <- unique(metaclusters$Sample)
medex <- medex[order(medex$Sample),]
t.comp <- t.comp[order(t.comp$sample),]

for (i in medex$Sample){
  patient <- metaclusters[metaclusters$Sample==i,]
  row <- which(medex$Sample==i)
  
  for (m in names(medex)[2:9]){
    med <- median(patient[[m]])
    medex[row,m] <- med
  }
}

i.comp <- i.comp[order(i.comp$sample),]
medex <- cbind(medex, i.comp[,2:3])

#remember to check immune and tumor when using
plot3 <- function(a,b,marker){
  reg<-lm(b ~ a, data = imedex)
  co=coefficients(reg)
  
  c1 <- cor.test(a, b, method=c("spearman"))
  c1s <- signif(c1$p.value,2)
  c1r <- signif(c1$estimate,2)
  
  c <- paste0("p=",as.character(c1s),", R=",as.character(c1r))
  d <- paste0(marker," median expression in immune cells")
  
  ggplot(imedex, aes(x=a, y=b))+
    geom_point()+
    geom_abline(intercept = co[1], slope = co[2], color="blue")+
    labs(subtitle = c, title=paste0("Immune heterogeneity and ",marker, " in immune"),
         x="Simpson heterogeneity index in immune cells", y=d)+
    theme(plot.subtitle=element_text(size=5), plot.title=element_text(size=7),
          axis.title=element_text(size=5), axis.text=element_text(size=4))
}

#immune hg
p1 <- plot3(medex$simpson.i, medex$BP53_1,marker="BP53")
p2 <- plot3(medex$simpson.i, medex$cCasp3,marker="cCasp3")
p3 <- plot3(medex$simpson.i, medex$pSTAT1,marker="pSTAT1")
p4 <- plot3(medex$simpson.i, medex$yH2AX,marker="yH2AX")
p5 <- plot3(medex$simpson.i, medex$Ki67,marker="Ki67")
p6 <- plot3(medex$simpson.i, medex$PDL1,marker="PDL1")
p7 <- plot3(medex$simpson.i, medex$P21,marker="P21")
p8 <- plot3(medex$simpson.i, medex$PD1,marker="PD1")
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 3)


#tumor hg
p1 <- plot3(medex$simpson, medex$BP53_1,marker="BP53")
p2 <- plot3(medex$simpson, medex$cCasp3,marker="cCasp3")
p3 <- plot3(medex$simpson, medex$pSTAT1,marker="pSTAT1")
p4 <- plot3(medex$simpson, medex$yH2AX,marker="yH2AX")
p5 <- plot3(medex$simpson, medex$Ki67,marker="Ki67")
p6 <- plot3(medex$simpson, medex$PDL1,marker="PDL1")
p7 <- plot3(medex$simpson, medex$P21,marker="P21")
p8 <- plot3(medex$simpson, medex$PD1,marker="PD1")
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 3)


#immune marker activity

imedex <- immuneclusters[1:44,c(3,11,19:21,23:24,27,31)]
imedex$Sample <- unique(immuneclusters$Sample)
imedex <- imedex[order(imedex$Sample),]

for (i in imedex$Sample){
  patient <- immuneclusters[immuneclusters$Sample==i,]
  row <- which(imedex$Sample==i)
  
  for (m in names(imedex)[2:9]){
    med <- median(patient[[m]])
    imedex[row,m] <- med
  }
}

i.comp <- i.comp[order(i.comp$sample),]
imedex <- cbind(imedex, i.comp[,2:3])

#tumor hg
p1 <- plot3(imedex$simpson, imedex$BP53_1,marker="BP53")
p2 <- plot3(imedex$simpson, imedex$cCasp3,marker="cCasp3")
p3 <- plot3(imedex$simpson, imedex$pSTAT1,marker="pSTAT1")
p4 <- plot3(imedex$simpson, imedex$yH2AX,marker="yH2AX")
p5 <- plot3(imedex$simpson, imedex$Ki67,marker="Ki67")
p6 <- plot3(imedex$simpson, imedex$PDL1,marker="PDL1")
p7 <- plot3(imedex$simpson, imedex$P21,marker="P21")
p8 <- plot3(imedex$simpson, imedex$PD1,marker="PD1")
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 3)

#immune hg
p1 <- plot3(imedex$simpson.i, imedex$BP53_1,marker="BP53")
p2 <- plot3(imedex$simpson.i, imedex$cCasp3,marker="cCasp3")
p3 <- plot3(imedex$simpson.i, imedex$pSTAT1,marker="pSTAT1")
p4 <- plot3(imedex$simpson.i, imedex$yH2AX,marker="yH2AX")
p5 <- plot3(imedex$simpson.i, imedex$Ki67,marker="Ki67")
p6 <- plot3(imedex$simpson.i, imedex$PDL1,marker="PDL1")
p7 <- plot3(imedex$simpson.i, imedex$P21,marker="P21")
p8 <- plot3(imedex$simpson.i, imedex$PD1,marker="PD1")
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 3)