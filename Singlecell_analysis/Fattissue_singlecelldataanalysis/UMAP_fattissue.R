#R code to analyze the Single cell RNA-seq data analysis
#Data: fat tissue
#This script does the following: UMAP and feature dot plot and markers
#loading packages
source("singlecellpackages.R")
library(ggplot2)
library(RColorBrewer)


#Reading the data
fhl <- readRDS(file = "fattissue.rds")
fhl.markers <- FindAllMarkers(object = fhl, min.pct = 0.25,logfc.threshold = 0.25,test.use = "wilcox")
lapply(fhl.markers, head)

#write out markers list
write.csv(fhl.markers,file="fatallmarkers.csv",sep = ",")
#Features
fh<-c("SLC4A10","CCR6","KLRB1","IL2RB","RGS1","ITGA1","CD69","TNFRSF9","CRTAM",
"LAG3","PDCD1","CXCR6","CCL5","CCL4","IFNG","TNF","PRF1","GZMK","GZMB","KLRG1", 
"FOXP3","KLF2","IL7R","CCR7","SELL","CD27", "CD4","CD8B","CD8A","TRGC2","TRGC1",
"TRDC","TRBC2","TRBC1")
#Checkduplidates

sum(duplicated(fh))
fh[fh %in% fh[(duplicated(fh))]]
#Create unique list of markers
f1<-unique(fh)

#UMAP analysis Cluster annotation
current.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10","11")
new.cluster.ids <- c("CD8 CD7","CD8 GZMK","CD4 SOCS3", "CD8 CD74","CD4 IL7R",
                     "CD8 RUNX3","MAIT ","Cycling CD8",
                     "CD8 HSPA1A","Gamma delta tyrobp","CD8 effector","CD4 Treg")
names(new.cluster.ids) <- levels(fhl)
fat<- RenameIdents(fhl, new.cluster.ids)

#UMAP
DimPlot(fat, cols = c("#a6cee3","#FD1212","#b2df8a","#fb9a99","#8073ac","#1a1a1a","#33a02c","#ff7f00","#cab2d6","#878787",
                       "#4A85AE","#FFFF00"))
#Feature dot plot
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100))
#plot
DotPlot(fat, features = f1,dot.scale = 5) + coord_flip() +
  theme(axis.text.x = element_text(angle=90)) + sc



                  
