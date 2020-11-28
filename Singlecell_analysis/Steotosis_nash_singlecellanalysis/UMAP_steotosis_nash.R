#R code to analyze the Single cell RNA-seq data analysis
#Data: Steotosis vs nash Tcells 
#This script does the following: UMAP and feature dot plot (Tcells)
#loading packages
source("singlecellpackages.R")
library(ggplot2)
library(RColorBrewer)
#Reading the .rds data
hln1 <- readRDS(file = "steotosisandnash_Tcellsonly_1res.rds")
hln2<-GetAssayData(hln1,assay = "RNA")

#Removal of Interferon genes
counts <- hln2[-(which(rownames(hln2) %in% c('ABCE1','ADAR','BST2','CACTIN','CDC37','CNOT7','DCST1','EGR1','FADD','GBP2','HLA-E','HLA-F','HLA-G','HLA-H','HSP90AB1','IFI27','IFI35','IFI6','IFIT1','IFIT2','IFIT3','IFITM1','IFITM2','IFITM3','IFNA1','IFNA10','IFNA13','IFNA14','IFNA16','IFNA17','IFNA2','IFNA21','IFNA4','IFNA5','IFNA6','IFNA7','IFNA8','IFNAR1','IFNAR2','IFNB1','IKBKE','IP6K2','IRAK1','IRF1','IRF2','IRF3','IRF4','IRF5','IRF6','IRF7','IRF8','IRF9','ISG15','ISG20','JAK1','LSM14A','MAVS','METTL3','MIR21','MMP12','MUL1','MX1','MX2','MYD88','NLRC5','OAS1','OAS2','OAS3','OASL','PSMB8','PTPN1','PTPN11','PTPN2','PTPN6','RNASEL','RSAD2','SAMHD1','SCRIB','SETD2','SHFL','SHMT2','SP100','STAT1','STAT2','TBK1','TREX1','TRIM56','TRIM6','TTLL12','TYK2','UBE2K','USP18','WNT5A','XAF1','YTHDF2','YTHDF3','ZBP1'))),]
hln3 <- subset(hln1, features = rownames(counts))

#Annotating clusters in Umap
new.cluster.ids<-c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19")
new.cluster.ids <- c("CD8 effector 1","CD8 RGS1","CD8 GzmK","NKT 1","CD4 memory 1","CD4 memory 2","NKT 2",
                     "MAIT 1", "CD8 Trm","CD4 RGS1","MAIT 2","Gamma delta 1","MAIT 3","CD8 effector 2",
                     "CD4 naive","MAIT 4","CD8 naive","CD8 memory","Gamma delta 2", "CD4 Treg")
names(new.cluster.ids) <- levels(hln2)
hln4 <- RenameIdents(hln3, new.cluster.ids)

#Markers
f<-c("FOXP3","TRGC2","TRGC1","TRDC","TRBC2","TRBC1", "KLF2","LAG3","CCL5","IFNG","PDCD1",
     "IL2RB","SELL","CD69","CXCR6","TNF","PRF1","GZMK","GZMB","KLRG1","IL7R","CCR7",
     "CD4","CD27","CD8B","CRTAM","CD8A","CCL4","RGS1","SLC4A10","CCR6","KLRB1","ITGA1","TNFRSF9")

#UMAP
DimPlot(hln4, cols = c("4A85AE","#FF3200","#EE8F8F",
                                              "#814D06","#B3F182",
                                              "#52A213","#BB9246","#02053E",
                                              "#8C1803","#034D06","#5F64BF",
                                              "#C6C5C5","#A0A4EC","445C6D",
                                              "#FFD300","#874CCF","#78D9EE","#FF33CE","#000000","#FFFF00"))

#Checkduplidates
sum(duplicated(f))
f[f %in% f[(duplicated(f))]]
#Create unique list of markers
f1<-unique(f)

#Feature dot plot
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100))
DotPlot(hln3, features = f1,dot.scale = 5) + coord_flip() +
theme(axis.text.x = element_text(angle=90)) + sc

