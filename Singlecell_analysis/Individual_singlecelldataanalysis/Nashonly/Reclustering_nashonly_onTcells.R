#R code to analyze the Single cell RNA-seq data analysis
#Data Nash only dataset
#This script does the following: Removal of clusters, reanalyzing on Tcells clusters, Find Tcells markers, cluster annotation, Feature plots, dot plots
#loading packages
library(Seurat)
library(Rtsne)
library(dplyr)
library(cowplot)
library(DESeq2)
library(MAST)
library(ggplot2)
library(readr)
library(patchwork)
library(data.table)


nash <- readRDS(file = "nashonly.rds")
#Removing specific clusters
nash1<-subset(nash, idents = c("0","1","2","6","8","10","11","12","13","14"), invert = TRUE)

#Running PCA 
nash <- RunPCA(nash1, features = VariableFeatures(object=nash1))

#Elbowplot analysis
ElbowPlot(nash)

#Elbowplot analysis
ElbowPlot(nashonlyliver)

#Clustering
nashonlyliver <- FindNeighbors(nashonlyliver , dims = 1:15)
# Finding neighbors and clusters
nash <- FindNeighbors(nashonlyliver, reduction = "pca",dims = 1:15)
nashonlyliver <- FindClusters(nashonlyliver, resolution = 0)
nashonlyliver <- FindClusters(nashonlyliver, resolution = 0.1)
nashonlyliver <- FindClusters(nashonlyliver, resolution = 0.2)
nashonlyliver <- FindClusters(nashonlyliver, resolution = 0.3)
nashonlyliver <- FindClusters(nashonlyliver, resolution = 0.4)
nashonlyliver <- FindClusters(nashonlyliver, resolution = 0.5)
nashonlyliver <- FindClusters(nashonlyliver, resolution = 0.6)
nashonlyliver<- FindClusters(nashonlyliver, resolution = 0.7)
nashonlyliver <- FindClusters(nashonlyliver, resolution = 0.8)
nashonlyliver <- FindClusters(nashonlyliver, resolution = 0.9)
nashonlyliver <- FindClusters(nashonlyliver, resolution = 1)

#Clustree analysis
clustree(nashonlyliver.combined)
nashonlyliverliver <- FindClusters(nashonlyliver, resolution = 1)
head(Idents(nashonlyliverliver), 5)

#UMAP analysis
nashonlylivercluster <- RunUMAP(nashonlyliver, dims = 1:15)
DimPlot(nashonlylivercluster, reduction = "umap")

#SaveRDS
saveRDS(nashonlylivercluster, file = "Tcells_steotisis_1res.rds")

#Reading object
nash <- readRDS(file = "Tcells_steotisis_1res.rds")

#Feature plot
fe<-c("CD4","CD8B","IL2RB","SLC4A10","CXCR6","RGS1","PDCD1","TRDC","FOXP3","CCR7","TNF","IFNG","GZMK","CCL4","NCAM1","FGFBP2")
FeaturePlot(nash, features = fe) 

#Markers identification
nashonly.markers <- FindAllMarkers(object = nash,logfc.threshold = 0.25,test.use = "wilcox") 
nashonly.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(nash.markers, 'Nashonly_Tcellsmarkers_1res.csv') 

#Annotating clusters

#UMAP analysis Cluster annotation
current.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14")
new.cluster.ids <- c("C0 CD4 TOB1","C1 CD8 RGS1","C2 gd GNLY","C3 MAIT KLF2",
                     "C4 CD8 GNLY FGFBP2","C5 CD8 PLCG2","C6 gd CMC","C7 MAIT NFKBIA",
                     "C8 CD8 CCR7","C9 MAIT PLCG2","C10 NKT KLRC2","C11 CD4 CCR7","C12 CD4 RGS1","C13 gd CMC1 AREG","C14 CD4 FOXP3")
names(new.cluster.ids) <- levels(nash)
nash<- RenameIdents(nash, new.cluster.ids)
#UMAP
DimPlot(nash, cols = c("#228B22","#FD1212","#C6C5C5","#02053E","#4A85AE","#F59191","#cfb53b","#5F64BF","#78D9EE","#874CCF",
                       "#814D06","#FFD300","#034D06","#000000","#FFFF00"))


#Making the feature dot plot for nash only
f<-c("FOXP3","TRGC2","TRGC1","TRDC","TRBC2","TRBC1", "KLF2","LAG3","CCL5","IFNG","PDCD1",
     "IL2RB","SELL","CD69","CXCR6","TNF","PRF1","GZMK","GZMB","KLRG1","IL7R","CCR7",
     "CD4","CD27","CD8B","CRTAM","CD8A","CCL4","RGS1","SLC4A10","CCR6","KLRB1","ITGA1","TNFRSF9")

#Checkduplidates
sum(duplicated(f))
f[f %in% f[(duplicated(f))]]
#Create unique list of markers
f1<-unique(f)

# #Feature dot plot
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100))
DotPlot(nash, features = f1,dot.scale = 5) + coord_flip() + 
  theme(axis.text.x = element_text(angle=90)) + sc


