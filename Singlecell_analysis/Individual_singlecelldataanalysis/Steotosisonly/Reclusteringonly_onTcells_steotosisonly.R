#R code to analyze the Single cell RNA-seq data analysis only Tcell clusters
#Data steotosis only
#This script does the following: Removal of clusters, reanalyzing, Find markers only tcells
#loading packages
source("singlecellpackages.R")
#Reading the data
steotosisonlyliver <- readRDS(file = "steotosisonlyliver.rds")
#Removing the clusters
steotosisonlyiver1<-subset(steotosisonlyliver, idents = c("1","2","4","5","10","11","14","15","16","17","19"), invert = TRUE)

#Running PCA 
steotosisonlyliver <- RunPCA(steotosisonlyliver1, features = VariableFeatures(object=steotosisonlyliverliver1))

#Elbowplot analysis
ElbowPlot(steotosisonlyliver)

#Clustering
steotosisonlyliver <- FindNeighbors(steotosisonlyliver , dims = 1:15)
# Finding neighbors and clusters
steotosis <- FindNeighbors(steotosisonlyliver, reduction = "pca",dims = 1:15)
steotosisonlyliver <- FindClusters(steotosisonlyliver, resolution = 0)
steotosisonlyliver <- FindClusters(steotosisonlyliver, resolution = 0.1)
steotosisonlyliver <- FindClusters(steotosisonlyliver, resolution = 0.2)
steotosisonlyliver <- FindClusters(steotosisonlyliver, resolution = 0.3)
steotosisonlyliver <- FindClusters(steotosisonlyliver, resolution = 0.4)
steotosisonlyliver <- FindClusters(steotosisonlyliver, resolution = 0.5)
steotosisonlyliver <- FindClusters(steotosisonlyliver, resolution = 0.6)
steotosisonlyliver<- FindClusters(steotosisonlyliver, resolution = 0.7)
steotosisonlyliver <- FindClusters(steotosisonlyliver, resolution = 0.8)
steotosisonlyliver <- FindClusters(steotosisonlyliver, resolution = 0.9)
steotosisonlyliver <- FindClusters(steotosisonlyliver, resolution = 1)

#Clustree analysis
clustree(steotosisonlyliver.combined)
steotosisonlyliverliver <- FindClusters(steotosisonlyliver, resolution = 1)
head(Idents(steotosisonlyliverliver), 5)

#UMAP analysis
steotosisonlylivercluster <- RunUMAP(steotosisonlyliver, dims = 1:15)
DimPlot(steotosisonlylivercluster, reduction = "umap")


#SaveRDS
saveRDS(steotosisonlylivercluster, file = "Tcells_steotisis_1res.rds")

#Reading RDS
steotosis <- readRDS(file = "Tcells_steotisis_1res.rds")

#Markers identification
steotosisonly.markers <- FindAllMarkers(object = steotosis,logfc.threshold = 0.25,test.use = "wilcox") 
steotosisonly.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(nash.markers, 'steotisis_Tcellsmarkers_1res.csv') 

#UMAP analysis Cluster annotation
current.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14")
new.cluster.ids <- c("C0 CD8 PLCG2","C1 CD4 PLCG2","C2 NKT FCER1G", "C3 CD4 TOB1","C4 CD8 RGS1",
                     "C5 MAIT NFKBIA","C6 CD8 gd GNLY","C7 MAIT CCL20",
                     "C8 CD8 GNLY FGFBP2","C9 CD8 ITGA1","C10 CD4 CCR7","C11 MAIT PLCG2","C12 CD8 KLRC1 FGFBP2","C13 CD4 FOXP3","C14 gd CMC1 AREG")
names(new.cluster.ids) <- levels(steotosis)
steotosis<- RenameIdents(steotosis, new.cluster.ids)
DimPlot(steotosis, cols = c("#F59191","#99d8c9","#ec7014","#228B22","#FD1212","#5F64BF","#C6C5C5","#c51b7d","#4A85AE","#9B1E04","#FFD300","#874CCF","#81D8F6","#FFFF00","#000000"))




