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
hln1 <- readRDS(file = "Tcells_steotisis_1res.rds")

#Markers identification
steotosisonly.markers <- FindAllMarkers(object = hln1,logfc.threshold = 0.25,test.use = "wilcox") 
steotosisonly.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(nash.markers, 'steotisis_Tcellsmarkers_1res.csv') 



