#R code to analyze the Single cell RNA-seq data analysis
#Input: Data steotosis vs Nash 

#This script does the following: Removal of clusters, reclustering, Find markers, cluster annotation, Feature plots
#loading packages
source("singlecellpackages.R")

#Reading the data
steotosisandnash <- readRDS(file = "steotosis_nash.rds")
#Removing specific clusters other than Tcells 
steotosisandnash1<-subset(steotosisandnash, idents = c("0","1","3","7","10","12","14","15","16","18","19","20"), invert = TRUE)

#Running PCA 
steotosisandnash <- RunPCA(steotosisandnash1, features = VariableFeatures(object=steotosisandnash1))

#Elbowplot analysis
ElbowPlot(steotosisandnash)

# Finding neighbors and clusters
steotosisandnashcluster <- FindNeighbors(steotosisandnash, reduction = "pca",dims = 1:15)
# Finding clusters: we run multiple permutations to allow clustree to analyze optimal clustering resulition
steotosisandnashcluster <- FindClusters(steotosisandnashcluster, resolution = 0)
steotosisandnashcluster <- FindClusters(steotosisandnashcluster, resolution = 0.1)
steotosisandnashcluster <- FindClusters(steotosisandnashcluster, resolution = 0.2)
steotosisandnashcluster <- FindClusters(steotosisandnashcluster, resolution = 0.3)
steotosisandnashcluster <- FindClusters(steotosisandnashcluster, resolution = 0.4)
steotosisandnashcluster <- FindClusters(steotosisandnashcluster, resolution = 0.5)
steotosisandnashcluster <- FindClusters(steotosisandnashcluster, resolution = 0.6)
steotosisandnashcluster <- FindClusters(steotosisandnashcluster, resolution = 0.7)
steotosisandnashcluster <- FindClusters(steotosisandnashcluster, resolution = 0.8)
steotosisandnashcluster <- FindClusters(steotosisandnashcluster, resolution = 0.9)
steotosisandnashcluster <- FindClusters(steotosisandnashcluster, resolution = 1)
#Clustering tree analysis to find resolution
clustree(steotosisandnashcluster)

#Clustering resolution
steotosisandnashcluster <- FindClusters(steotosisandnashcluster, resolution = 1)
head(Idents(steotosisandnashcluster), 5)

#UMAP analysis
steotosisandnashcluster <- RunUMAP(steotosisandnashcluster, dims = 1:15)

#UMAP plot
DimPlot(steotosisandnashcluster, reduction = "umap",label = TRUE)

#Save the data
saveRDS(steotosisandnashcluster, file = "steotosisandnash_Tcellsonly_1res.rds")

#Reading the object for finding markers
steotosisandnash2 <- readRDS(file = "steotosisandnash_Tcellsonly_1res.rds")
steotosisandnash3<-GetAssayData(steotosisandnash2)
#Raw count data
steotosisandnash.c <- as.data.frame(as.matrix(steotosisandnash3))
#Extract the nash raw expression data
fwrite(x = steotosisandnash.c, row.names = TRUE,file = "steotosisandnash_expressionmatrix_Tcells_1res.txt",sep = "\t")

#Find markers
steotosisandnash.markers <- FindAllMarkers(object = steotosisandnash2,logfc.threshold = 0.25,test.use = "wilcox") 
steotosisandnash.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(steotosisandnash.markers, 'steotosisandnash_Allmarkers_Tcellsonly_1res.csv')                                                  







                  
