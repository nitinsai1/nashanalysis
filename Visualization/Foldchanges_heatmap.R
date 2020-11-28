#!/usr/bin/env Rscript

#R Code to plot heatmap using pheatmap package
#Data: fold changes of genes
#Load  packages
library(pheatmap)
library(RColorBrewer)

#Reading data check path
JI<-read.csv("data.csv",header=TRUE,sep=",")

#Row names
row.names(JI)<-JI[,1]

#Delete first column
JI<-JI[,2:dim(JI)[2]]

#Create matrix
JI_matrix<-data.matrix(JI)

#Color coding parameter
#Up and down
colcod<-colorRampPalette(c("blue","ghostwhite","red"))(50)
#Up regulated color codes
#colcod <- colorRampPalette(c("ghostwhite","lightsalmon","lightsalmon","lightcoral","salmon", "red"))(50)
#Downregulated regulated color codes
#colcod <- colorRampPalette(c("blue","dodgerblue","skyblue","lightblue","ghostwhite"))(50)


#Create heatmap
pheatmap(JI_matrix, dendogram='none', legend=TRUE,
         cluster_rows = FALSE, cluster_cols = FALSE,
         clustering_distance_rows = "FALSE",
         color = colcod,
         cellwidth = 18, cellheight=14,fontsize_row=12, fontsize_col=12, 
         border_color="black")



         
         
         
         
         
