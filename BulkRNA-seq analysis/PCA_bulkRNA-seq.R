#!/usr/bin/env Rscript
#R code to perform PCA analysis on all samples

#----Readme of the script----#
#Data: Transformed data from bulk RNA seq data

#PCA 2Dplot

#load neccesary R packages for analysis
library(ggplot2)
library(rgl)
library(FactoMineR)
library(factoextra)
library(scatterplot3d) 
library(devtools)
library(wesanderson)
library(gplots)
library(ggfortify)
library(dplyr)

#Load transformed data of all samples 
data<-read.table("Transformeddata.txt",head=TRUE, sep = "\t")
#row names
row.names(data)<-data[,1]
#uni<-data %>% distinct(Genes,.keep_all = TRUE)
#delete first column
pc<-data[,2:dim(data)[2]]

# create matrix
pc_matrix<-data.matrix(pc)
 
#PCA calculation
pca<-prcomp(t(pc_matrix),scale = FALSE,center = TRUE )
summary(pca)

#Percentage of variance 
percent <- 100*pca$sdev^2/sum(pca$sdev^2)
percent

#Plotting percentage of variance
perc_data <- data.frame(percent=percent, PC=1:length(percent))
ggplot(perc_data, aes(x=PC, y=percent)) + 
  geom_bar(stat="identity") + 
  geom_text(aes(label=round(percent, 2)), size=4, vjust=-.5) + 
  ylim(0, 100) +
  scale_x_continuous(breaks=1:35) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border=element_blank(),
      text = element_text(size=16.8),axis.title.x = element_text(size=23),
      axis.text.x = element_text(size=20),
      axis.title.y = element_text(size=23),
      axis.text.y=element_text(size=20))

#2D PCA Ggplot
#load data 
scores = as.data.frame(pca$x)
names=rownames(scores)
#Plotting 2D scatter plot for PCA 
ggplot(data = scores, aes(x = PC1, y = PC2)) + geom_point(aes(color=names,shape=names),size = 3.9) + 
  scale_shape_manual(values=c(15,15,15,15,15,15,17,17,17,17,17,17,16,16,16,16,16,16,0,0,0,0,0,0,
                              2,2,2,2,2,2,1,1,1,1,1,1))+
  scale_color_manual(values=c('grey0','grey0', 'grey0','grey0','grey0','grey0','blue','blue','blue','blue','blue','blue',
                             'red2','red2','red2','red2','red2','red2','grey0','grey0', 'grey0','grey0','grey0','grey0',
                             'blue','blue','blue','blue','blue','blue','red2','red2','red2','red2','red2','red2')) +
theme_bw()+theme_bw() + theme(legend.position="bottom", legend.box = "horizontal")  +
theme(axis.line = element_line(colour = "black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
      text = element_text(size=16.8),axis.title.x = element_text(size=16),
      axis.text.x = element_text(size=16),
      axis.title.y = element_text(size=16),
      axis.text.y=element_text(size=16),legend.position="none")


