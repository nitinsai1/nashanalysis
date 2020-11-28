#!/usr/bin/env Rscript
#R script to compute in-degree and out-degree of a transcriptional factor (TF) target network
#Data: .CSV file with two columns TF and target
library(igraph)

#Load network
net<-read.csv("data.csv",sep=",",header=F)
net1<-graph.data.frame(net)
dnet<-as.directed(net1)

#Compute the in-degree

indegree<-degree(dnet,mode = "in")

write.table(indegree,file="indegree.txt",sep="\t")

#Compute the out-degree

outdegree<-degree(dnet,mode = "out")

write.table(outdegree,file="outdegree.txt",sep="\t")