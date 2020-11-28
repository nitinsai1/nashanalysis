#Rscript to plot GSEA enrichment summary plot

#Data input: Gene set name,no.of genes, normalized Enrichment Scores, FDR q val, 

#Load packages
library(ggplot2) 

#Loading the data
data<-read.table("GSEA.txt",sep = "\t", header = T)


# Plotting
p <- ggplot(data, aes(Normalized_enrichment_score, Functions))
p + geom_point(aes(colour=FDR.qvalue, size=Genes)) +
  scale_color_gradientn(colours=rainbow(4), limits=c(0, 1)) +
  geom_vline(xintercept=0, size=0.5, colour="black") +
  theme_bw() + theme(axis.title.x = element_text(size=14,colour = "black"),
                     axis.text.x = element_text(size=14,colour = "black"),
                     panel.border = element_rect(colour = "black", size=1.5),
                     axis.title.y = element_text(size=14,colour = "black"),
                     axis.text.y=element_text(size=14,colour = "black"))+
  #This paramter changes according to the input
  expand_limits(x=c(-2,-1.5,-1,1.5,1,2)) +
  scale_y_discrete(limits=rev(data$Functions))


