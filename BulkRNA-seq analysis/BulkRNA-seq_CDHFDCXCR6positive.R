#!/usr/bin/env Rscript
#Program to identify the differentially expressed genes from RNA Seq Count data
#Data: Bulk RNA-seq count data CDHFD liver CXCR6neg vs CDHFD liver CXCR6pos
#Output: Differentially expressed genes, transformed data, volcanoplot

#Load packages
library(DESeq2)
library(BiocGenerics)
library(parallel)
library(RColorBrewer)
library(grDevices)
library(tidyverse)
library(gplots)
library(heatmap3)
library(genefilter)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(ggrepel)



#Load the Count Data
countData<-read.csv("countdata.csv",header=T,row.names=1)
dim(countData)
head(countData)


#Match the Count data and ColData
colData<-DataFrame(condition=factor(c("LiverCDHFDCXCR6-","LiverCDHFDCXCR6-","LiverCDHFDCXCR6-",
                                      "LiverCDHFDCXCR6-","LiverCDHFDCXCR6-","LiverCDHFDCXCR6-",
                                      "LiverCDHFDCXCR6+","LiverCDHFDCXCR6+",
                                      "LiverCDHFDCXCR6+","LiverCDHFDCXCR6+","LiverCDHFDCXCR6+")))
                                     

#Run DESeq2
dds<-DESeqDataSetFromMatrix(countData,colData,formula(~ condition))


#Pre-filtering
dds <- dds[ rowSums(counts(dds)) >= 10,]

#factor levels
colData(dds)$condition<-relevel(colData(dds)$condition,"LiverCDHFDCXCR6-")
dds<-DESeq(dds)

# get differentially expressed genes
res <- results(dds)
# order by BH adjusted p-value
resOrdered <- res[order(res$padj),]
# top of ordered matrix
head(resOrdered)

#Extract DEGS based on with FDR p-value 0.05 and Log fold change 2 (up and down)
sig <- resOrdered[!is.na(resOrdered$padj) & resOrdered$padj<0.05 &
                    abs(resOrdered$log2FoldChange)>=0.6,]
                    
#Write the output: All DEGS
write.table(resOrdered,file="ALL_DEGs_CDHFDCXCR6Pos_vs_CDHFDCXCR6neg.txt",sep="\t")

#Write the output: Significant genes with FDR p-value 0.05 and Log fold change 2 (up and down)
write.table(sig,file="Sig_DEGs_CDHFDCXCR6Pos_vs_CDHFDCXCR6neg.txt",sep="\t")

#Volcanoplot
#Genes highlighted
genes<-c('Ly6c2','Tyrobp','Cd74','Ccr7','S1pr1','Il7r','Siglech','Sell','Tcf7','Ly6d','Klf2','Cxcr6',
         'Nkg7','Sh2d2a','Ccl5','Pdcd1','Gzmb','Cd160','Adgrg1','Tcrg-C4','Osgin1','Chn2','Csf1','Rgs1')
EnhancedVolcano(res,lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 0.6,
                selectLab = genes,
                xlab = bquote(~Log[2]~ 'fold change'),
                pointSize = 1.3,
                labSize = 2.9,
                labhjust = 1.5,
                boxedLabels = TRUE,
                colCustom = keyvals.colour,
                colAlpha = 1,
                legendPosition = 'top',
                legendLabSize = 15,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 0.7,
                colConnectors = 'grey30',
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                border = 'full',
                borderWidth = 1.0,
                borderColour = 'black')

#Transformation
rld <- rlog(dds,blind = TRUE)
de <- rownames(resOrdered)
tmat <- assay(rld)[de,]



#Create the heatmap for matrix of significant genes
de <- rownames(sig)
sigmat <- assay(rld)[de,]
smat<-as.matrix(sigmat)
#Export the transformed all differentially expressed genes
write.table(sigmat,file="Transformeddata_DEGsCDHFDCXCR6Pos_vs_CDHFDCXCR6neg.txt",sep="\t")
