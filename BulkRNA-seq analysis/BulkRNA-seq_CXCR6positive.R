#!/usr/bin/env Rscript
#Script to identify the differentially expressed genes from RNA Seq Count data
#Data: Bulk RNA-seq count data ND Liver CXCR6positive vs CDHFD liver CXCR6pos
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
colData<-DataFrame(condition=factor(c("LiverNDCXCR6+","LiverNDCXCR6+","LiverNDCXCR6+","LiverNDCXCR6+",
                                      "LiverNDCXCR6+","LiverNDCXCR6+","LiverCDHFDCXCR6+","LiverCDHFDCXCR6+",
                                      "LiverCDHFDCXCR6+","LiverCDHFDCXCR6+","LiverCDHFDCXCR6+")))
                                     

#Run DESeq2
dds<-DESeqDataSetFromMatrix(countData,colData,formula(~ condition))


#Pre-filtering
dds <- dds[ rowSums(counts(dds)) >= 10,]

#factor levels
colData(dds)$condition<-relevel(colData(dds)$condition,"LiverNDCXCR6+")
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
write.table(resOrdered,file="ALL_DEGs_CDHFDCXCR6pos_vs_NDCXCR6pos.txt",sep="\t")

#Write the output: Significant genes with FDR p-value 0.05 and Log fold change 2 (up and down)
write.table(sig,file="Sig_DEGs_CDHFDCXCR6pos_vs_NDCXCR6pos.txt",sep="\t")

#Volcanoplot
#Genes highlighted
genes<-c('Ly6c2','Klf2','Il7r','Dnajc15','S1pr1','Tnf','Adgrg5','Sell','Ccr7','Pdcd1',
         'Bhlhe40','Gzmk','Tox','Maf','Ctla4','Adgrg1','Il10','Batf','Tnfrdf4')
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
write.table(sigmat,file="Transformeddata_DEGsCDHFDCXCR6pos_vs_NDCXCR6pos.txt",sep="\t")
