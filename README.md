# nashanalysis
Readme 

Transcriptome analysis 

This repository contains various R scripts used to perform Bulk RNA-seq and single cell RNA-seq datasets presented in Dudek.M et al.

All the scripts are organized in five subfolders

1.	Bulk RNA-seq analysis

This folder contains scripts used to perform mouse bulk RNA-seq data and PCA script. The subfolder (GSEA genesets human to mouse) contains two shell scripts used to convert human geneset.gmt file to mouse geneset.gmt. To execute these commands the human .gmt  file should be present in the folder. 
           1) ConvertingHuman_mouse_GMT.sh (no changes required) 2) Directly execute     following command Convertgmtcommand.sh.

2.	 Transcription factor analysis

This folder contains scripts used to extract the transcription factor binding motifs from promoter sequences (python script). Rscript to compute the in and out degree of TF-target networks.

3.	Visualization 

This folder contains scripts used to plot: 1) GSEA dot plot and 2) heatmap

4.	Single cell analysis

This folder contains collection of scripts to analyze the human single cell data.
1)	Steotosis and nash single cell analysis, 2) Individual single cell nashonly data, 3) fat tissue




 
