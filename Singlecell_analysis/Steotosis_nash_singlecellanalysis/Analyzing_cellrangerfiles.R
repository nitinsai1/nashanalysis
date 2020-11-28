#R code to analyze the Single cell RNA-seq data analysis
#Input:Cell ranger files 
#Steotosis vs nash data
#This script does the following
#Loading data,integrating dataset, Filtering, QC, Dimensionality, Clustering using UMAP, Marker identification,Feature plots
#loading packages
source("singlecellpackages.R")

#Optional to use all RAM this depends on system
options (future.globals.maxSize = 4000 * 1024^5)

#Reading the data cell ranger output files for merging nash data
nash.data1 <- Read10X(data.dir = "../NASH21_nash/")
nash.data2 <- Read10X(data.dir = "../NASH8_nash/")
nash.data3 <- Read10X(data.dir = "../NASH9_nash/")


#Create data objects
nash.data1 <- CreateSeuratObject(counts = nash.data1,project = "nashdata1")
nash.data2 <- CreateSeuratObject(counts = nash.data2,project = "nashdata2")
nash.data3 <- CreateSeuratObject(counts = nash.data3,project = "nashdata3")

#Combine the files and create the object
nash.combined <- merge(nash.data1, y = c(nash.data2,nash.data3), add.cell.ids = c("8K","8K","8K"), project = "nash8K")
table(nash.combined$orig.ident)
nash1<-GetAssayData(nash.combined)
#Raw count data
nash.c <- as.data.frame(as.matrix(nash1))
#Extract the nash raw expression data
fwrite(x = nash.c, row.names = TRUE,file = "nash_expressionmatrix.txt",sep = "\t")

#Reading the data cell ranger output files for merging data
steotosis.data1 <- Read10X(data.dir = "../steotosis1/")
steotosis.data2 <- Read10X(data.dir = "../steotosis2/")
steotosis.data3 <- Read10X(data.dir = "../steotosis3/")
steotosis.data4 <- Read10X(data.dir = "../steotosis4/")


#Create data objects
steotosis.data1 <- CreateSeuratObject(counts = steotosis.data1,project = "steotosisdata1")
steotosis.data2 <- CreateSeuratObject(counts = steotosis.data2,project = "steotosisdata2")
steotosis.data3 <- CreateSeuratObject(counts = steotosis.data3,project = "steotosisdata3")
steotosis.data4 <- CreateSeuratObject(counts = steotosis.data4,project = "steotosisdata4")


#Combine the files and create the object
steotosis.combined <- merge(steotosis.data1, y = c(steotosis.data2,steotosis.data3,steotosis.data4), add.cell.ids = c("8K","8K","8K","8K"), project = "steotosis8K")
table(steotosis.combined$orig.ident)
steotosis1<-GetAssayData(steotosis.combined)
#Raw count data
control.c <- as.data.frame(as.matrix(steotosis1))
#Extract the bnash expression data
fwrite(x = control.c, row.names = TRUE,file = "steotosis_expressionmatrix.txt",sep = "\t")


#Reading the raw expression matrix
steotosisliver.data <- read.table(file = "steotosis_expressionmatrix.txt", header = TRUE, row.names = 1, sep = "\t", as.is = TRUE)
nash.data <- read.table(file = "nash_expressionmatrix.txt", header = TRUE, row.names = 1, sep = "\t", as.is = TRUE)

#steotosis liver object (Initialize the Seurat object with the raw combined file (non-normalized data))
steotosisliver <- CreateSeuratObject(counts = steotosisliver.data, project = "steotosis_liver", min.cells = 3,min.features = 200)

#QC analysis steotosis liver
steotosisliver[["percent.mt"]] <- PercentageFeatureSet(steotosisliver, pattern = "^MT-")
VlnPlot(steotosisliver, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
plot1 <- FeatureScatter(steotosisliver, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(steotosisliver, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#Subset the RNA features and mitochondria
steotosisliver <- subset(steotosisliver, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 20)
countmath<-GetAssayData(steotosisliver,assay = "RNA")
#Load black listed genes
#Oxidative stress
countsh <- countmath[-(which(rownames(countmath) %in% c('ABCC2','ABCD1','ABL1','ACOX2','ADA','ADAM9','ADIPOQ','ADNP2','ADPRHL2','AGAP3','AIF1','AIFM1','AKR1C3','AKT1','ALAD','ALDH3B1','ALS2','ANGPTL7','ANKRD2','ANKZF1','ANXA1','APEX1','APOA4','APOD','APOE','APP','APTX','AQP1','AREG','ARG1','ARL6IP5','ARNT','ARNTL','ATF4','ATG7','ATOX1','ATP13A2','ATP2A2','ATP7A','ATRN','AXL','BAD','BAG5','BAK1','BCL2','BECN1','BMP7','BNIP3','BRF2','BTG1','BTK','C19orf12','CA3','CAMKK2','CAPN2','CASP3','CAT','CBX8','CCL19','CCNA2','CCR7','CCS','CD36','CD38','CDK1','CDK2','CFLAR','CHD6','CHRNA4','CHUK','CLN8','COA8','COL1A1','CPEB2','CRK','CRYAB','CRYGD','CTNNB1','CYBA','CYBB','CYCS','CYGB','CYP1B1','CYP2E1','DAPK1','DGKK','DHCR24','DHFR','DHFRP1','DHRS2','DIABLO','DNM2','DPEP1','DUOX1','DUOX2','DUSP1','ECT2','EDN1','EEF2','EGFR','EGLN1','EIF2S1','ENDOG','EPAS1','EPX','ERCC1','ERCC2','ERCC3','ERCC6','ERCC6L2','ERCC8','ERO1A','ETFDH','ETS1','ETV5','EZH2','FABP1','FANCC','FANCD2','FBLN5','FBXO7','FBXW7','FER','FGF8','FKBP1B','FOS','FOSL1','FOXO1','FOXO3','FUT8','FXN','FYN','FZD1','G6PD','GATA4','GCH1','GCLC','GCLM','GGT7','GJB2','GLRX2','GNAO1','GPR37','GPR37L1','GPX1','GPX2','GPX3','GPX4','GPX5','GPX6','GPX7','GPX8','GSKIP','GSR','GSS','GSTP1','GUCY1B1','HAO1','HBA1','HBA2','HBB','HDAC2','HDAC6','HGF','HIF1A','HMOX1','HMOX2','HNRNPD','HNRNPM','HP','HSF1','HSPA1A','HSPA1B','HSPB1','HTRA2','HYAL1','HYAL2','IDH1','IL10','IL18RAP','IL6','IMPACT','INS','IPCEF1','JAK2','JUN','KCNA5','KCNC2','KDM6B','KEAP1','KLF2','KLF4','KPNA4','KRT1','LANCL1','LCN2','LDHA','LIAS','LONP1','LPO','LRRK2','MACROH2A1','MAP1LC3A','MAP3K5','MAPK1','MAPK13','MAPK3','MAPK7','MAPK8','MAPK9','MAPKAP1','MAPT','MB','MBL2','MCL1','MCTP1','MDM2','MEAK7','MELK','MET','MGAT3','MGST1','MICB','MIR103A1','MIR103A2','MIR107','MIR132','MIR133A1','MIR133A2','MIR17','MIR195','MIR19A','MIR21','MIR29B1','MIR29B2','MIR34A','MIR92A1','MIR92A2','MIRLET7B','MMP14','MMP2','MMP3','MMP9','MPO','MPV17','MSRA','MSRB2','MSRB3','MT-CO1','MT-ND1','MT-ND3','MT-ND5','MT-ND6','MT3','MTF1','MTR','MYB','MYEF2','NAPRT','NCF1','NCF2','NCF4','NCOA7','NDUFA12','NDUFA6','NDUFB4','NDUFS2','NDUFS8','NEIL1','NET1','NFE2L1','NFE2L2','NME2','NME5','NME8','NOL3','NONO','NOS3','NOX1','NOX4','NOX5','NQO1','NR4A2','NR4A3','NUDT1','NUDT15','NUDT2','OGG1','OSER1','OXR1','OXSR1','P4HB','PAGE4','PARK7','PARP1','PAWR','PAX2','PCGF2','PCNA','PDCD10','PDE8A','PDGFD','PDGFRA','PDGFRB','PDK1','PDK2','PDLIM1','PENK','PINK1','PKD2','PLA2R1','PLEKHA1','PLK3','PML','PNKP','PNPT1','PON2','PPARGC1A','PPARGC1B','PPIF','PPP1R15B','PPP2CB','PPP5C','PRDX1','PRDX2','PRDX3','PRDX4','PRDX5','PRDX6','PRKAA1','PRKAA2','PRKCD','PRKD1','PRKN','PRKRA','PRNP','PRODH','PRR5L','PSAP','PSEN1','PSIP1','PSMB5','PTGS1','PTGS2','PTK2B','PTPRK','PTPRN','PXDN','PXDNL','PXN','PYCR1','PYCR2','PYROXD1','RACK1','RAD52','RBM11','RBPMS','RELA','REST','RGN','RGS14','RHOB','RIPK1','RIPK3','RNF112','ROMO1','RPS3','RRM2B','RWDD1','S100A7','SCARA3','SCGB1A1','SDC1','SELENOK','SELENON','SELENOP','SELENOS','SESN1','SESN2','SESN3','SETX','SFPQ','SGK2','SIGMAR1','SIN3A','SIRPA','SIRT1','SIRT2','SLC23A2','SLC25A24','SLC7A11','SLC8A1','SMPD3','SNCA','SOD1','SOD2','SOD3','SP1','SRC','SRXN1','STAR','STAT1','STAT6','STAU1','STC2','STK24','STK25','STK26','STOX1','STX2','STX4','SZT2','TAT','TBC1D24','TLR4','TLR6','TMEM161A','TNF','TNFAIP3','TOR1A','TP53','TP53INP1','TPM1','TPO','TRAF2','TRAP1','TREX1','TRPA1','TRPC6','TRPM2','TSC1','TXN','TXN2','TXNIP','TXNRD1','TXNRD2','UBE3A','UBQLN1','UCN','UCP1','UCP3','VKORC1L1','VNN1','VRK2','WNT1','WNT16','WRN','XPA','XRCC1','ZC3H12A','ZNF277','ZNF580','ZNF622'))),]
#mitochondria genes
mito.genes <- grep(pattern = "^MT-", rownames(x = steotosisliver), value = TRUE)
#ribosomal genes
ribo.genes <- grep(pattern = "^RP[SL]", rownames(x = steotosisliver), value = TRUE)
#mito genes removal
countsh <- countsh[-(which(rownames(countsh) %in% c('MT-ND1','MT-ND2','MT-CO1','MT-CO2','MT-ATP8','MT-ATP6','MT-CO3','MT-ND3','MT-ND4L','MT-ND4','MT-ND5','MT-ND6','MT-CYB'))),]
#ribo genes removal
countsh <- countsh[-(which(rownames(countsh) %in% c('RPL22','RPL11','RPS6KA1','RPS8','RPL5','RPS27','RPS6KC1','RPS7','RPS27A','RPL31','RPL37A','RPL32','RPL15','RPSA','RPL14','RPL29','RPL24','RPL22L1','RPL39L','RPL35A','RPL9','RPL34-AS1','RPL34','RPS3A','RPL37','RPS23','RPS14','RPL26L1','RPS18','RPS10','RPL10A','RPL7L1','RPS12','RPS6KA2','RPS6KA3','RPS4X','RPS6KA6','RPL36A','RPL39','RPL10','RPS20','RPL7','RPL30','RPL8','RPS6','RPL35','RPL12','RPL7A','RPLP2','RPL27A','RPS13','RPS6KA4','RPS6KB2','RPS3','RPS25','RPS24','RPS26','RPL41','RPL6','RPLP0','RPL21','RPL10L','RPS29','RPL36AL','RPS6KL1','RPS6KA5','RPS27L','RPL4','RPLP1','RPS17','RPL3L','RPS2','RPS15A','RPL13','RPL26','RPL23A','RPL23','RPL19','RPL27','RPS6KB1','RPL38','RPL17-C18orf32','RPL17','RPS21','RPS15','RPL36','RPS28','RPL18A','RPS16','RPS19','RPL18','RPL13A','RPS11','RPS9','RPL28','RPS5','RPL3','RPS19BP1'))),]
#Subsetting data after removal of black listed genes
steotosisliver <- subset(steotosisliver, features = rownames(countsh))


#NASH object
nash <- CreateSeuratObject(counts = nash.data, project = "Nash", min.cells = 3,min.features = 200)
#QC analysis steotosis liver
nash[["percent.mt"]] <- PercentageFeatureSet(nash, pattern = "^MT-")
VlnPlot(nash, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
plot1 <- FeatureScatter(nash, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(nash, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#Subset the RNA features and mitochondria
nash <- subset(nash, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mt < 20)
countnash<-GetAssayData(nash,assay = "RNA")
#Load black listed genes
#Oxidative stress
countnash <- countnash[-(which(rownames(countnash) %in% c('ABCC2','ABCD1','ABL1','ACOX2','ADA','ADAM9','ADIPOQ','ADNP2','ADPRHL2','AGAP3','AIF1','AIFM1','AKR1C3','AKT1','ALAD','ALDH3B1','ALS2','ANGPTL7','ANKRD2','ANKZF1','ANXA1','APEX1','APOA4','APOD','APOE','APP','APTX','AQP1','AREG','ARG1','ARL6IP5','ARNT','ARNTL','ATF4','ATG7','ATOX1','ATP13A2','ATP2A2','ATP7A','ATRN','AXL','BAD','BAG5','BAK1','BCL2','BECN1','BMP7','BNIP3','BRF2','BTG1','BTK','C19orf12','CA3','CAMKK2','CAPN2','CASP3','CAT','CBX8','CCL19','CCNA2','CCR7','CCS','CD36','CD38','CDK1','CDK2','CFLAR','CHD6','CHRNA4','CHUK','CLN8','COA8','COL1A1','CPEB2','CRK','CRYAB','CRYGD','CTNNB1','CYBA','CYBB','CYCS','CYGB','CYP1B1','CYP2E1','DAPK1','DGKK','DHCR24','DHFR','DHFRP1','DHRS2','DIABLO','DNM2','DPEP1','DUOX1','DUOX2','DUSP1','ECT2','EDN1','EEF2','EGFR','EGLN1','EIF2S1','ENDOG','EPAS1','EPX','ERCC1','ERCC2','ERCC3','ERCC6','ERCC6L2','ERCC8','ERO1A','ETFDH','ETS1','ETV5','EZH2','FABP1','FANCC','FANCD2','FBLN5','FBXO7','FBXW7','FER','FGF8','FKBP1B','FOS','FOSL1','FOXO1','FOXO3','FUT8','FXN','FYN','FZD1','G6PD','GATA4','GCH1','GCLC','GCLM','GGT7','GJB2','GLRX2','GNAO1','GPR37','GPR37L1','GPX1','GPX2','GPX3','GPX4','GPX5','GPX6','GPX7','GPX8','GSKIP','GSR','GSS','GSTP1','GUCY1B1','HAO1','HBA1','HBA2','HBB','HDAC2','HDAC6','HGF','HIF1A','HMOX1','HMOX2','HNRNPD','HNRNPM','HP','HSF1','HSPA1A','HSPA1B','HSPB1','HTRA2','HYAL1','HYAL2','IDH1','IL10','IL18RAP','IL6','IMPACT','INS','IPCEF1','JAK2','JUN','KCNA5','KCNC2','KDM6B','KEAP1','KLF2','KLF4','KPNA4','KRT1','LANCL1','LCN2','LDHA','LIAS','LONP1','LPO','LRRK2','MACROH2A1','MAP1LC3A','MAP3K5','MAPK1','MAPK13','MAPK3','MAPK7','MAPK8','MAPK9','MAPKAP1','MAPT','MB','MBL2','MCL1','MCTP1','MDM2','MEAK7','MELK','MET','MGAT3','MGST1','MICB','MIR103A1','MIR103A2','MIR107','MIR132','MIR133A1','MIR133A2','MIR17','MIR195','MIR19A','MIR21','MIR29B1','MIR29B2','MIR34A','MIR92A1','MIR92A2','MIRLET7B','MMP14','MMP2','MMP3','MMP9','MPO','MPV17','MSRA','MSRB2','MSRB3','MT-CO1','MT-ND1','MT-ND3','MT-ND5','MT-ND6','MT3','MTF1','MTR','MYB','MYEF2','NAPRT','NCF1','NCF2','NCF4','NCOA7','NDUFA12','NDUFA6','NDUFB4','NDUFS2','NDUFS8','NEIL1','NET1','NFE2L1','NFE2L2','NME2','NME5','NME8','NOL3','NONO','NOS3','NOX1','NOX4','NOX5','NQO1','NR4A2','NR4A3','NUDT1','NUDT15','NUDT2','OGG1','OSER1','OXR1','OXSR1','P4HB','PAGE4','PARK7','PARP1','PAWR','PAX2','PCGF2','PCNA','PDCD10','PDE8A','PDGFD','PDGFRA','PDGFRB','PDK1','PDK2','PDLIM1','PENK','PINK1','PKD2','PLA2R1','PLEKHA1','PLK3','PML','PNKP','PNPT1','PON2','PPARGC1A','PPARGC1B','PPIF','PPP1R15B','PPP2CB','PPP5C','PRDX1','PRDX2','PRDX3','PRDX4','PRDX5','PRDX6','PRKAA1','PRKAA2','PRKCD','PRKD1','PRKN','PRKRA','PRNP','PRODH','PRR5L','PSAP','PSEN1','PSIP1','PSMB5','PTGS1','PTGS2','PTK2B','PTPRK','PTPRN','PXDN','PXDNL','PXN','PYCR1','PYCR2','PYROXD1','RACK1','RAD52','RBM11','RBPMS','RELA','REST','RGN','RGS14','RHOB','RIPK1','RIPK3','RNF112','ROMO1','RPS3','RRM2B','RWDD1','S100A7','SCARA3','SCGB1A1','SDC1','SELENOK','SELENON','SELENOP','SELENOS','SESN1','SESN2','SESN3','SETX','SFPQ','SGK2','SIGMAR1','SIN3A','SIRPA','SIRT1','SIRT2','SLC23A2','SLC25A24','SLC7A11','SLC8A1','SMPD3','SNCA','SOD1','SOD2','SOD3','SP1','SRC','SRXN1','STAR','STAT1','STAT6','STAU1','STC2','STK24','STK25','STK26','STOX1','STX2','STX4','SZT2','TAT','TBC1D24','TLR4','TLR6','TMEM161A','TNF','TNFAIP3','TOR1A','TP53','TP53INP1','TPM1','TPO','TRAF2','TRAP1','TREX1','TRPA1','TRPC6','TRPM2','TSC1','TXN','TXN2','TXNIP','TXNRD1','TXNRD2','UBE3A','UBQLN1','UCN','UCP1','UCP3','VKORC1L1','VNN1','VRK2','WNT1','WNT16','WRN','XPA','XRCC1','ZC3H12A','ZNF277','ZNF580','ZNF622'))),]
#mitochondria genes
mito.genes <- grep(pattern = "^MT-", rownames(x = nash), value = TRUE)
#ribosomal genes
ribo.genes <- grep(pattern = "^RP[SL]", rownames(x = nash), value = TRUE)
#mito genes removal
countnash <- countnash[-(which(rownames(countnash) %in% c('MT-ND1','MT-ND2','MT-CO1','MT-CO2','MT-ATP8','MT-ATP6','MT-CO3','MT-ND3','MT-ND4L','MT-ND4','MT-ND5','MT-ND6','MT-CYB'))),]
#ribo genes removal
countnash <- countnash[-(which(rownames(countnash) %in% c('RPL22','RPL11','RPS6KA1','RPS8','RPL5','RPS27','RPS6KC1','RPS7','RPS27A','RPL31','RPL37A','RPL32','RPL15','RPSA','RPL14','RPL29','RPL24','RPL22L1','RPL39L','RPL35A','RPL9','RPL34-AS1','RPL34','RPS3A','RPL37','RPS23','RPS14','RPL26L1','RPS18','RPS10','RPL10A','RPL7L1','RPS12','RPS6KA2','RPS6KA3','RPS4X','RPS6KA6','RPL36A','RPL39','RPL10','RPS20','RPL7','RPL30','RPL8','RPS6','RPL35','RPL12','RPL7A','RPLP2','RPL27A','RPS13','RPS6KA4','RPS6KB2','RPS3','RPS25','RPS24','RPS26','RPL41','RPL6','RPLP0','RPL21','RPL10L','RPS29','RPL36AL','RPS6KL1','RPS6KA5','RPS27L','RPL4','RPLP1','RPS17','RPL3L','RPS2','RPS15A','RPL13','RPL26','RPL23A','RPL23','RPL19','RPL27','RPS6KB1','RPL38','RPL17-C18orf32','RPL17','RPS21','RPS15','RPL36','RPS28','RPL18A','RPS16','RPS19','RPL18','RPL13A','RPS11','RPS9','RPL28','RPS5','RPL3','RPS19BP1'))),]
#Subsetting data after removal of black listed genes
nash <- subset(nash, features = rownames(countnash))
#merging
steotosisliver@meta.data$treatment<-"STEOTOSIS"
nash@meta.data$treatment<-"NASH"
steotosisliverandnash <- merge(steotosisliver, y = c(nash), add.cell.ids = c("steotosis", "nash"), project = "steotosislivernash")

#Splitting the data
steotosisliverandnash <- SplitObject (steotosisliverandnash, split.by = "treatment")
steotosisliverandnash <- steotosisliverandnash[c("STEOTOSIS", "NASH")]

#Normalization and finding variable features
steotosisliverandnash <- lapply(X = steotosisliverandnash, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

#Integration analysis
steotosisandnash.anchors <- FindIntegrationAnchors(object.list = steotosisliverandnash, dims = 1:20)
steotosisandnash.combined <- IntegrateData(anchorset = steotosisandnash.anchors, dims = 1:20)
DefaultAssay(steotosisandnash.combined) <- "integrated"

# Run the standard workflow scaling and cell cycle regresion
steotosisandnash.combined <- ScaleData(steotosisandnash.combined, verbose = FALSE)
#Remove cell cycle effect
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#Cell cycle scoring
steotosisandnash.combined <- CellCycleScoring(steotosisandnash.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#Regressing cell cycle
steotosisandnash.combined <- ScaleData(steotosisandnash.combined, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(steotosisandnash.combined))

# Run PCA and UMAPclustering analysis
steotosisandnash.combined <- RunPCA(steotosisandnash.combined, npcs = 30, verbose = FALSE)

# Elbowplot
ElbowPlot(steotosisandnash.combined)
# Finding neighbors and clusters
steotosisandnash.combined <- FindNeighbors(steotosisandnash.combined, reduction = "pca",dims = 1:15)
steotosisandnash.combined <- FindClusters(steotosisandnash.combined, resolution = 0)
steotosisandnash.combined <- FindClusters(steotosisandnash.combined, resolution = 0.1)
steotosisandnash.combined <- FindClusters(steotosisandnash.combined, resolution = 0.2)
steotosisandnash.combined <- FindClusters(steotosisandnash.combined, resolution = 0.3)
steotosisandnash.combined <- FindClusters(steotosisandnash.combined, resolution = 0.4)
steotosisandnash.combined <- FindClusters(steotosisandnash.combined, resolution = 0.5)
steotosisandnash.combined <- FindClusters(steotosisandnash.combined, resolution = 0.6)
steotosisandnash.combined <- FindClusters(steotosisandnash.combined, resolution = 0.7)
steotosisandnash.combined <- FindClusters(steotosisandnash.combined, resolution = 0.8)
steotosisandnash.combined <- FindClusters(steotosisandnashcluster, resolution = 0.9)
steotosisandnash.combined <- FindClusters(steotosisandnash.combined, resolution = 1)

#Clustree analysis
clustree(steotosisandnash.combined)
#Find clusters
steotosisandnash.combined <- FindClusters(steotosisandnash.combined, resolution = 0.6)
#UMAP analysis
steotosisandnash.combined <- RunUMAP(steotosisandnash.combined, dims = 1:15)
#UMAP splitting according to condtions
p1 <- DimPlot(steotosisandnash.combined, reduction = "umap", group.by = "treatment")
p2 <- DimPlot(steotosisandnash.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
#UMAP
DimPlot(steotosisandnash.combined, reduction = "umap", split.by = "treatment")

#Marker identification
steotosisandnash.combined.markers <- FindAllMarkers(object = steotosisandnash.combined,logfc.threshold = 0.25,test.use = "wilcox") 
steotosisandnash.combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(steotosisandnash.combined.markers, 'steotosisandnash.combined_Allmarkers.csv')      

#feature plot with specific markers
f<-c("FOXP3","TRGC2","TRGC1","TRDC","TRBC2","TRBC1", "KLF2","LAG3","CCL5","IFNG","PDCD1",
     "IL2RB","SELL","CD69","CXCR6","TNF","PRF1","GZMK","GZMB","KLRG1","IL7R","CCR7",
     "CD4","CD27","CD8B","CRTAM","NCAM1")

#Checkduplidates
sum(duplicated(f))
f[f %in% f[(duplicated(f))]]
#Create unique list of markers
f1<-unique(f)

#Feature plot
DotPlot(hln, features = rev(f1), cols = c("blue", "red"), dot.scale = 4, 
        split.by = "treatment") + RotatedAxis()

#Save the data
saveRDS(steotosisandnash.combined, file = "steotosis_nash.rds")


