#install R packages
install.packages("BiocManager")
install("SingleCellExperiment")
BiocManager::install(c('scater', 'scran', 'uwot', 'Rtsne'))
install.packages('Seurat')
install.packages("devtools")

library(data.table)
library(BiocManager)
library(SingleCellExperiment)
library(scater)
library(scran)
library(uwot)
library(Rtsne)

#input CITE-seq-CBMC or CITE-seq-PBMC dataset
rna_CITE <- as.sparse(read.csv("E:\\raw_data\\CITE\\GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv", header = TRUE,row.names = 1))
rna_CITE=as.matrix(rna_CITE)
adt_CITE <- as.sparse(read.csv("E:\\raw_data\\CITE\\GSE100866_CBMC_8K_13AB_10X-ADT_clr-transformed.csv", header = TRUE,row.names = 1))
adt_CITE=as.matrix(adt_CITE)

#delete three low quality proteins, including CCR5,CCR7,CD10
cut_adt=c("CCR5","CCR7","CD10")
cut_adt=setdiff(rownames(adt_CITE),cut_adt)
adt_CITE=adt_CITE[cut_adt,]

# delete all mouse genes except the first 100
rna_CITE=CollapseSpeciesExpressionMatrix(rna_CITE)
#quality control
library(Seurat)
se <- CreateSeuratObject(counts = rna_CITE)
se[["percent.mt"]] <- PercentageFeatureSet(object = se, pattern = "^MT-") 
se <- subset(x = se, subset = nFeature_RNA > 250 & percent.mt < 20)

rna_qc=as.matrix(se[["RNA"]]@data)
cut_sample<-intersect(colnames(rna_CITE), colnames(rna_qc))
adt_qc=adt_CITE[,cut_sample]
adt_qc=as.matrix(adt_qc)
write.csv (rna_qc, file ="E:\\rna_qc.csv")

#denoised
library(devtools)
install_github("jingshuw/SAVERX")
library(SAVERX)
file <- saverx("E:\\rna_qc.csv")
denoised<- readRDS(file)
rna_denoised=denoised$estimate
write.csv (rna_denoised, file ="E:\\rna_denoised_CITE_CBMC.csv")

#output the normalized adt
write.csv (adt_qc, file ="E:\\adt_clr_CITE.csv")

