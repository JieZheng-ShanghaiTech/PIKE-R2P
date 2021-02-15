#install R packages
install.packages("BiocManager")
install("SingleCellExperiment")
BiocManager::install(c('scater', 'scran', 'uwot', 'Rtsne'))
install.packages('Seurat')
install.packages("compositions")
install.packages("devtools")

library(data.table)
library(BiocManager)
library(SingleCellExperiment)
library(scater)
library(scran)
library(uwot)
library(Rtsne)

#input REAP-seq dataset
rna1<-read.table("E:\\raw_data\\REAP_PBMC\\GSM2685238_mRNA_2_PBMCs_matrix.txt\\GSM2685238_mRNA_2_PBMCs_matrix.txt", header = TRUE, 
               sep = "\t", quote = "\"",
               dec = ".", fill = TRUE, comment.char = "")
adt1<-read.table("E:\\raw_data\\REAP_PBMC\\GSM2685243_protein_2_PBMCs_matrix.txt\\GSM2685243_protein_2_PBMCs_matrix.txt", header = TRUE, 
                sep = "\t", quote = "\"",
                dec = ".", fill = TRUE, comment.char = "")
rna2<-read.table("E:\\raw_data\\REAP_PBMC\\GSM2685239_mRNA_3_PBMCs_matrix.txt\\GSM2685239_mRNA_3_PBMCs_matrix.txt", header = TRUE, 
               sep = "\t", quote = "\"",
               dec = ".", fill = TRUE, comment.char = "")
adt2<-read.table("E:\\raw_data\\REAP_PBMC\\GSM2685244_protein_3_PBMCs_matrix.txt\\GSM2685244_protein_3_PBMCs_matrix.txt", header = TRUE, 
                sep = "\t", quote = "\"",
                dec = ".", fill = TRUE, comment.char = "")
rownames(adt1)=adt1[,1]
adt1=adt1[,-1]
rownames(adt2)=adt2[,1]
adt2=adt2[,-1]

#delete genes that express 0 in all samples
rna1=rna1[which(rowSums(rna1) >0),]
rna2=rna2[which(rowSums(rna2) >0),]
adt1=adt1[which(rowSums(adt1) > 0),]
adt2=adt2[which(rowSums(adt2) > 0),]
cut_rna<-intersect(rownames(rna1), rownames(rna2))
cut_adt<-intersect(rownames(adt1),rownames(adt2))
rna1=rna1[cut_rna,]
rna2=rna2[cut_rna,]
adt1=adt1[cut_adt,]
adt2=adt2[cut_adt,]

#combine two RNA datasets and two ADT datasets
rna=cbind(rna1,rna2)
adt=cbind(adt1,adt2)

#quality control
library(Seurat)
se <- CreateSeuratObject(counts = rna)
se[["percent.mt"]] <- PercentageFeatureSet(object = se, pattern = "^MT-") 
se <- subset(x = se, subset = nFeature_RNA > 250 & percent.mt < 20)

rna_qc=as.matrix(se[["RNA"]]@data)
cut_sample<-intersect(colnames(rna), colnames(rna_qc))
adt_qc=adt[,cut_sample]
adt_qc=as.matrix(adt_qc)
write.csv (rna_qc, file ="E:\\rna_qc.csv")

#denoised
library(devtools)
install_github("jingshuw/SAVERX")
library(SAVERX)
file <- saverx("E:\\rna_qc.csv")
denoised<- readRDS(file)
rna_denoised=denoised$estimate
write.csv (rna_denoised, file ="E:\\rna_denoised.csv")

#output the normalized adt
library(compositions)
adt_clr=clr(adt_qc)
adt_clr=as.data.frame(adt_clr)
write.csv (adt_clr, file ="E:\\adt_clr.csv")




