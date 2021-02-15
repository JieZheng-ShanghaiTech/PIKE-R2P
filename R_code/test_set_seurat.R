install.packages("BiocManager")
install.packages("rlang")
install.packages("ggplot2")
install.packages("SingleR")
install.packages("dplyr")
install.packages("tibble")

library(rlang)
library(Seurat)
library(SingleR)
library(dplyr)
library(tibble)
library(ggplot2)
hpca.se <- HumanPrimaryCellAtlasData()
bpe.se <- BlueprintEncodeData()

#input CITE-seq-CBMC dataset or CITE-seq-PBMC dataset to generate CITE_CBMC_test_set.rds or CITE_PBMC_test_set.rds 
#the test set of CITE-seq dataset
rna_log <- as.sparse(read.csv("E:\\ML_py\\data\\dev_CITE_rna_denoised.csv", header = TRUE,row.names = 1))
rna_log <- t(rna_log)
adt_clr <- as.sparse(read.csv("E:\\ML_py\\data\\dev_CITE_adt_clr.csv", header = TRUE,row.names = 1))
adt_clr <- t(adt_clr)
adt_GNN <- as.sparse(read.csv("E:\\ML_py\\data\\y_pred_GNN_CITE(5).csv", header = TRUE,row.names = 1))
adt_ctp <- as.sparse(read.csv("E:\\ML_py\\data\\y_pred_ctpnet_CITE.csv", header = TRUE,row.names = 1))
#the test set of CITE-seq PBMC dataset
rna_log <- as.sparse(read.csv("E:\\ML_py\\data\\dev_CITE_PBMC_rna_denoised.csv", header = TRUE,row.names = 1))
rna_log <- t(rna_log)
adt_clr <- as.sparse(read.csv("E:\\ML_py\\data\\dev_CITE_PBMC_adt_clr.csv", header = TRUE,row.names = 1))
adt_clr <- t(adt_clr)
adt_GNN <- as.sparse(read.csv("E:\\ML_py\\data\\y_pred_GNN_CITE_PBMC(5).csv", header = TRUE,row.names = 1))
adt_ctp <- as.sparse(read.csv("E:\\ML_py\\data\\y_pred_ctpnet_CITE_PBMC.csv", header = TRUE,row.names = 1))

#RNA normalization
se <- CreateSeuratObject(counts = rna_log)
se <- NormalizeData(se)
se <- FindVariableFeatures(se)
se <- ScaleData(se)
#RNA Run PCA/tsne/umap
se <- RunPCA(se, verbose = FALSE)
se <- RunTSNE(se, dims = 1:25, method = "FIt-SNE")
se <- RunUMAP(object=se,dims=1:25)
#RNA clustering/marker genes
se <- FindNeighbors(se, dims = 1:25)
se <- FindClusters(se, resolution = 0.8)
rna_marker <- FindAllMarkers(se, max.cells.per.ident = 100, min.diff.pct = 0.3, only.pos = TRUE)
new.cluster.ids <- c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16")
names(new.cluster.ids) <- levels(se)
#RNA annotation
se@meta.data$cell.type <- Idents(se)
test <- as.SingleCellExperiment(se)
Anno <- SingleR(test = test,
                ref = list(HP = hpca.se , BP = bpe.se),
                labels = list(hpca.se$label.main , bpe.se$label.main),
                method = "cluster",
                cluster = test$cell.type)
Anno$cluster <- rownames(Anno)
fin <- Anno %>% dplyr::tbl_df() %>% dplyr::select(cluster,labels)
new.cluster.ids <- fin$labels
names(new.cluster.ids) <- levels(se)
se <- RenameIdents(se, new.cluster.ids)
se[["rnaClusterID"]] <- Idents(se)
DimPlot(se,reduction = "umap",label = TRUE) +ggtitle("cell cluster based on RNA")

#add groud_truth
se[["groud_truth"]] <- CreateAssayObject(counts = adt_clr)
se <- ScaleData(se, assay = "groud_truth")
DefaultAssay(se) <- "groud_truth"
adt.data <- GetAssayData(se, slot = "data")
adt.dist <- dist(t(adt.data))
se[["tsne_clr"]] <- RunTSNE(adt.dist, assay = "groud_truth", reduction.name = "UMAP_")
se[["tsne_clr"]] <- RunTSNE(adt.dist, assay = "groud_truth", reduction.name = "UMAP_")
se[["umap_clr"]] <- RunUMAP(object=adt.dist, assay = "groud_truth", reduction.name = "UMAP_")
se[["umap_clr"]] <- RunUMAP(object=adt.dist, assay = "groud_truth", reduction.name = "UMAP_")
se[["clr_snn"]] <- FindNeighbors(adt.dist)$snn
se <- FindClusters(se, resolution = 0.8, graph.name = "clr_snn")
#groud_truth clustering/markers
clustering.table <- table(Idents(se), se$rnaClusterID)
clustering.table
new.cluster.ids <- colnames(clustering.table)[apply(clustering.table, 1, function(x){which.max(x)})]
names(new.cluster.ids) <- levels(se)
se <- RenameIdents(se, new.cluster.ids)
se[["clrClusterID"]] <- Idents(se)
DimPlot(se,reduction = "umap_clr",label = TRUE)+ggtitle("cell cluster based ground turth")

#add PIKE_R2P
se[["PIKE_R2P"]] <- CreateAssayObject(counts = adt_GNN)
se <- ScaleData(se, assay = "PIKE_R2P")
DefaultAssay(se) <- "PIKE_R2P"
adt.data <- GetAssayData(se, slot = "data")
adt.dist <- dist(t(adt.data))
se[["tsne_PIKE_R2P"]] <- RunTSNE(adt.dist, assay = "PIKE_R2P", reduction.name = "UMAP_")
se[["tsne_PIKE_R2P"]] <- RunTSNE(adt.dist, assay = "PIKE_R2P", reduction.name = "UMAP_")
se[["umap_PIKE_R2P"]] <- RunUMAP(object=adt.dist, assay = "PIKE_R2P", reduction.name = "UMAP_")
se[["umap_PIKE_R2P"]] <- RunUMAP(object=adt.dist, assay = "PIKE_R2P", reduction.name = "UMAP_")
se[["PIKE_R2P_snn"]] <- FindNeighbors(adt.dist)$snn
se <- FindClusters(se, resolution = 0.8, graph.name = "PIKE_R2P_snn")
#PIKE_R2P clustering/markers
clustering.table <- table(Idents(se), se$rnaClusterID)
clustering.table
new.cluster.ids <- colnames(clustering.table)[apply(clustering.table, 1, function(x){which.max(x)})]
names(new.cluster.ids) <- levels(se)
se <- RenameIdents(se, new.cluster.ids)
se[["PIKER2PClusterID"]] <- Idents(se)
DimPlot(se,reduction = "umap_PIKE_R2P",label = TRUE)+ggtitle("cell cluster based on PIKE-R2P prediction")

#add cTP-net
se[["CTP"]] <- CreateAssayObject(counts = adt_ctp)
se <- ScaleData(se, assay = "CTP")
DefaultAssay(se) <- "CTP"
adt.data <- GetAssayData(se, slot = "data")
adt.dist <- dist(t(adt.data))
se[["tsne_cTP"]] <- RunTSNE(adt.dist, assay = "CTP", reduction.name = "UMAP_")
se[["tsne_cTP"]] <- RunTSNE(adt.dist, assay = "CTP", reduction.name = "UMAP_")
se[["umap_cTP"]] <- RunUMAP(object=adt.dist, assay = "CTP", reduction.name = "UMAP_")
se[["umap_cTP"]] <- RunUMAP(object=adt.dist, assay = "CTP", reduction.name = "UMAP_")
se[["cTP_snn"]] <- FindNeighbors(adt.dist)$snn
se <- FindClusters(se, resolution = 0.8, graph.name = "cTP_snn")
#cTP-net clustering/markers
clustering.table <- table(Idents(se), se$rnaClusterID)
clustering.table
new.cluster.ids <- colnames(clustering.table)[apply(clustering.table, 1, function(x){which.max(x)})]
names(new.cluster.ids) <- levels(se)
se <- RenameIdents(se, new.cluster.ids)
se[["cTPClusterID"]] <- Idents(se)
DimPlot(se,reduction = "umap_cTP",label = TRUE)+ggtitle("cell cluster based on cTP-Net")
saveRDS(se, file = "E:\\output\\CITE_CBMC_test_set.rds") #or CITE_PBMC_test_set.rds