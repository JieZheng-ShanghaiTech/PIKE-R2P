#CITE-seq CBMC 
#clustering 
library(Seurat)
library(ggplot2)
se <- readRDS("E:\\output\\CITE_CBMC_test_set.rds")
umap_rnaClusters <- DimPlot(se, reduction = "umap", group.by = "rnaClusterID")
umap_rnaClusters <- umap_rnaClusters + ggtitle("cell cluster based on RNA") + theme(plot.title = element_text(hjust = 0.5))
umap_rnaClusters <- LabelClusters(plot = umap_rnaClusters, id = "rnaClusterID", size = 4)
ggsave("E:\\ML_py\\fig\\fig3(a).png", width = 225, height = 110, units = "mm",dpi = 300)

umap_adtclrClusters <- DimPlot(se, reduction = "umap_clr", group.by = "clrClusterID")
umap_adtclrClusters <- umap_adtclrClusters + ggtitle("cell cluster based ground truth") + theme(plot.title = element_text(hjust = 0.5))
umap_adtclrClusters <- LabelClusters(plot = umap_adtclrClusters, id = "clrClusterID", size = 4)
ggsave("E:\\ML_py\\fig\\fig3(b).png", width = 225, height = 110, units = "mm",dpi = 300)

umap_adtGNNClusters <- DimPlot(se, reduction = "umap_PIKE_R2P", group.by = "PIKER2PClusterID")
umap_adtGNNClusters <- umap_adtGNNClusters + ggtitle("cell cluster based on PIKE-R2P prediction") + theme(plot.title = element_text(hjust = 0.5))
umap_adtGNNClusters <- LabelClusters(plot = umap_adtGNNClusters, id = "PIKER2PClusterID", size = 4)
ggsave("E:\\ML_py\\fig\\fig3(c).png", width = 225, height = 110, units = "mm",dpi = 300)
#CombinePlots(plots = list(umap_rnaClusters,umap_adtclrClusters,umap_adtGNNClusters,umap_adtctpClusters), ncol=2)

#FeaturePlot
DefaultAssay(se) <- "RNA"
titles <- c("groud_truth_CD3","groud_truth_CD4","groud_truth_CD8","groud_truth_CD45RA","groud_truth_CD56","groud_truth_CD16","groud_truth_CD11c","groud_truth_CD14","groud_truth_CD19","groud_truth_CD34")
gene_ids <- c("groudtruth_CD3","groudtruth_CD4","groudtruth_CD8","groudtruth_CD45RA","groudtruth_CD56","groudtruth_CD16","groudtruth_CD11c","groudtruth_CD14","groudtruth_CD19","groudtruth_CD34")
gg_Fig <- FeaturePlot(se, features=gene_ids,min.cutoff = "q05", max.cutoff = "q95")
gg_Fig <- lapply(1:length(gene_ids), function(x) { gg_Fig[[x]] + labs(title=titles[x])})
gg_Fig <- CombinePlots(gg_Fig,ncol=5)
ggsave("E:\\ML_py\\fig\\fig4(a).png", width =400, height = 225, units = "mm",dpi = 300)

titles <- c("PIKE-R2P_CD3","PIKE-R2P_CD4","PIKE-R2P_CD8","PIKE-R2P_CD45RA","PIKE-R2P_CD56","PIKE-R2P_CD16","PIKE-R2P_CD11c","PIKE-R2P_CD14","PIKE-R2P_CD19","PIKE-R2P_CD34")
gene_ids <- c("piker2p_CD3","piker2p_CD4","piker2p_CD8","piker2p_CD45RA","piker2p_CD56","piker2p_CD16","piker2p_CD11c","piker2p_CD14","piker2p_CD19","piker2p_CD34")
gg_Fig <- FeaturePlot(se, features=gene_ids, min.cutoff = "q05", max.cutoff = "q95")
gg_Fig <- lapply(1:length(gene_ids), function(x) { gg_Fig[[x]] + labs(title=titles[x])})
gg_Fig <- CombinePlots(gg_Fig,ncol=5)
ggsave("E:\\ML_py\\fig\\fig4(b).png", width =400, height = 225, units = "mm",dpi = 300)

# CITE-seq PBMC 
#clustering 
library(Seurat)
library(ggplot2)
se <- readRDS("E:\\output\\CITE_PBMC_RNA+clr+gnn+ctp(5)_edit title.rds")
umap_rnaClusters <- DimPlot(se, reduction = "umap", group.by = "rnaClusterID")
umap_rnaClusters <- umap_rnaClusters + ggtitle("cell cluster based on RNA") + theme(plot.title = element_text(hjust = 0.5))
umap_rnaClusters <- LabelClusters(plot = umap_rnaClusters, id = "rnaClusterID", size = 4)
ggsave("E:\\ML_py\\fig\\fig7(a).png", width = 225, height = 110, units = "mm",dpi = 300)

umap_adtclrClusters <- DimPlot(se, reduction = "umap_clr", group.by = "clrClusterID")
umap_adtclrClusters <- umap_adtclrClusters + ggtitle("cell cluster based ground truth") + theme(plot.title = element_text(hjust = 0.5))
umap_adtclrClusters <- LabelClusters(plot = umap_adtclrClusters, id = "clrClusterID", size = 4)
ggsave("E:\\ML_py\\fig\\fig7(b).png", width = 225, height = 110, units = "mm",dpi = 300)

umap_adtGNNClusters <- DimPlot(se, reduction = "umap_PIKE_R2P", group.by = "PIKER2PClusterID")
umap_adtGNNClusters <- umap_adtGNNClusters + ggtitle("cell cluster based on PIKE-R2P prediction") + theme(plot.title = element_text(hjust = 0.5))
umap_adtGNNClusters <- LabelClusters(plot = umap_adtGNNClusters, id = "PIKER2PClusterID", size = 4)
ggsave("E:\\ML_py\\fig\\fig7(c).png", width = 225, height = 110, units = "mm",dpi = 300)

umap_adtCTPClusters <- DimPlot(se, reduction = "umap_cTP", group.by = "cTPClusterID")
umap_adtCTPClusters <- umap_adtCTPClusters + ggtitle("cell cluster based on cTP-net prediction") + theme(plot.title = element_text(hjust = 0.5))
umap_adtCTPClusters <- LabelClusters(plot = umap_adtCTPClusters, id = "cTPClusterID", size = 4)
ggsave("E:\\ML_py\\fig\\fig7(d).png", width = 225, height = 110, units = "mm",dpi = 300)
#CombinePlots(plots = list(umap_rnaClusters,umap_adtclrClusters,umap_adtCTPClusters,umap_adtctpClusters), ncol=2)

#FeaturePlot
library(Seurat)
library(ggplot2)
se <- readRDS("E:\\output\\CITE_PBMC_test_set.rds")
DefaultAssay(se) <- "RNA"
a=FeaturePlot(se, features = c("groudtruth_CD14"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1)+labs(title="groud_truth_CD14")
ggsave("E:\\ML_py\\fig\\fig6(1).png", width =85, height = 75, units = "mm",dpi = 300)
a=FeaturePlot(se, features = c("piker2p_CD14"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1)+labs(title="PIKE-R2P_CD14")
ggsave("E:\\ML_py\\fig\\fig6(2).png", width =85, height = 75, units = "mm",dpi = 300)
a=FeaturePlot(se, features = c("ctp_CD14"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1)+labs(title="cTP-net_CD14")
ggsave("E:\\ML_py\\fig\\fig6(3).png", width =85, height = 75, units = "mm",dpi = 300)
a=FeaturePlot(se, features = c("groudtruth_CD11c"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1)+labs(title="groud_truth_CD11c")
ggsave("E:\\ML_py\\fig\\fig6(4).png", width =85, height = 75, units = "mm",dpi = 300)
a=FeaturePlot(se, features = c("piker2p_CD11c"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1)+labs(title="PIKE-R2P_CD11c")
ggsave("E:\\ML_py\\fig\\fig6(5).png", width =85, height = 75, units = "mm",dpi = 300)
a=FeaturePlot(se, features = c("ctp_CD11c"), min.cutoff = "q05", max.cutoff = "q95", ncol = 1)+labs(title="cTP-net_CD11c")
ggsave("E:\\ML_py\\fig\\fig6(6).png", width =85, height = 75, units = "mm",dpi = 300)




