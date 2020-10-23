#pipeline prepared by Jessica S. Sadick
#majority of pipeline based on code originally deposited by https://github.com/HelenaLC and published in doi.org/10.1101/713412

#Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

#Load astrocyte-only SeuratObjects for each dataset
##Comparing current study with AD mice models
so_Hasel_astro <- readRDS(file.path("file_path", "Hasel_so_astro.rds"))
so_Zhou_7mo_astro <- readRDS(file.path("file_path", "Zhou_mouse_CTX_7mo_so_astro_final_15PC.rds"))
so_Zhou_15mo_astro <- readRDS(file.path("file_path", "Zhou_mouse_CTX_15mo_so_astro_final_10PC.rds"))
so_Habib_astro <- readRDS(file.path("file_path", "Habib_Ms_CTX_so_astro_PC10.rds"))

#Add dataset variable in metadata
Meta_Hasel <- so_Hasel_astro@meta.data
Meta_Hasel["dataset"] <- c("Hasel")
NewMeta_Hasel <- subset(Meta_Hasel, select = c("dataset"))
so_Hasel_astro <- AddMetaData(so_Hasel_astro, NewMeta_Hasel)
head(x = so_Hasel_astro[[]])

Meta_Zhou_7mo <- so_Zhou_7mo_astro@meta.data
Meta_Zhou_7mo["dataset"] <- c("Zhou_7mo")
NewMeta_Zhou_7mo <- subset(Meta_Zhou_7mo, select = c("dataset"))
so_Zhou_7mo_astro <- AddMetaData(so_Zhou_7mo_astro, NewMeta_Zhou_7mo)
head(x = so_Zhou_7mo_astro[[]])

Meta_Zhou_15mo <- so_Zhou_15mo_astro@meta.data
Meta_Zhou_15mo["dataset"] <- c("Zhou_15mo")
NewMeta_Zhou_15mo <- subset(Meta_Zhou_15mo, select = c("dataset"))
so_Zhou_15mo_astro <- AddMetaData(so_Zhou_15mo_astro, NewMeta_Zhou_15mo)
head(x = so_Zhou_15mo_astro[[]])

Meta_Habib <- so_Habib_astro@meta.data
Meta_Habib["dataset"] <- c("Habib")
NewMeta_Habib <- subset(Meta_Habib, select = c("dataset"))
so_Habib_astro <- AddMetaData(so_Habib_astro, NewMeta_Habib)
head(x = so_Habib_astro[[]])

#Merge SeuratObjects
so_astro_merge <- merge(x = so_Hasel_astro, y = c(so_Zhou_7mo_astro, so_Zhou_15mo_astro, so_Habib_astro))

#Remove original integrated data associated with object
DefaultAssay(so_astro_merge) <- "RNA"
so_astro_merge[['integrated']] <- NULL 

#RECLUSTER MULTI-DATASET ASTROCYTE-SPECIFIC SeuratObject

#Increase future's maximum allowed size of objects
options(future.globals.maxSize = 2048 * 1024 ^20)
memory.limit(size = 1000000)

#Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

#Convert SeuratObject to SCE
sce <- as.SingleCellExperiment(so_astro_merge, assay = "RNA")

#INTEGRATE
#Split by sample
cells_by_sample <- split(colnames(sce), sce$sample_id)
so_astro_merge <- lapply(cells_by_sample, function(i)
  SubsetData(so_astro_merge, cells = i))

#Normalize, find variable genes, and scale
so_astro_merge <- lapply(so_astro_merge, NormalizeData, verbose = FALSE)
so_astro_merge <- lapply(so_astro_merge, FindVariableFeatures, nfeatures = 2e3,
                   selection.method = "vst", do.plot = FALSE, verbose = FALSE)
so_astro_merge <- lapply(so_astro_merge, ScaleData, verbose = FALSE)

#Find anchors and integrate
##Decrease k.filter to minimize number of astrocytes identified per sample across all datasets
as <- FindIntegrationAnchors(so_astro_merge, verbose = FALSE, k.filter = 90)
so_astro_merge <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = FALSE)

#Scale integrated data
DefaultAssay(so_astro_merge) <- "integrated"
so_astro_merge <- ScaleData(so_astro_merge, display.progress = FALSE)

#DIMENSION REDUCTION
so_astro_merge <- RunPCA(so_astro_merge, npcs = 50, verbose = FALSE)
ElbowPlot(so_astro_merge, ndims = 50)
##Update number of PCs used
so_astro_merge <- RunTSNE(so_astro_merge, reduction = "pca", dims = seq_len(20),
                    seed.use = 1, do.fast = TRUE, verbose = FALSE, check_duplicates = FALSE)
so_astro_merge <- RunUMAP(so_astro_merge, reduction = "pca", dims = seq_len(20),
                    seed.use = 1, verbose = FALSE)

#CLUSTERING
so_astro_merge <- FindNeighbors(so_astro_merge, reduction = "pca", dims = seq_len(20), verbose = FALSE)
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
  so_astro_merge <- FindClusters(so_astro_merge, resolution = res, random.seed = 1, verbose = FALSE)

#DR COLORED BY SAMPLE, GROUP, AND CLUSTER ID
thm <- theme(aspect.ratio = 1, legend.position = "none")
ps <- lapply(c("sample_id", "group_id", "ident"), function(u) {
  p1 <- DimPlot(so_astro_merge, reduction = "tsne", group.by = u) + thm
  p2 <- DimPlot(so_astro_merge, reduction = "umap", group.by = u)
  lgd <- get_legend(p2)
  p2 <- p2 + thm
  list(p1, p2, lgd)
  plot_grid(p1, p2, lgd, nrow = 1,
            rel_widths = c(1, 1, 0.5))
})
plot_grid(plotlist = ps, ncol = 1)

#SAVE SeuratObject
saveRDS(so_astro_merge, file.path("file_path", "so_astro_merge_HHZ_PC20.rds"))

#---------------------------------------------------------------------------------------------------
#CLUSTER ANNOTATION FOR MULTI-DATASET ASTROCYTE ONLY SeuratObject

#Load packages
library(ComplexHeatmap)
library(cowplot)
library(ggplot2)
library(dplyr)
library(purrr)
library(RColorBrewer)
library(viridis)
library(scran)
library(Seurat)
library(SingleCellExperiment)

#Load data and convert to SCE
so_astro_merge <- readRDS(file.path("file_path", "so_astro_merge_HHZ_PC20.rds"))
sce <- as.SingleCellExperiment(so_astro_merge, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>% 
  mutate_if(is.character, as.factor) %>% 
  DataFrame(row.names = colnames(sce))

#Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)

so_astro_merge <- SetIdent(so_astro_merge, value = "integrated_snn_res.0.3")
so_astro_merge@meta.data$cluster_id <- Idents(so_astro_merge)
sce$cluster_id <- Idents(so_astro_merge)
(n_cells <- table(sce$cluster_id, sce$sample_id))
write.csv(table(sce$cluster_id, sce$sample_id), "file_path/so_astro_merge_HHZ_numbers.csv")

nk <- length(kids <- set_names(levels(sce$cluster_id)))
ns <- length(sids <- set_names(levels(sce$sample_id)))
ng <- length(gids <- set_names(levels(sce$group_id)))

#Choose color palettes for cluster, sample, group IDs, and # cells
pal <- CATALYST:::.cluster_cols
cluster_id_pal <- set_names(pal[seq_len(nk)], kids)
sample_id_pal <- set_names(pal[seq_len(ns) + nk], sids)
group_id_pal <- set_names(c("royalblue", "orange", "red", "green", "purple"), gids)

#Generate relative cluster abundances
fqs <- prop.table(n_cells, margin = 2)
mat <- as.matrix(unclass(fqs))
Heatmap(mat,
        col = rev(brewer.pal(11, "RdGy")[-6]),
        name = "Frequency",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_side = "left",
        row_title = "cluster_id",
        column_title = "sample_id",
        column_title_side = "bottom",
        rect_gp = gpar(col = "white"),
        cell_fun = function(i, j, x, y, width, height, fill)
          grid.text(round(mat[j, i] * 100, 2), x = x, y = y, 
                    gp = gpar(col = "white", fontsize = 8)))

#DR colored by cluster ID
cs <- sample(colnames(so_astro_merge), 5e3)
.plot_dr <- function(so_astro_merge, dr, id)
  DimPlot(so_astro_merge, cells = cs, group.by = id, reduction = dr, pt.size = 0.4) +
  scale_color_manual(id, values = get(paste0(id, "_pal"))) + 
  guides(col = guide_legend(nrow = 10, 
                            override.aes = list(size = 3, alpha = 1))) +
  theme_void() + theme(aspect.ratio = 1)

ids <- c("cluster_id", "group_id", "sample_id")
for (id in ids) {
  cat("## ", id, "\n")
  p1 <- .plot_dr(so_astro_merge, "tsne", id)
  lgd <- get_legend(p1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- .plot_dr(so_astro_merge, "umap", id) + theme(legend.position = "none")
  ps <- plot_grid(plotlist = list(p1, p2), nrow = 1)
  p <- plot_grid(ps, lgd, nrow = 1, rel_widths = c(1, 0.2))
  print(p)
  cat("\n\n")
}

#QC METRICS CHECK
mito.genes <- grep(pattern = "^mt-", x = rownames(so_astro_merge@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(so_astro_merge@assays[["RNA"]][mito.genes, ])/Matrix::colSums(so_astro_merge@assays[["RNA"]])
so_astro_merge$percent.mito <- percent.mito

#Reorder levels for cluster_id
so_astro_merge_copy <- so_astro_merge
Idents(so_astro_merge_copy) <- "cluster_id"
designated_levels <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
Idents(so_astro_merge_copy) <- factor(Idents(so_astro_merge_copy), levels= designated_levels)
DefaultAssay(so_astro_merge_copy) <- "RNA"
#Visualize QC metric features/counts/mitochondrial contamination across clusters
VlnPlot(object = so_astro_merge_copy, features = c("nFeature_RNA"), ncol = 1, pt.size = 0, cols = c("purple2", "orangered2", "goldenrod2", "royalblue"))
VlnPlot(object = so_astro_merge_copy, features = c("nCount_RNA"), ncol = 1, pt.size = 0, cols = c("purple2", "orangered2", "goldenrod2", "royalblue"))
VlnPlot(object = so_astro_merge_copy, features = c("percent.mito"), ncol = 1, pt.size = 0, cols = c("purple2", "orangered2", "goldenrod2", "royalblue"))

#Reorder levels for dataset
Idents(so_astro_merge_copy) <- "dataset"
designated_levels <- c("Habib", "Zhou_7mo", "Zhou_15mo", "Hasel")
Idents(so_astro_merge_copy) <- factor(Idents(so_astro_merge_copy), levels= designated_levels)
DefaultAssay(so_astro_merge_copy) <- "RNA"
#Visualize QC metric features/counts/mitochondrial contamination across datasets
VlnPlot(object = so_astro_merge_copy, features = c("nFeature_RNA"), ncol = 1, pt.size = 0, cols = c("purple2", "orangered2", "goldenrod2", "royalblue"))
VlnPlot(object = so_astro_merge_copy, features = c("nCount_RNA"), ncol = 1, pt.size = 0, cols = c("purple2", "orangered2", "goldenrod2", "royalblue"))
VlnPlot(object = so_astro_merge_copy, features = c("percent.mito"), ncol = 1, pt.size = 0, cols = c("purple2", "orangered2", "goldenrod2", "royalblue"))

#Generate summary statistics per sample
library(data.table)
library(psych)
feature_by_sample <- as.data.frame(so_astro_merge_copy$nFeature_RNA, row.names = so_astro_merge_copy$sample_id)
feature_by_sample_table <- describeBy(feature_by_sample, group = so_astro_merge_copy$sample_id, mat = TRUE)
write.csv(feature_by_sample_table, "file_path/so_astro_merge_HHZ_feature_by_sample.csv")

count_by_sample <- as.data.frame(so_astro_merge_copy$nCount_RNA, row.names = so_astro_merge_copy$sample_id)
count_by_sample_table <- describeBy(count_by_sample, group = so_astro_merge_copy$sample_id, mat = TRUE)
write.csv(count_by_sample_table, "file_path/so_astro_merge_HHZ_count_by_sample.csv")

#Generate summary statistics per cluster
feature_by_cluster <- as.data.frame(so_astro_merge_copy$nFeature_RNA, row.names = so_astro_merge_copy$cluster_id)
feature_by_cluster_table <- describeBy(feature_by_cluster, group = so_astro_merge_copy$cluster_id, mat = TRUE)
write.csv(feature_by_cluster_table, "file_path/so_astro_merge_HHZ_feature_by_cluster.csv")

count_by_cluster <- as.data.frame(so_astro_merge_copy$nCount_RNA, row.names = so_astro_merge_copy$cluster_id)
count_by_cluster_table <- describeBy(count_by_cluster, group = so_astro_merge_copy$cluster_id, mat = TRUE)
write.csv(count_by_cluster_table, "file_path/so_astro_merge_HHZ_count_by_cluster.csv")

#DETERMINE WHAT DEFINES EACH ASTROCYTE CLUSTER
#Find all markers
DefaultAssay(so_astro_merge_copy) <- "RNA"
so_astro_merge_copy.markers <- FindAllMarkers(so_astro_merge_copy, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so_astro_merge_copy.markers, "file_path/so_astro_merge_HHZ_genes.csv")

#Visualize data (examples)
DefaultAssay(so_astro_merge_copy) <- "RNA"
Idents(so_astro_merge_copy) <- "cluster_id"
DimPlot(so_astro_merge_copy, reduction = "tsne", pt.size = 0.001) + theme(aspect.ratio = 1) + scale_color_manual(values = CATALYST:::.cluster_cols)

Idents(so_astro_merge_copy) <- "dataset"
DimPlot(so_astro_merge_copy, reduction = "tsne", pt.size = 0.001, cols = c("purple2", "orangered2", "goldenrod2", "royalblue")) + theme(aspect.ratio = 1)

#COMPARE WHITE MATTER-SPECIFIC (current study only cluster 4) AND IFN RESPONSE-SPECIFIC (current study only cluster 8) ASTROCYTE FEATURES IN MULTI-DATASET ASTROCYTE ONLY SeuratObject
#Add combined clusterXdataset ident
Idents(so_astro_merge_copy) <- "cluster_id"
so_astro_merge_copy$cluster.dataset_id <- paste(Idents(so_astro_merge_copy), so_astro_merge_copy$dataset, sep = "_")
so_astro_merge_copy$cluster <- Idents(so_astro_merge_copy)
Idents(so_astro_merge_copy) <- "cluster.dataset_id"

#Organize levels of clusterXdonor ident
designated_levels <- c("0_Hasel", "0_Zhou_15mo", "0_Zhou_7mo", "0_Habib", "1_Hasel", "1_Zhou_15mo", "1_Zhou_7mo", "1_Habib", "2_Hasel", "2_Zhou_15mo", "2_Zhou_7mo", "2_Habib", "3_Hasel", "3_Zhou_15mo", "3_Zhou_7mo", "3_Habib", "4_Hasel", "4_Zhou_15mo", "4_Zhou_7mo", "4_Habib", "5_Hasel", "5_Zhou_15mo", "5_Zhou_7mo", "5_Habib", "6_Hasel", "6_Zhou_15mo", "6_Zhou_7mo", "6_Habib", "7_Hasel", "7_Zhou_15mo", "7_Zhou_7mo", "7_Habib", "8_Hasel", "8_Zhou_15mo", "8_Zhou_7mo", "8_Habib", "9_Hasel", "9_Zhou_15mo", "9_Zhou_7mo", "9_Habib", "10_Hasel", "10_Zhou_15mo", "10_Zhou_7mo", "10_Habib", "11_Hasel", "11_Zhou_15mo", "11_Zhou_7mo", "11_Habib")
Idents(so_astro_merge_copy) <- factor(Idents(so_astro_merge_copy), levels= designated_levels)
DefaultAssay(so_astro_merge_copy) <- "RNA"

#Visualize white matter-specific astrocyte top features
DotPlot(so_astro_merge_copy, features = c("Lcn2", "Ifitm3", "H2-D1", "Timp1", "B2m", "Vim", "Gfap", "Psmb8", "Bst2", "H2-T23", "Cd63", "Serping1", "Tspo", "Rps4x", "S100a6", "Tuba1a", "H2-K1", "Mt1", "Gpx1", "S100a11", "Fabp5", "Aqp4", "Eif1b", "Pfn1", "Hspb1", "Gap43", "Cfl1", "Eif1")) + theme(axis.text.x = element_text(angle = 45, hjust=1))

#Visualize IFN response-specific astrocyte top features
DotPlot(so_astro_merge_copy, features = c("Igtp", "Gm4951", "Iigp1", "Ifit3", "Gbp2", "Cxcl10", "Irgm1", "Tap1", "Isg15", "Ifit1", "H2-T23", "Stat1", "F830016B08Rik", "Psmb8", "Tap2", "B2m", "Ifit3b", "Usp18", "Ifitm3", "Gbp7", "Oasl2", "Bst2", "Rtp4", "Rsad2", "H2-D1", "Psmb10", "Serping1", "Stat2", "Gbp3", "H2-Ab1", "H2-K1", "Herc6", "Ifi47", "Irf7", "Lgals3bp", "Rnf213", "Ddx60", "Psmb9", "H2-Q4", "Xaf1", "Lcn2", "Irf1", "Trim30a", "Tapbp", "Tapbpl")) + theme(axis.text.x = element_text(angle = 45, hjust=1))

#SUBSET OUT WHITE MATTER-SPECIFIC OR IFN RESPONSE-SPECIFIC ASTROCYTE CLUSTERS IN MULTI-DATASET ASTROCYTE ONLY SeuratObject
#Subset only white matter-specific astrocyte cluster (i.e., multi-dataset cluster 7)
Idents(so_astro_merge_copy) <- "cluster_id"
so_astro_merge_WM <- subset(x = so_astro_merge_copy, idents = c("7"), invert = FALSE)

#Add combined datasetXgroup ident
Idents(so_astro_merge_WM) <- "group_id"
so_astro_merge_WM$group.dataset_id <- paste(Idents(so_astro_merge_WM), so_astro_merge_WM$dataset, sep = "_")
so_astro_merge_WM$group <- Idents(so_astro_merge_WM)
Idents(so_astro_merge_WM) <- "group.dataset_id"

#Organize levels of groupXdonor ident
designated_levels <- c("CNT_Habib", "5XFAD_Habib", "CNT_Zhou_7mo", "5XFAD_Zhou_7mo", "TREM2_KO_Zhou_7mo", "TREM2_KO_5XFAD_Zhou_7mo", "CNT_Zhou_15mo", "5XFAD_Zhou_15mo", "CNT_Hasel", "LPS_Hasel")
Idents(so_astro_merge_WM) <- factor(Idents(so_astro_merge_WM), levels= designated_levels)
DefaultAssay(so_astro_merge_WM) <- "RNA"

#Plot white matter-specific astrocyte transcripts of interest
VlnPlot(so_astro_merge_WM, features = c("Timp1"), pt.size = 0)

#Subset only IFN response-specific astrocyte cluster (i.e., multi-dataset cluster 9)
Idents(so_astro_merge_copy) <- "cluster_id"
so_astro_merge_doomer <- subset(x = so_astro_merge_copy, idents = c("9"), invert = FALSE)

#Add combined datasetXgroup ident
Idents(so_astro_merge_doomer) <- "group_id"
so_astro_merge_doomer$group.dataset_id <- paste(Idents(so_astro_merge_doomer), so_astro_merge_doomer$dataset, sep = "_")
so_astro_merge_doomer$group <- Idents(so_astro_merge_doomer)
Idents(so_astro_merge_doomer) <- "group.dataset_id"

#Organize levels of groupXdonor ident
designated_levels <- c("CNT_Habib", "5XFAD_Habib", "CNT_Zhou_7mo", "5XFAD_Zhou_7mo", "TREM2_KO_Zhou_7mo", "TREM2_KO_5XFAD_Zhou_7mo", "CNT_Zhou_15mo", "5XFAD_Zhou_15mo", "CNT_Hasel", "LPS_Hasel")
Idents(so_astro_merge_doomer) <- factor(Idents(so_astro_merge_doomer), levels= designated_levels)
DefaultAssay(so_astro_merge_doomer) <- "RNA"

#Plot IFN response-specific astrocyte transcripts of interest
VlnPlot(so_astro_merge_doomer, features = c("Igtp"), pt.size = 0)
