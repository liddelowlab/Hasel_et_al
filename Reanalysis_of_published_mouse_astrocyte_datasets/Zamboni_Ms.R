#Pipeline prepared and Zamboni et al. (2020) data reanalyzed by Jessica S. Sadick
#Majority of pipeline based on code originally deposited by https://github.com/HelenaLC and published in doi.org/10.1101/713412

#---------------------------------------------------------------------------------------------------
#DEMULTIPLEX AGGREGATED FILES
##NOTE: Used filtered dataset, which contains astrocytes and progenitor cells

#Load libraries
library(data.table)
library(Matrix)
library(dplyr)

setwd("file_path")

#Read in aggregated matrix file
mat <- fread("zcat < GSE139842_filtered_10X_gbm.csv.gz")

#Isolate barcodes for each sample, and write to a tsv file
barcodes_healthy_wt <- mat[, grep("healthy_wt_", colnames(mat[1]), value = TRUE)]
write.table(barcodes_healthy_wt, file='./barcodes_healthy_wt.tsv', quote=FALSE,
            sep='\t', col.names = FALSE, row.names = FALSE)

barcodes_inj_wt <- mat[, grep("inj_wt_", colnames(mat[1]), value = TRUE)]
write.table(barcodes_inj_wt, file='./barcodes_inj_wt.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcodes_healthy_hom <- mat[, grep("healthy_hom_", colnames(mat[1]), value = TRUE)]
write.table(barcodes_healthy_hom, file='./barcodes_healthy_hom.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcodes_inj_hom <- mat[, grep("inj_hom_", colnames(mat[1]), value = TRUE)]
write.table(barcodes_inj_hom, file='./barcodes_inj_hom.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

barcodes_inj_hom_2 <- mat[, grep('^[A-Z]', colnames(mat[1]), value = TRUE)]
barcodes_inj_hom_2 <- barcodes_inj_hom_2[2:2642]
write.table(barcodes_inj_hom_2, file='./barcodes_inj_hom_2.tsv', quote=FALSE, 
            sep='\t', col.names = FALSE, row.names = FALSE)

#Isolate gene information from original matrix
genes <- mat[,1][[1]]
#Isolate gene name
gene_name <- gsub(".+[|]", "", genes)
gene_name <- as.data.frame(gene_name)
#Add gene ENSEMBL ID
library(biomaRt)
ensembl <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
gene_id <- getBM(attributes=c("ensembl_gene_id"),
                 filters = c("mgi_symbol"),
                 values = gene_name,
                 mart = ensembl)
gene_id <- as.data.frame(gene_id)
##Concatenate gene ID and gene name in one dataframe
##NOTE: For downstream processing, feature file must be organized in a 3 column format. Can concatenate gene_name twice as a place holder if EMSEMBL ID conversion is not ideal
features <- cbind(gene_id,gene_name)
##Add "Gene Expression" column
features["new.col"] <- c("Gene Expression")
##Write to tsv file
write.table(features, file='./features.tsv', quote=FALSE,
            sep='\t', col.names = FALSE, row.names = FALSE)

#Subset matrix by sample
mat_healthy_wt <- select(mat, contains("healthy_wt_"))
mat_inj_wt <- select(mat, contains("inj_wt_"))
mat_healthy_hom <- select(mat, contains("healthy_hom_"))
mat_inj_hom <- select(mat, contains("inj_hom_"))
mat_inj_hom_2 <- select(mat, -contains(c("healthy_wt_", "inj_wt_", "healthy_hom_", "inj_hom_", "V1")))

#Write each new subset to new matrix file
##Remove column and row names
names(mat_healthy_wt) = NULL
##Convert dataframe to sparse matrix
mat_healthy_wt <- Matrix(as.matrix(mat_healthy_wt), sparse=TRUE)
##Write sparse matrix to file
writeMM(obj = mat_healthy_wt, file="./mat_healthy_wt.mtx")

names(mat_inj_wt) = NULL
mat_inj_wt <- Matrix(as.matrix(mat_inj_wt), sparse=TRUE)
writeMM(obj = mat_inj_wt, file="./mat_inj_wt.mtx")

names(mat_healthy_hom) = NULL
mat_healthy_hom <- Matrix(as.matrix(mat_healthy_hom), sparse=TRUE)
writeMM(obj = mat_healthy_hom, file="./mat_healthy_hom.mtx")

names(mat_inj_hom) = NULL
mat_inj_hom <- Matrix(as.matrix(mat_inj_hom), sparse=TRUE)
writeMM(obj = mat_inj_hom, file="./mat_inj_hom.mtx")

names(mat_inj_hom_2) = NULL
mat_inj_hom_2 <- Matrix(as.matrix(mat_inj_hom_2), sparse=TRUE)
writeMM(obj = mat_inj_hom_2, file="./mat_inj_hom_2.mtx")

#---------------------------------------------------------------------------------------------------
#SCE OBJECT GENERATION AND QC

#Load libraries
library(cowplot)
library(ggplot2)
library(scater)
library(scds)
library(SingleCellExperiment)

#LOAD AND REFORMAT DATA
#Load raw counts
fastq_dirs <- list.dirs("file_path", recursive = FALSE, full.names = TRUE)
names(fastq_dirs) <- basename(fastq_dirs)
sce <- DropletUtils::read10xCounts(fastq_dirs)

#Rename row/colData colnames and SCE dimnames
names(rowData(sce)) <- c("ENSEMBL", "SYMBOL")
names(colData(sce)) <- c("sample_id", "barcode")
sce$sample_id <- factor(basename(sce$sample_id))
dimnames(sce) <- list(
  with(rowData(sce), paste(SYMBOL, sep = ".")),
  with(colData(sce), paste(barcode, sample_id, sep = ".")))

#Load metadata
md_dir <- file.path("file_path", "metadata_Zamboni_Ms.xlsx")
md <- readxl::read_excel(md_dir)
m <- match(sce$sample_id, md$`Sample ID`)
sce$group_id <- md$Characteristics[m]
sce$sex_id <- md$Sex[m]
sce$age_id <- md$Age[m]
sce$split_id <- md$Split[m]

#Remove undetected genes
sce <- sce[Matrix::rowSums(counts(sce) > 0) > 0, ]
dim(sce)

#DOUBLET REMOVAL
#Split SCE by sample
cs_by_s <- split(colnames(sce), sce$sample_id)
sce_by_s <- lapply(cs_by_s, function(cs) sce[, cs])

#Run 'scds'
sce_by_s <- lapply(sce_by_s, function(u) 
  cxds_bcds_hybrid(bcds(cxds(u))))

#Remove doublets
sce_by_s <- lapply(sce_by_s, function(u) {
  #Compute expected nb. of doublets (10x)
  n_dbl <- ceiling(0.01 * ncol(u)^2 / 1e3)
  #Remove 'n_dbl' cells w/ highest doublet score
  o <- order(u$hybrid_score, decreasing = TRUE)
  u[, -o[seq_len(n_dbl)]]
})

#Merge back into single SCE
sce <- do.call(cbind, sce_by_s)

#CALCULATE QC METRICS
(mito <- grep("mt-", rownames(sce), value = TRUE))
sce <- calculateQCMetrics(sce, feature_controls = list(Mt = mito))
plotHighestExprs(sce, n = 20)

#FILTERING
#Get sample-specific outliers
cols <- c("total_counts", "total_features_by_counts", "pct_counts_Mt")
log <- c(TRUE, TRUE, FALSE)
type <- c("both", "both", "higher")

drop_cols <- paste0(cols, "_drop")
for (i in seq_along(cols))
  colData(sce)[[drop_cols[i]]] <- isOutlier(sce[[cols[i]]], 
                                            nmads = 2.5, type = type[i], log = log[i], batch = sce$sample_id)

sapply(drop_cols, function(i) 
  sapply(drop_cols, function(j)
    sum(sce[[i]] & sce[[j]])))

cd <- data.frame(colData(sce))
ps <- lapply(seq_along(cols), function (i) {
  p <- ggplot(cd, aes_string(x = cols[i], alpha = drop_cols[i])) +
    geom_histogram(bins = 100, show.legend = FALSE) +
    scale_alpha_manual(values = c("FALSE" = 1, "TRUE" = 0.4)) +
    facet_wrap(~sample_id, ncol = 1, scales = "free") + 
    theme_classic() + theme(strip.background = element_blank())
  if (log[i]) 
    p <- p + scale_x_log10()
  return(p)
})
plot_grid(plotlist = ps, ncol = 3)

layout(matrix(1:2, nrow = 1))
ol <- Matrix::rowSums(as.matrix(colData(sce)[drop_cols])) != 0
x <- sce$total_counts
y <- sce$total_features_by_counts
LSD::heatscatter(x, y, log="xy", main = "unfiltered", 
                 xlab = "Total counts", ylab = "Non-zero features")
LSD::heatscatter(x[!ol], y[!ol], log="xy", main = "filtered", 
                 xlab = "Total counts", ylab = "Non-zero features")

#Summary of cells kept
ns <- table(sce$sample_id)
ns_fil <- table(sce$sample_id[!ol])
print(rbind(
  unfiltered = ns, filtered = ns_fil, 
  "%" = ns_fil / ns * 100), digits = 0)

#Drop outlier cells
sce <- sce[, !ol]
dim(sce)

#Require count > 1 in at least 20 cells
sce <- sce[Matrix::rowSums(counts(sce) > 1) >= 20, ]
dim(sce)

#SAVE SCE
saveRDS(sce, file.path("file_path", "Zamboni_Ms_SCE.rds"))

#---------------------------------------------------------------------------------------------------
#CLUSTERING

#Increase future's maximum allowed size of objects
options(future.globals.maxSize = 2048 * 1024 ^20)
memory.limit(size = 1000000)

#Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

#Load SCE
sce <- readRDS(file.path("file_path", "Zamboni_Ms_SCE.rds"))

#INTEGRATE
#Create SeuratObject
so_astro <- CreateSeuratObject(
  counts = counts(sce),
  meta.data = data.frame(colData(sce)),
  project = "Zamboni_Ms_10x_data")

#Split by sample
cells_by_sample <- split(colnames(sce), sce$sample_id)
so_astro <- lapply(cells_by_sample, function(i)
  SubsetData(so_astro, cells = i))

#Normalize, find variable genes, and scale
so_astro <- lapply(so_astro, NormalizeData, verbose = FALSE)
so_astro <- lapply(so_astro, FindVariableFeatures, nfeatures = 2e3,
             selection.method = "vst", do.plot = FALSE, verbose = FALSE)
so_astro <- lapply(so_astro, ScaleData, verbose = FALSE)

#Find anchors and integrate
as <- FindIntegrationAnchors(so_astro, verbose = FALSE)
so_astro <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = FALSE)

#Scale integrated data
DefaultAssay(so_astro) <- "integrated"
so_astro <- ScaleData(so_astro, display.progress = FALSE)

#DIMENSION REDUCTION 
so_astro <- RunPCA(so_astro, npcs = 100, verbose = FALSE)
ElbowPlot(so_astro, ndims = 100)
##Update number of PCs used
so_astro <- RunTSNE(so_astro, reduction = "pca", dims = seq_len(15),
              seed.use = 1, do.fast = TRUE, verbose = FALSE)
so_astro <- RunUMAP(so_astro, reduction = "pca", dims = seq_len(15),
              seed.use = 1, verbose = FALSE)

#CLUSTERING 
so_astro <- FindNeighbors(so_astro, reduction = "pca", dims = seq_len(15), verbose = FALSE)
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
  so_astro <- FindClusters(so_astro, resolution = res, random.seed = 1, verbose = FALSE)

#DR COLORED BY SAMPLE, GROUP, AND CLUSTER ID
thm <- theme(aspect.ratio = 1, legend.position = "none")
ps <- lapply(c("sample_id", "group_id", "ident"), function(u) {
  p1 <- DimPlot(so_astro, reduction = "tsne", group.by = u) + thm
  p2 <- DimPlot(so_astro, reduction = "umap", group.by = u)
  lgd <- get_legend(p2)
  p2 <- p2 + thm
  list(p1, p2, lgd)
  plot_grid(p1, p2, lgd, nrow = 1,
            rel_widths = c(1, 1, 0.5))
})
plot_grid(plotlist = ps, ncol = 1)

#SAVE SeuratObject
saveRDS(so_astro, file.path("file_path", "Zamboni_so_astro_PC15.rds"))

#---------------------------------------------------------------------------------------------------
#CLUSTER ANNOTATION

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
so_astro <- readRDS(file.path("file_path", "Zamboni_so_astro_PC15.rds"))
sce <- as.SingleCellExperiment(so_astro, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>% 
  mutate_if(is.character, as.factor) %>% 
  DataFrame(row.names = colnames(sce))

#Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)

so_astro <- SetIdent(so_astro, value = "integrated_snn_res.0.1")
so_astro@meta.data$cluster_id <- Idents(so_astro)
sce$cluster_id <- Idents(so_astro)
(n_cells <- table(sce$cluster_id, sce$sample_id))
write.csv(table(sce$cluster_id, sce$sample_id), "file_path/Zamboni_so_astro_numbers.csv")

nk <- length(kids <- set_names(levels(sce$cluster_id)))
ns <- length(sids <- set_names(levels(sce$sample_id)))
ng <- length(gids <- set_names(levels(sce$group_id)))

#Choose color palettes for cluster, sample, group IDs, and # cells
pal <- CATALYST:::.cluster_cols
cluster_id_pal <- set_names(pal[seq_len(nk)], kids)
sample_id_pal <- set_names(pal[seq_len(ns) + nk], sids)
group_id_pal <- set_names(c("royalblue", "orange", "red", "green"), gids)

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
cs <- sample(colnames(so_astro), 5e3)
.plot_dr <- function(so_astro, dr, id)
  DimPlot(so_astro, cells = cs, group.by = id, reduction = dr, pt.size = 0.4) +
  scale_color_manual(id, values = get(paste0(id, "_pal"))) + 
  guides(col = guide_legend(nrow = 10, 
                            override.aes = list(size = 3, alpha = 1))) +
  theme_void() + theme(aspect.ratio = 1)

ids <- c("cluster_id", "group_id", "sample_id")
for (id in ids) {
  cat("## ", id, "\n")
  p1 <- .plot_dr(so_astro, "tsne", id)
  lgd <- get_legend(p1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- .plot_dr(so_astro, "umap", id) + theme(legend.position = "none")
  ps <- plot_grid(plotlist = list(p1, p2), nrow = 1)
  p <- plot_grid(ps, lgd, nrow = 1, rel_widths = c(1, 0.2))
  print(p)
  cat("\n\n")
}

#QC METRICS CHECK
mito.genes <- grep(pattern = "^mt-", x = rownames(so_astro@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(so_astro@assays[["RNA"]][mito.genes, ])/Matrix::colSums(so_astro@assays[["RNA"]])
so_astro$percent.mito <- percent.mito

VlnPlot(object = so_astro, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

#Generate summary statistics per sample
library(data.table)
library(psych)
feature_by_sample <- as.data.frame(so_astro$nFeature_RNA, row.names = so_astro$sample_id)
feature_by_sample_table <- describeBy(feature_by_sample, group = so_astro$sample_id, mat = TRUE)
write.csv(feature_by_sample_table, "file_path/Zamboni_so_asto_feature_by_sample.csv")

count_by_sample <- as.data.frame(so_astro$nCount_RNA, row.names = so_astro$sample_id)
count_by_sample_table <- describeBy(count_by_sample, group = so_astro$sample_id, mat = TRUE)
write.csv(count_by_sample_table, "file_path/Zamboni_so_asto_count_by_sample.csv")

#summary statistics per cluster
feature_by_cluster <- as.data.frame(so_astro$nFeature_RNA, row.names = so_astro$cluster_id)
feature_by_cluster_table <- describeBy(feature_by_cluster, group = so_astro$cluster_id, mat = TRUE)
write.csv(feature_by_cluster_table, "file_path/Zamboni_so_asto_feature_by_cluster.csv")

count_by_cluster <- as.data.frame(so_astro$nCount_RNA, row.names = so_astro$cluster_id)
count_by_cluster_table <- describeBy(count_by_cluster, group = so_astro$cluster_id, mat = TRUE)
write.csv(count_by_cluster_table, "file_path/Zamboni_so_asto_count_by_cluster.csv")

#DETERMINE WHAT DEFINES EACH ASTROCYTE CLUSTER
#Find all markers
DefaultAssay(so_astro) <- "RNA"
so_astro.markers <- FindAllMarkers(so_astro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so_astro.markers, "file_path/Zamboni_so_astro_genes.csv")

#Visualize data (examples)
DefaultAssay(so_astro) <- "RNA"
DimPlot(so_astro, reduction = "tsne") + theme(aspect.ratio = 1) + scale_color_manual(values = CATALYST:::.cluster_cols)

VlnPlot(so_astro, features = c("Aqp4", "Gfap", "Fgfr3", "Cldn10", "Gja1", "Aldh1l1", "Slc1a3", "Slc1a2"), pt.size = 0)

FeaturePlot(so_astro, features = c("Igtp"), reduction = "tsne", pt.size = 0.001) + theme(aspect.ratio = 1)

#---------------------------------------------------------------------------------------------------
#MAKE FINAL ASTROCYTE-SPECIFIC SeuratObject
##Original astrocyte object had neural stem/progenitor cell contamination
#Create SeuratObject with only astrocyte clusters
so_astro_final <- subset(x = so_astro, idents = c("2"), invert = TRUE)

#RECLUSTER FINAL ASTROCYTE-SPECIFIC SeuratObject

#Increase future's maximum allowed size of objects
options(future.globals.maxSize = 2048 * 1024 ^20)
memory.limit(size = 1000000)

#Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

#Convert SeuratObject to SCE
sce <- as.SingleCellExperiment(so_astro_final, assay = "RNA")

#INTEGRATE
#Split by sample
cells_by_sample <- split(colnames(sce), sce$sample_id)
so_astro_final <- lapply(cells_by_sample, function(i)
  SubsetData(so_astro_final, cells = i))

#Normalize, find variable genes, and scale
so_astro_final <- lapply(so_astro_final, NormalizeData, verbose = FALSE)
so_astro_final <- lapply(so_astro_final, FindVariableFeatures, nfeatures = 2e3,
                         selection.method = "vst", do.plot = FALSE, verbose = FALSE)
so_astro_final <- lapply(so_astro_final, ScaleData, verbose = FALSE)

#Find anchors and integrate
##Decrease k.filter to minimize number of astrocytes identified per sample
as <- FindIntegrationAnchors(so_astro_final, k.filter = 185, verbose = FALSE)
so_astro_final <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = FALSE)

#Scale integrated data
DefaultAssay(so_astro_final) <- "integrated"
so_astro_final <- ScaleData(so_astro_final, display.progress = FALSE)

#DIMENSION REDUCTION
so_astro_final <- RunPCA(so_astro_final, npcs = 50, verbose = FALSE)
ElbowPlot(so_astro_final, ndims = 50)
##Update number of PCs used
so_astro_final <- RunTSNE(so_astro_final, reduction = "pca", dims = seq_len(12),
                          seed.use = 1, do.fast = TRUE, verbose = FALSE)
so_astro_final <- RunUMAP(so_astro_final, reduction = "pca", dims = seq_len(12),
                          seed.use = 1, verbose = FALSE)

#CLUSTERING
so_astro_final <- FindNeighbors(so_astro_final, reduction = "pca", dims = seq_len(12), verbose = FALSE)
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
  so_astro_final <- FindClusters(so_astro_final, resolution = res, random.seed = 1, verbose = FALSE)

#DR COLORED BY SAMPLE, GROUP, AND CLUSTER ID
thm <- theme(aspect.ratio = 1, legend.position = "none")
ps <- lapply(c("sample_id", "group_id", "ident"), function(u) {
  p1 <- DimPlot(so_astro_final, reduction = "tsne", group.by = u) + thm
  p2 <- DimPlot(so_astro_final, reduction = "umap", group.by = u)
  lgd <- get_legend(p2)
  p2 <- p2 + thm
  list(p1, p2, lgd)
  plot_grid(p1, p2, lgd, nrow = 1,
            rel_widths = c(1, 1, 0.5))
})
plot_grid(plotlist = ps, ncol = 1)

#SAVE SeuratObject
saveRDS(so_astro_final, file.path("file_path", "Zamboni_Ms_so_astro_final_12PC.rds"))

#---------------------------------------------------------------------------------------------------
#CLUSTER ANNOTATION FOR FINAL ASTROCYTE ONLY SeuratObject

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
so_astro_final <- readRDS(file.path("file_path", "Zamboni_Ms_so_astro_final_12PC.rds"))
sce <- as.SingleCellExperiment(so_astro_final, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>% 
  mutate_if(is.character, as.factor) %>% 
  DataFrame(row.names = colnames(sce))

#Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)

so_astro_final <- SetIdent(so_astro_final, value = "integrated_snn_res.0.2")
so_astro_final@meta.data$cluster_id <- Idents(so_astro_final)
sce$cluster_id <- Idents(so_astro_final)
(n_cells <- table(sce$cluster_id, sce$sample_id))
write.csv(table(sce$cluster_id, sce$sample_id), "file_path/Zamboni_so_astro_final_numbers.csv")

nk <- length(kids <- set_names(levels(sce$cluster_id)))
ns <- length(sids <- set_names(levels(sce$sample_id)))
ng <- length(gids <- set_names(levels(sce$group_id)))

#Choose color palettes for cluster, sample, group IDs, and # cells
pal <- CATALYST:::.cluster_cols
cluster_id_pal <- set_names(pal[seq_len(nk)], kids)
sample_id_pal <- set_names(pal[seq_len(ns) + nk], sids)
group_id_pal <- set_names(c("royalblue", "orange", "red", "green"), gids)

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

#DR colored by cluster ID #changed 5e3 to 5e2
cs <- sample(colnames(so_astro_final), 5e2)
.plot_dr <- function(so_astro_final, dr, id)
  DimPlot(so_astro_final, cells = cs, group.by = id, reduction = dr, pt.size = 0.4) +
  scale_color_manual(id, values = get(paste0(id, "_pal"))) + 
  guides(col = guide_legend(nrow = 10, 
                            override.aes = list(size = 3, alpha = 1))) +
  theme_void() + theme(aspect.ratio = 1)

ids <- c("cluster_id", "group_id", "sample_id")
for (id in ids) {
  cat("## ", id, "\n")
  p1 <- .plot_dr(so_astro_final, "tsne", id)
  lgd <- get_legend(p1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- .plot_dr(so_astro_final, "umap", id) + theme(legend.position = "none")
  ps <- plot_grid(plotlist = list(p1, p2), nrow = 1)
  p <- plot_grid(ps, lgd, nrow = 1, rel_widths = c(1, 0.2))
  print(p)
  cat("\n\n")
}

#QC METRICS CHECK
mito.genes <- grep(pattern = "^mt-", x = rownames(so_astro_final@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(so_astro_final@assays[["RNA"]][mito.genes, ])/Matrix::colSums(so_astro_final@assays[["RNA"]])
so_astro_final$percent.mito <- percent.mito

VlnPlot(object = so_astro_final, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

#DETERMINE WHAT DEFINES EACH ASTROCYTE CLUSTER
#Find all markers
DefaultAssay(so_astro_final) <- "RNA"
so_astro_final.markers <- FindAllMarkers(so_astro_final, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so_astro_final.markers, "file_path/Zamboni_so_astro_final_genes.csv")

#Visualize data (examples)
DefaultAssay(so_astro_final) <- "RNA"
DimPlot(so_astro_final, reduction = "tsne", pt.size = 0.001) + theme(aspect.ratio = 1) + scale_color_manual(values = CATALYST:::.cluster_cols)

VlnPlot(so_astro_final, features = c("Aqp4", "Gfap", "Fgfr3", "Cldn10", "Gja1", "Aldh1l1", "Slc1a3", "Slc1a2"), pt.size = 0)

FeaturePlot(so_astro_final, features = c("Igtp"), reduction = "tsne", pt.size = 0.001) + theme(aspect.ratio = 1)
