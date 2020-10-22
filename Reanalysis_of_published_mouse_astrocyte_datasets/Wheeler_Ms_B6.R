#Load relevant libraries
library(cowplot)
library(ggplot2)
library(scater)
library(scds)
library(SingleCellExperiment)
library(data.table)
library(Matrix)
library(dplyr)
library(Seurat)

#LOAD AND REFORMAT DATA
#Load DGE files for each sample
setwd("file_path")
mat_GSM3717013 <- read.table(file = "GSM3717013_1_HG37GDMXX.1.TAAGGCGA.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM3717014 <- read.table(file = "GSM3717014_1_HG37GDMXX.1.CGTACTAG.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM3717015 <- read.table(file = "GSM3717015_1_HG37GDMXX.1.AGGCAGAA.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM3717016 <- read.table(file = "GSM3717016_1_HG37GDMXX.1.TCCTGAGC.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM3717017 <- read.table(file = "GSM3717017_1_HG37GDMXX.1.GGACTCCT.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM3717018 <- read.table(file = "GSM3717018_1_HG37GDMXX.1.TAGGCATG.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM3717019 <- read.table(file = "GSM3717019_1_HG37GDMXX.1.CTCTCTAC.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM3717020 <- read.table(file = "GSM3717020_1_HG37GDMXX.1.CAGAGAGG.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM3717021 <- read.table(file = "GSM3717021_1_HG37GDMXX.1.GCTACGCT.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM3717022 <- read.table(file = "GSM3717022_1_HG37GDMXX.1.CGAGGCTG.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM3717023 <- read.table(file = "GSM3717023_2_HG37GDMXX.2.AAGAGGCA.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM3717024 <- read.table(file = "GSM3717024_2_HG37GDMXX.2.GTAGAGGA.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM3717025 <- read.table(file = "GSM3717025_2_HG37GDMXX.2.GCTCATGA.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM3717026 <- read.table(file = "GSM3717026_2_HG37GDMXX.2.ATCTCAGG.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM3717027 <- read.table(file = "GSM3717027_2_HG37GDMXX.2.ACTCGCTA.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM4002707 <- read.table(file = "GSM4002707_2_H7NF2DRXX.2.TAAGGCGA.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM4002708 <- read.table(file = "GSM4002708_2_H7NF2DRXX.2.CGTACTAG.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM4002709 <- read.table(file = "GSM4002709_2_H7NF2DRXX.2.AGGCAGAA.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM4002710 <- read.table(file = "GSM4002710_2_H7NF2DRXX.2.TCCTGAGC.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM4002711 <- read.table(file = "GSM4002711_2_H7NF2DRXX.2.GCTACGCT.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM4002712 <- read.table(file = "GSM4002712_2_H7NF2DRXX.2.CGAGGCTG.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM4002713 <- read.table(file = "GSM4002713_2_H7NF2DRXX.2.AAGAGGCA.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM4002714 <- read.table(file = "GSM4002714_1_H7NF2DRXX.1.GTAGAGGA.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM4002715 <- read.table(file = "GSM4002715_1_H7NF2DRXX.1.GCTCATGA.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM4002716 <- read.table(file = "GSM4002716_1_H7NF2DRXX.1.ATCTCAGG.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM4002717 <- read.table(file = "GSM4002717_1_H7NF2DRXX.1.ACTCGCTA.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))
mat_GSM4002718 <- read.table(file = "GSM4002718_1_H7NF2DRXX.1.GGAGCTAC.dge.txt", header = TRUE, row.names = 1, colClasses =c("character"))

#Convert each DGE to a SeuratObject
GSM3717013 <- CreateSeuratObject(counts = mat_GSM3717013, project = "GSM3717013")
Meta_GSM3717013 <- GSM3717013@meta.data
Meta_GSM3717013["sample_id"] <- c("GSM3717013")
NewMeta_GSM3717013 <- subset(Meta_GSM3717013, select = c("sample_id"))
GSM3717013 <- AddMetaData(GSM3717013, NewMeta_GSM3717013)
head(x = GSM3717013[[]])

GSM3717014 <- CreateSeuratObject(counts = mat_GSM3717014, project = "GSM3717014")
Meta_GSM3717014 <- GSM3717014@meta.data
Meta_GSM3717014["sample_id"] <- c("GSM3717014")
NewMeta_GSM3717014 <- subset(Meta_GSM3717014, select = c("sample_id"))
GSM3717014 <- AddMetaData(GSM3717014, NewMeta_GSM3717014)
head(x = GSM3717014[[]])

GSM3717015 <- CreateSeuratObject(counts = mat_GSM3717015, project = "GSM3717015")
Meta_GSM3717015 <- GSM3717015@meta.data
Meta_GSM3717015["sample_id"] <- c("GSM3717015")
NewMeta_GSM3717015 <- subset(Meta_GSM3717015, select = c("sample_id"))
GSM3717015 <- AddMetaData(GSM3717015, NewMeta_GSM3717015)
head(x = GSM3717015[[]])

GSM3717016 <- CreateSeuratObject(counts = mat_GSM3717016, project = "GSM3717016")
Meta_GSM3717016 <- GSM3717016@meta.data
Meta_GSM3717016["sample_id"] <- c("GSM3717016")
NewMeta_GSM3717016 <- subset(Meta_GSM3717016, select = c("sample_id"))
GSM3717016 <- AddMetaData(GSM3717016, NewMeta_GSM3717016)
head(x = GSM3717016[[]])

GSM3717017 <- CreateSeuratObject(counts = mat_GSM3717017, project = "GSM3717017")
Meta_GSM3717017 <- GSM3717017@meta.data
Meta_GSM3717017["sample_id"] <- c("GSM3717017")
NewMeta_GSM3717017 <- subset(Meta_GSM3717017, select = c("sample_id"))
GSM3717017 <- AddMetaData(GSM3717017, NewMeta_GSM3717017)
head(x = GSM3717017[[]])

GSM3717018 <- CreateSeuratObject(counts = mat_GSM3717018, project = "GSM3717018")
Meta_GSM3717018 <- GSM3717018@meta.data
Meta_GSM3717018["sample_id"] <- c("GSM3717018")
NewMeta_GSM3717018 <- subset(Meta_GSM3717018, select = c("sample_id"))
GSM3717018 <- AddMetaData(GSM3717018, NewMeta_GSM3717018)
head(x = GSM3717018[[]])

GSM3717019 <- CreateSeuratObject(counts = mat_GSM3717019, project = "GSM3717019")
Meta_GSM3717019 <- GSM3717019@meta.data
Meta_GSM3717019["sample_id"] <- c("GSM3717019")
NewMeta_GSM3717019 <- subset(Meta_GSM3717019, select = c("sample_id"))
GSM3717019 <- AddMetaData(GSM3717019, NewMeta_GSM3717019)
head(x = GSM3717019[[]])

GSM3717020 <- CreateSeuratObject(counts = mat_GSM3717020, project = "GSM3717020")
Meta_GSM3717020 <- GSM3717020@meta.data
Meta_GSM3717020["sample_id"] <- c("GSM3717020")
NewMeta_GSM3717020 <- subset(Meta_GSM3717020, select = c("sample_id"))
GSM3717020 <- AddMetaData(GSM3717020, NewMeta_GSM3717020)
head(x = GSM3717020[[]])

GSM3717021 <- CreateSeuratObject(counts = mat_GSM3717021, project = "GSM3717021")
Meta_GSM3717021 <- GSM3717021@meta.data
Meta_GSM3717021["sample_id"] <- c("GSM3717021")
NewMeta_GSM3717021 <- subset(Meta_GSM3717021, select = c("sample_id"))
GSM3717021 <- AddMetaData(GSM3717021, NewMeta_GSM3717021)
head(x = GSM3717021[[]])

GSM3717022 <- CreateSeuratObject(counts = mat_GSM3717022, project = "GSM3717022")
Meta_GSM3717022 <- GSM3717022@meta.data
Meta_GSM3717022["sample_id"] <- c("GSM3717022")
NewMeta_GSM3717022 <- subset(Meta_GSM3717022, select = c("sample_id"))
GSM3717022 <- AddMetaData(GSM3717022, NewMeta_GSM3717022)
head(x = GSM3717022[[]])

GSM3717023 <- CreateSeuratObject(counts = mat_GSM3717023, project = "GSM3717023")
Meta_GSM3717023 <- GSM3717023@meta.data
Meta_GSM3717023["sample_id"] <- c("GSM3717023")
NewMeta_GSM3717023 <- subset(Meta_GSM3717023, select = c("sample_id"))
GSM3717023 <- AddMetaData(GSM3717023, NewMeta_GSM3717023)
head(x = GSM3717023[[]])

GSM3717024 <- CreateSeuratObject(counts = mat_GSM3717024, project = "GSM3717024")
Meta_GSM3717024 <- GSM3717024@meta.data
Meta_GSM3717024["sample_id"] <- c("GSM3717024")
NewMeta_GSM3717024 <- subset(Meta_GSM3717024, select = c("sample_id"))
GSM3717024 <- AddMetaData(GSM3717024, NewMeta_GSM3717024)
head(x = GSM3717024[[]])

GSM3717025 <- CreateSeuratObject(counts = mat_GSM3717025, project = "GSM3717025")
Meta_GSM3717025 <- GSM3717025@meta.data
Meta_GSM3717025["sample_id"] <- c("GSM3717025")
NewMeta_GSM3717025 <- subset(Meta_GSM3717025, select = c("sample_id"))
GSM3717025 <- AddMetaData(GSM3717025, NewMeta_GSM3717025)
head(x = GSM3717025[[]])

GSM3717026 <- CreateSeuratObject(counts = mat_GSM3717026, project = "GSM3717026")
Meta_GSM3717026 <- GSM3717026@meta.data
Meta_GSM3717026["sample_id"] <- c("GSM3717026")
NewMeta_GSM3717026 <- subset(Meta_GSM3717026, select = c("sample_id"))
GSM3717026 <- AddMetaData(GSM3717026, NewMeta_GSM3717026)
head(x = GSM3717026[[]])

GSM3717027 <- CreateSeuratObject(counts = mat_GSM3717027, project = "GSM3717027")
Meta_GSM3717027 <- GSM3717027@meta.data
Meta_GSM3717027["sample_id"] <- c("GSM3717027")
NewMeta_GSM3717027 <- subset(Meta_GSM3717027, select = c("sample_id"))
GSM3717027 <- AddMetaData(GSM3717027, NewMeta_GSM3717027)
head(x = GSM3717027[[]])

GSM4002707 <- CreateSeuratObject(counts = mat_GSM4002707, project = "GSM4002707")
Meta_GSM4002707 <- GSM4002707@meta.data
Meta_GSM4002707["sample_id"] <- c("GSM4002707")
NewMeta_GSM4002707 <- subset(Meta_GSM4002707, select = c("sample_id"))
GSM4002707 <- AddMetaData(GSM4002707, NewMeta_GSM4002707)
head(x = GSM4002707[[]])

GSM4002708 <- CreateSeuratObject(counts = mat_GSM4002708, project = "GSM4002708")
Meta_GSM4002708 <- GSM4002708@meta.data
Meta_GSM4002708["sample_id"] <- c("GSM4002708")
NewMeta_GSM4002708 <- subset(Meta_GSM4002708, select = c("sample_id"))
GSM4002708 <- AddMetaData(GSM4002708, NewMeta_GSM4002708)
head(x = GSM4002708[[]])

GSM4002709 <- CreateSeuratObject(counts = mat_GSM4002709, project = "GSM4002709")
Meta_GSM4002709 <- GSM4002709@meta.data
Meta_GSM4002709["sample_id"] <- c("GSM4002709")
NewMeta_GSM4002709 <- subset(Meta_GSM4002709, select = c("sample_id"))
GSM4002709 <- AddMetaData(GSM4002709, NewMeta_GSM4002709)
head(x = GSM4002709[[]])

GSM4002710 <- CreateSeuratObject(counts = mat_GSM4002710, project = "GSM4002710")
Meta_GSM4002710 <- GSM4002710@meta.data
Meta_GSM4002710["sample_id"] <- c("GSM4002710")
NewMeta_GSM4002710 <- subset(Meta_GSM4002710, select = c("sample_id"))
GSM4002710 <- AddMetaData(GSM4002710, NewMeta_GSM4002710)
head(x = GSM4002710[[]])

GSM4002711 <- CreateSeuratObject(counts = mat_GSM4002711, project = "GSM4002711")
Meta_GSM4002711 <- GSM4002711@meta.data
Meta_GSM4002711["sample_id"] <- c("GSM4002711")
NewMeta_GSM4002711 <- subset(Meta_GSM4002711, select = c("sample_id"))
GSM4002711 <- AddMetaData(GSM4002711, NewMeta_GSM4002711)
head(x = GSM4002711[[]])

GSM4002712 <- CreateSeuratObject(counts = mat_GSM4002712, project = "GSM4002712")
Meta_GSM4002712 <- GSM4002712@meta.data
Meta_GSM4002712["sample_id"] <- c("GSM4002712")
NewMeta_GSM4002712 <- subset(Meta_GSM4002712, select = c("sample_id"))
GSM4002712 <- AddMetaData(GSM4002712, NewMeta_GSM4002712)
head(x = GSM4002712[[]])

GSM4002713 <- CreateSeuratObject(counts = mat_GSM4002713, project = "GSM4002713")
Meta_GSM4002713 <- GSM4002713@meta.data
Meta_GSM4002713["sample_id"] <- c("GSM4002713")
NewMeta_GSM4002713 <- subset(Meta_GSM4002713, select = c("sample_id"))
GSM4002713 <- AddMetaData(GSM4002713, NewMeta_GSM4002713)
head(x = GSM4002713[[]])

GSM4002714 <- CreateSeuratObject(counts = mat_GSM4002714, project = "GSM4002714")
Meta_GSM4002714 <- GSM4002714@meta.data
Meta_GSM4002714["sample_id"] <- c("GSM4002714")
NewMeta_GSM4002714 <- subset(Meta_GSM4002714, select = c("sample_id"))
GSM4002714 <- AddMetaData(GSM4002714, NewMeta_GSM4002714)
head(x = GSM4002714[[]])

GSM4002715 <- CreateSeuratObject(counts = mat_GSM4002715, project = "GSM4002715")
Meta_GSM4002715 <- GSM4002715@meta.data
Meta_GSM4002715["sample_id"] <- c("GSM4002715")
NewMeta_GSM4002715 <- subset(Meta_GSM4002715, select = c("sample_id"))
GSM4002715 <- AddMetaData(GSM4002715, NewMeta_GSM4002715)
head(x = GSM4002715[[]])

GSM4002716 <- CreateSeuratObject(counts = mat_GSM4002716, project = "GSM4002716")
Meta_GSM4002716 <- GSM4002716@meta.data
Meta_GSM4002716["sample_id"] <- c("GSM4002716")
NewMeta_GSM4002716 <- subset(Meta_GSM4002716, select = c("sample_id"))
GSM4002716 <- AddMetaData(GSM4002716, NewMeta_GSM4002716)
head(x = GSM4002716[[]])

GSM4002717 <- CreateSeuratObject(counts = mat_GSM4002717, project = "GSM4002717")
Meta_GSM4002717 <- GSM4002717@meta.data
Meta_GSM4002717["sample_id"] <- c("GSM4002717")
NewMeta_GSM4002717 <- subset(Meta_GSM4002717, select = c("sample_id"))
GSM4002717 <- AddMetaData(GSM4002717, NewMeta_GSM4002717)
head(x = GSM4002717[[]])

GSM4002718 <- CreateSeuratObject(counts = mat_GSM4002718, project = "GSM4002718")
Meta_GSM4002718 <- GSM4002718@meta.data
Meta_GSM4002718["sample_id"] <- c("GSM4002718")
NewMeta_GSM4002718 <- subset(Meta_GSM4002718, select = c("sample_id"))
GSM4002718 <- AddMetaData(GSM4002718, NewMeta_GSM4002718)
head(x = GSM4002718[[]])

#Merge all SeuratObjects
Wheeler_merge <- merge(GSM3717013, y = c(GSM3717014, GSM3717015, GSM3717016, GSM3717017, 
                                         GSM3717018, GSM3717019, GSM3717020, GSM3717021, 
                                         GSM3717022, GSM3717023, GSM3717024, GSM3717025, 
                                         GSM3717026, GSM3717027, GSM4002707, GSM4002708, 
                                         GSM4002709, GSM4002710, GSM4002711, GSM4002712, 
                                         GSM4002713, GSM4002714, GSM4002715, GSM4002716, 
                                         GSM4002717, GSM4002718), add.cell.ids = 
                         c("GSM3717013", "GSM3717014", "GSM3717015", "GSM3717016", "GSM3717017", 
                                         "GSM3717018", "GSM3717019", "GSM3717020", "GSM3717021", 
                                         "GSM3717022", "GSM3717023", "GSM3717024", "GSM3717025", 
                                         "GSM3717026", "GSM3717027", "GSM4002707", "GSM4002708", 
                                         "GSM4002709", "GSM4002710", "GSM4002711", "GSM4002712", 
                                         "GSM4002713", "GSM4002714", "GSM4002715", "GSM4002716", 
                                         "GSM4002717", "GSM4002718"))

#Convert merged SeuratObject to SCE
sce <- as.SingleCellExperiment(Wheeler_merge, assay = "RNA")

#Load metadata
md_dir <- file.path("file_path", "metadata_Wheeler_Ms_B6.xlsx")
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
saveRDS(sce, file.path("file_path", "Wheeler_Ms_B6_SCE.rds"))

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
sce <- readRDS(file.path("file_path", "Wheeler_Ms_B6_SCE.rds"))

#INTEGRATE
#Create SeuratObject
so <- CreateSeuratObject(
  counts = counts(sce),
  meta.data = data.frame(colData(sce)),
  project = "Wheeler_Ms_B6_10x_data")

#Split by sample
cells_by_sample <- split(colnames(sce), sce$sample_id)
so <- lapply(cells_by_sample, function(i)
  SubsetData(so, cells = i))

#Normalize, find variable genes, and scale
so <- lapply(so, NormalizeData, verbose = FALSE)
so <- lapply(so, FindVariableFeatures, nfeatures = 2e3,
             selection.method = "vst", do.plot = FALSE, verbose = FALSE)
so <- lapply(so, ScaleData, verbose = FALSE)

#Find anchors and integrate
as <- FindIntegrationAnchors(so, verbose = FALSE)
so <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = FALSE)

#Scale integrated data
DefaultAssay(so) <- "integrated"
so <- ScaleData(so, display.progress = FALSE)

#DIMENSION REDUCTION 
so <- RunPCA(so, npcs = 100, verbose = FALSE)
ElbowPlot(so, ndims = 100)
##Update number of PCs used
so <- RunTSNE(so, reduction = "pca", dims = seq_len(30),
              seed.use = 1, do.fast = TRUE, verbose = FALSE)
so <- RunUMAP(so, reduction = "pca", dims = seq_len(30),
              seed.use = 1, verbose = FALSE)

#CLUSTERING 
so <- FindNeighbors(so, reduction = "pca", dims = seq_len(30), verbose = FALSE)
for (res in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))
  so <- FindClusters(so, resolution = res, random.seed = 1, verbose = FALSE)

#DR COLORED BY SAMPLE, GROUP, AND CLUSTER ID
thm <- theme(aspect.ratio = 1, legend.position = "none")
ps <- lapply(c("sample_id", "group_id", "ident"), function(u) {
  p1 <- DimPlot(so, reduction = "tsne", group.by = u) + thm
  p2 <- DimPlot(so, reduction = "umap", group.by = u)
  lgd <- get_legend(p2)
  p2 <- p2 + thm
  list(p1, p2, lgd)
  plot_grid(p1, p2, lgd, nrow = 1,
            rel_widths = c(1, 1, 0.5))
})
plot_grid(plotlist = ps, ncol = 1)

#SAVE SeuratObject
saveRDS(so, file.path("file_path", "Wheeler_B6_so_PC30.rds"))

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
so <- readRDS(file.path("file_path", "Wheeler_B6_so_PC30.rds"))
sce <- as.SingleCellExperiment(so, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>% 
  mutate_if(is.character, as.factor) %>% 
  DataFrame(row.names = colnames(sce))

#Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)

so <- SetIdent(so, value = "integrated_snn_res.0.2")
so@meta.data$cluster_id <- Idents(so)
sce$cluster_id <- Idents(so)
(n_cells <- table(sce$cluster_id, sce$sample_id))
write.csv(table(sce$cluster_id, sce$sample_id), "file_path/Wheeler_Ms_B6_so_numbers.csv")

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
cs <- sample(colnames(so), 5e3)
.plot_dr <- function(so, dr, id)
  DimPlot(so, cells = cs, group.by = id, reduction = dr, pt.size = 0.4) +
  scale_color_manual(id, values = get(paste0(id, "_pal"))) + 
  guides(col = guide_legend(nrow = 10, 
                            override.aes = list(size = 3, alpha = 1))) +
  theme_void() + theme(aspect.ratio = 1)

ids <- c("cluster_id", "group_id", "sample_id")
for (id in ids) {
  cat("## ", id, "\n")
  p1 <- .plot_dr(so, "tsne", id)
  lgd <- get_legend(p1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- .plot_dr(so, "umap", id) + theme(legend.position = "none")
  ps <- plot_grid(plotlist = list(p1, p2), nrow = 1)
  p <- plot_grid(ps, lgd, nrow = 1, rel_widths = c(1, 0.2))
  print(p)
  cat("\n\n")
}

#QC METRICS CHECK
mito.genes <- grep(pattern = "^mt-", x = rownames(so@assays[["RNA"]]), value = TRUE)
percent.mito <- Matrix::colSums(so@assays[["RNA"]][mito.genes, ])/Matrix::colSums(so@assays[["RNA"]])
so$percent.mito <- percent.mito

VlnPlot(object = so, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)

#Generate summary statistics per sample
library(data.table)
library(psych)
feature_by_sample <- as.data.frame(so$nFeature_RNA, row.names = so$sample_id)
feature_by_sample_table <- describeBy(feature_by_sample, group = so$sample_id, mat = TRUE)
write.csv(feature_by_sample_table, "file_path/Wheeler_Ms_B6_so_feature_by_sample.csv")

count_by_sample <- as.data.frame(so$nCount_RNA, row.names = so$sample_id)
count_by_sample_table <- describeBy(count_by_sample, group = so$sample_id, mat = TRUE)
write.csv(count_by_sample_table, "file_path/Wheeler_Ms_B6_so_count_by_sample.csv")

#Generate summary statistics per cluster
feature_by_cluster <- as.data.frame(so$nFeature_RNA, row.names = so$cluster_id)
feature_by_cluster_table <- describeBy(feature_by_cluster, group = so$cluster_id, mat = TRUE)
write.csv(feature_by_cluster_table, "file_path/Wheeler_Ms_B6_so_feature_by_cluster.csv")

count_by_cluster <- as.data.frame(so$nCount_RNA, row.names = so$cluster_id)
count_by_cluster_table <- describeBy(count_by_cluster, group = so$cluster_id, mat = TRUE)
write.csv(count_by_cluster_table, "file_path/Wheeler_Ms_B6_sso_count_by_cluster.csv")

#DETERMINE WHAT CELL TYPES ARE PRESENT IN DATASET BEFORE MAKING ASTROCYTE-SPECIFIC SeuratObject
#Find all markers
DefaultAssay(so) <- "RNA"
so.markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so.markers, "file_path/Wheeler_Ms_B6_so_genes.csv")

#Visualize data (examples)
DefaultAssay(so) <- "RNA"
DimPlot(so, reduction = "tsne") + theme(aspect.ratio = 1) + scale_color_manual(values = CATALYST:::.cluster_cols)

DotPlot(so, features = c("Aldoc", "Aqp4", "Gfap", "Fgfr3", "Cldn10", "Gja1", "Aldh1l1", "Slc1a3", "Slc1a2", "Hist1h3e", "Tpx2", "Nusap1"))
DotPlot(so, features = c("Cldn5", "Nostrin", "Pdgfrb", "Anpep", "Flt1", "Acta2", "Pecam1", "Vwf"))
DotPlot(so, features = c("C1qb", "Tyrobp", "Ptprc"))
DotPlot(so, features = c("Snap25", "Stmn2", "Slc17a7", "Gad1"))
DotPlot(so, features = c("Pdgfra", "Cspg4", "Sox10"))
DotPlot(so, features = c("Plp1", "Mbp", "Mog", "Olig2"))

#---------------------------------------------------------------------------------------------------
#MAKE ASTROCYTE-SPECIFIC SeuratObject
#Assign cell type identity to clusters
so.renamed <- RenameIdents(so, `0` = "Macro_A", `1` = "Endo_A", `2` = "Mes-like_A", `3` = "Macro_B", `4` = "Astro_A", `5`= "Oligo", `6` = "Mes-like_B", `7`= "Endo_B", `8`= "Macro_C", `9` = "Neuro",`10` = "Macro_D",`11` = "Astro_B", `12` = "Peri", `13` = "DC", `14`= "OPC", `15` = "Astro_C", `16`= "Ependymal", `17`= "Macro_E")

#Create SeuratObject with only astrocyte clusters
so_astro <- subset(x = so.renamed, idents = c("Astro_A", "Astro_B", "Astro_C"), invert = FALSE)

#RECLUSTER ASTROCYTE-SPECIFIC SeuratObject

#Increase future's maximum allowed size of objects
options(future.globals.maxSize = 2048 * 1024 ^20)
memory.limit(size = 1000000)

#Load packages
library(cowplot)
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)

#Convert SeuratObject to SCE
sce <- as.SingleCellExperiment(so_astro, assay = "RNA")

#INTEGRATE
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
##Decrease k.filter to minimize number of astrocytes identified per sample
as <- FindIntegrationAnchors(so_astro, verbose = FALSE, k.filter = 60)
so_astro <- IntegrateData(anchorset = as, dims = seq_len(30), verbose = FALSE)

#Scale integrated data
DefaultAssay(so_astro) <- "integrated"
so_astro <- ScaleData(so_astro, display.progress = FALSE)

#DIMENSION REDUCTION
so_astro <- RunPCA(so_astro, npcs = 50, verbose = FALSE)
ElbowPlot(so_astro, ndims = 50)
##Update number of PCs used
so_astro <- RunTSNE(so_astro, reduction = "pca", dims = seq_len(10),
                    seed.use = 1, do.fast = TRUE, verbose = FALSE)
so_astro <- RunUMAP(so_astro, reduction = "pca", dims = seq_len(10),
                    seed.use = 1, verbose = FALSE)

#CLUSTERING
so_astro <- FindNeighbors(so_astro, reduction = "pca", dims = seq_len(10), verbose = FALSE)
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
saveRDS(so_astro, file.path("file_path", "Wheeler_B6_so_astro_PC10.rds"))

#---------------------------------------------------------------------------------------------------
#CLUSTER ANNOTATION FOR ASTROCYTE ONLY SeuratObject

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
so_astro <- readRDS(file.path("file_path", "Wheeler_B6_so_astro_PC10.rds"))
sce <- as.SingleCellExperiment(so_astro, assay = "RNA")
colData(sce) <- as.data.frame(colData(sce)) %>% 
  mutate_if(is.character, as.factor) %>% 
  DataFrame(row.names = colnames(sce))

#Define resolution
cluster_cols <- grep("res.[0-9]", colnames(colData(sce)), value = TRUE)
sapply(colData(sce)[cluster_cols], nlevels)

so_astro <- SetIdent(so_astro, value = "integrated_snn_res.0.3")
so_astro@meta.data$cluster_id <- Idents(so_astro)
sce$cluster_id <- Idents(so_astro)
(n_cells <- table(sce$cluster_id, sce$sample_id))
write.csv(table(sce$cluster_id, sce$sample_id), "file_path/Wheeler_B6_so_astro_numbers.csv")

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

#DETERMINE WHAT DEFINES EACH ASTROCYTE CLUSTER
#Find all markers
DefaultAssay(so_astro) <- "RNA"
so_astro.markers <- FindAllMarkers(so_astro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(so_astro.markers, "file_path/Wheeler_B6_so_astro_genes.csv")

#Visualize data (examples)
DefaultAssay(so_astro) <- "RNA"
DimPlot(so_astro, reduction = "tsne", pt.size = 0.001) + theme(aspect.ratio = 1) + scale_color_manual(values = CATALYST:::.cluster_cols)

VlnPlot(so_astro, features = c("Aqp4", "Gfap", "Fgfr3", "Cldn10", "Gja1", "Aldh1l1", "Slc1a3", "Slc1a2"), pt.size = 0)

FeaturePlot(so_astro, features = c("Igtp"), reduction = "tsne", pt.size = 0.001) + theme(aspect.ratio = 1)
