#loading appropriate libraries for job
library(dplyr)
library(Seurat)
library(patchwork)

setwd("~/Downloads/SchwannCellSingleCellAnalysis/")

# Loading P1 data
cts <- list()
cts[[1]] <- ReadMtx(mtx ="GSE138577_RAW/P1_1/GSM4113877_10X_P1_1_matrix.mtx.gz",
                    cells ="GSE138577_RAW/P1_1/GSM4113877_10X_P1_1_barcodes.tsv.gz",
                    features = "GSE138577_RAW/P1_1/GSM4113877_10X_P1_1_genes.tsv.gz")

cts[[2]] <- ReadMtx(mtx ="GSE138577_RAW/P1_2/GSM4113878_10X_P1_2_matrix.mtx.gz",
                    cells ="GSE138577_RAW/P1_2/GSM4113878_10X_P1_2_barcodes.tsv.gz",
                    features = "GSE138577_RAW/P1_2/GSM4113878_10X_P1_2_genes.tsv.gz")

cts[[3]] <- ReadMtx(mtx ="GSE138577_RAW/P1_3/GSM4113879_10X_P1_3_matrix.mtx.gz",
                    cells ="GSE138577_RAW/P1_3/GSM4113879_10X_P1_3_barcodes.tsv.gz",
                    features = "GSE138577_RAW/P1_3/GSM4113879_10X_P1_3_genes.tsv.gz")
#writing colnames for data identification
colnames(cts[[1]]) <- paste0("v1_", colnames(cts[[1]]))
colnames(cts[[2]]) <- paste0("v2_", colnames(cts[[2]]))
colnames(cts[[3]]) <- paste0("v3_", colnames(cts[[3]]))

cts <- cbind(cts[[1]], cts[[2]], cts[[3]])

# Initialize the Seurat object with the raw (non-normalized data).
P1data <- CreateSeuratObject(counts = cts, project = "p1Schwann", min.cells = 3, min.features = 200)


##### Preprocessing Data ########
P1data[["percent.mt"]] <- PercentageFeatureSet(P1data, pattern = "^MT-")

# Visualize QC metrics as a violin plot
p <- VlnPlot(P1data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

savePlot(p,filename = "ViolinPlot_P1.pdf")


P1data <- subset(P1data,
                 subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 25)

P1data <- NormalizeData(object = P1data, normalization.method = "LogNormalize",
                        scale.factor = 10000)

P1data <- FindVariableFeatures(P1data,
                               selection.method = "vst", nfeatures = 2000)



P1data <- ScaleData(P1data, features = rownames(P1data))
P1data <- RunPCA(P1data, features = VariableFeatures(object =P1data), print = TRUE,
                 pcs.print = 1:5, genes.print = 5)

#P1data <- ProjectPCA(object = P1data, do.print = FALSE)
P1data <- FindNeighbors(P1data, reduction = "pca", dims = 1:10,
                       resolution = 0.5, verbose = 0, prune.SNN = 0/100,
                       algorithm = 3, force.recalc=TRUE)

set.seed(100)
P1data <- RunTSNE(object = P1data, dims.use = 1:10, do.fast = TRUE, perplexity = 100)
saveRDS(P1data, file="P1data.rds")
