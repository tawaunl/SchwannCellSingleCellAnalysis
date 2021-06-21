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

savePlot(filename = file.path(res))
