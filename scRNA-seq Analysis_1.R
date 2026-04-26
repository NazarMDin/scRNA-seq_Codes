library(Seurat)

counts <- Read10X(data.dir = "C:/Users/Nazar/Downloads/Study/TBA/Assignment/scRNA/DS3")

head(counts)

seurat <- CreateSeuratObject(counts = counts, project = "DS3")

seurat
head(seurat@meta.data)

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT[-\\.]")
head(seurat@meta.data)


VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)

library(patchwork)
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")

plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

seurat
seurat <- subset(seurat, subset = nCount_RNA > 450 & nFeature_RNA > 700 & nFeature_RNA < 4500 & percent.mt < 5)

seurat
head(seurat@meta.data)


seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor =  10000)
seurat
head(seurat[["RNA"]]$data)


seurat <- FindVariableFeatures(seurat, nfeatures = 3000)
seurat
top_features <- head(VariableFeatures(seurat), 10)
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot=plot1, points = top_features, repel = TRUE)
plot1 + plot2


pdf("plot1.pdf", width = 10, height = 8)
plot1
dev.off

seurat <- ScaleData(seurat)
seurat

seurat <- RunPCA(seurat, features = VariableFeatures(seurat), npcs = 50)
seurat
ElbowPlot(seurat, ndims = 50, reduction = "pca")


seurat <- FindNeighbors(seurat, dims = 1:20)
seurat <- FindClusters(seurat, resolution = 0.5)


seurat <- RunTSNE(seurat, dims = 1:20)
seurat <- RunUMAP(seurat, dims = 1:20)

plot1 <- DimPlot(seurat, reduction = "tsne", label = TRUE)
plot2 <- DimPlot(seurat, reduction = "tsne", label = TRUE)
plot1 + plot2


seurat.markers <- FindAllMarkers(seurat, only.pos = TRUE)
write.table(seurat.markers,"C:/Users/Nazar/Downloads/Study/TBA/Assignment/scRNA/markers.txt")
head(seurat.markers)

library(dplyr)

data_markers <- seurat.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

dim(data_markers)

VlnPlot(seurat, features = c("SOX2", "TTYH1"))
VlnPlot(seurat, features = c("NKG7", "PF4"), slot = "counts")


plot2 <- FeaturePlot(seurat, c("MKI67", "NES", "DEX"), ncol=3, reduction = "umap")
plot2

