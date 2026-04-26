#Import library
library(Seurat)

#Read Count Data Sample 1
counts <- Read10X(data.dir = "C:/Users/Nazar/Downloads/Study/TBA/Assignment/scRNA/DS3")

#creating Seurat object
seurat <- CreateSeuratObject(counts, project = "DS3")
seurat
head(seurat@meta.data)

#Read Counts Data Sample 2
counts2 <- Read10X(data.dir = "C:/Users/Nazar/Downloads/Study/TBA/Assignment/scRNA/DS1")

#creating Seurat object 2
seurat2 <- CreateSeuratObject(counts2, project = "DS1")
seurat2
head(seurat2@meta.data)

#Merging Seurat And Seurat 2
seurat_data <- merge(seurat, seurat2, add.cells.ids=c("DS3","DS1"))
head(seurat_data@meta.data)
seurat_data <- JoinLayers(seurat_data)

#Calculate mitochondrial percentage
seurat_data[["percent.mt"]] <- PercentageFeatureSet(seurat_data, pattern = "^MT[-\\.]")
head(seurat_data@meta.data)

#Distribution of the metrics violin plot
VlnPlot(seurat_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(seurat_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)

#Plots
#Import library
library(patchwork)
plot1 <- FeatureScatter(seurat_data, feature1 = "nCount_RNA", feature2 = "percent.mt", split.by = "orig.ident")
plot2 <- FeatureScatter(seurat_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", split.by = "orig.ident")
plot1 + plot2

#Filtering Low Quality Cells
seurat_data <- subset(seurat_data, subset = nCount_RNA >200 & nFeature_RNA >500 & nFeature_RNA < 5000 & percent.mt < 5)
head(seurat_data)
seurat_data

#Normalization
seurat_data <- NormalizeData(seurat_data, normalization.method = "LogNormalize", scale.factor = 10000)
seurat <- NormalizeData(seurat_data)

#Feature Selection
seurat_data <- FindVariableFeatures(seurat_data, nfeatures = 3000)
top_features <- head(VariableFeatures(seurat_data), 10)

#Scaling
seurat_data <- ScaleData(seurat_data)

#Dimension Reduction
seurat_data <- RunPCA(seurat_data, features = VariableFeatures(seurat_data), npcs = 50)

#Elbow plot
ElbowPlot(seurat_data, ndims = 50, reduction = "pca")

seurat_data <- RunUMAP(seurat_data, dims = 1:20)

plot1 <- DimPlot(seurat_data, reduction = "umap", label = TRUE)
plot2 <- FeaturePlot(seurat_data, c("FOXG1", "DLX2", "LHX9"), ncol = 2, pt.size = 0.1)

plot1 + plot2 + plot_layout(widths = c(1.5, 2))

#Integration Using Harmony
#Import Library
library(harmony)

seurat_integrated <- RunHarmony(seurat_data, group.by.vars = "orig.ident")
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "harmony", dims = 1:20)

DimPlot(seurat_integrated, reduction = "umap", label = TRUE)

#Cluster Cells
seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "harmony", dims = 1:20)
seurat_integrated <- FindClusters(seurat_integrated, resolution = 0.2)
DimPlot(seurat_integrated, reduction = "umap", label = TRUE)

#Find Differential Expressed Markers
seurat.markers <- FindAllMarkers(seurat_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
write.table(seurat.markers, "C:/Users/Nazar/Downloads/Study/TBA/Assignment/scRNA/markers.txt", sep = "\t", row.names = FALSE)

saveRDS(seurat_integrated, "C:/Users/Nazar/Downloads/Study/TBA/Assignment/scRNA/integrated_obj.RDS")

head(seurat.markers)

library(dplyr)

data_markers <- seurat.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

head(seurat.markers)

dim(data_markers)

VlnPlot(seurat, features = c("HES11", "ZFP36L1", "FGFBP3", " TTYH1", "SFRP1"))

