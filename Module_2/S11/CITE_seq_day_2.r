

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

options(repr.plot.width=6, repr.plot.height=6)

#load in the counts matrices
obj.rna <- read.csv(file = './GEX_counts.csv', row.names = 1)
obj.ab <- read.csv(file = './ab_counts.csv', row.names = 1)

# creates a Seurat object based on the scRNA-seq data
obj <- CreateSeuratObject(counts = obj.rna)
obj
# The object contains an assay storing RNA measurement
Assays(obj)

# create a new assay to store ab information
ab_assay <- CreateAssayObject(counts = obj.ab)

# add this assay to the previously created Seurat object
obj[["ab"]] <- ab_assay

# Validate that the object now contains multiple assays
Assays(obj)

#look into the feature names in the AB-assay data
rownames(obj[["ab"]])

# The number of features and UMIs (nFeature_RNA and nCount_RNA) are automatically calculated for every object by Seurat.
# We calculate the percentage of mitochondrial features here and store it in object metadata as `percent.mito`.
# We use raw count data since this represents non-transformed and non-log-normalized counts
# The % of UMI mapping to MT-features is a common scRNA-seq QC metric.
mito.features <- grep(pattern = "^MT-", x = rownames(x = obj), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = obj, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = obj, slot = 'counts'))

# The [[ operator can add columns to object metadata, and is a great place to stash QC stats
obj[['percent.mito']] <- percent.mito
VlnPlot(object = obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0.01)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything 
# calculated by the object, i.e. columns in object metadata, PC scores etc.
FeatureScatter(object = obj, feature1 = "nCount_RNA", feature2 = "percent.mito")

#visualise the UMI and genes distribution
FeatureScatter(object = obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#First QC step done! 
obj

## Cluster and visualize cells using the usual scRNA-seq workflow, and examine for the potential presence of batch effects.

#with the help of this command you can switch into the AB assay
DefaultAssay(obj) <- 'ab'
#normalize and scale the AB assay
obj <- NormalizeData(object = obj, assay = "ab", normalization.method = "CLR")
obj <- ScaleData(obj, assay = "ab")

#Normalize your RNA assay
DefaultAssay(obj) <- 'RNA'
obj <- NormalizeData(object = obj, normalization.method = "LogNormalize", scale.factor = 1e4)

#FindVariable Features, remove the ribosomal genes from the variable features
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
length(VariableFeatures(object = obj))#length after removal

# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(obj), 20)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(obj)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot2

#Scale RNA assay, run PCA on the RNA assay
obj <- ScaleData(obj, features = VariableFeatures(obj))
obj <- RunPCA(obj, features = VariableFeatures(object = obj))

DimPlot(obj)

#look into the ellbow plot to determine the number of the PCs
ElbowPlot(object = obj, ndims = 50)

#look into the Heatmap to determine the number of the PCs
DimHeatmap(object = obj, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(object = obj, dims = 2, cells = 500, balanced = TRUE)
DimHeatmap(object = obj, dims = 3, cells = 500, balanced = TRUE)
DimHeatmap(object = obj, dims = 4, cells = 500, balanced = TRUE)
DimHeatmap(object = obj, dims = 5, cells = 500, balanced = TRUE)
DimHeatmap(object = obj, dims = 6, cells = 500, balanced = TRUE)
DimHeatmap(object = obj, dims = 7, cells = 500, balanced = TRUE)
DimHeatmap(object = obj, dims = 8, cells = 500, balanced = TRUE)
DimHeatmap(object = obj, dims = 9, cells = 500, balanced = TRUE)
DimHeatmap(object = obj, dims = 10, cells = 500, balanced = TRUE)
DimHeatmap(object = obj, dims = 11, cells = 500, balanced = TRUE)
DimHeatmap(object = obj, dims = 12, cells = 500, balanced = TRUE)
DimHeatmap(object = obj, dims = 13, cells = 500, balanced = TRUE)
DimHeatmap(object = obj, dims = 14, cells = 500, balanced = TRUE)
DimHeatmap(object = obj, dims = 15, cells = 500, balanced = TRUE)
DimHeatmap(object = obj, dims = 16, cells = 500, balanced = TRUE)
DimHeatmap(object = obj, dims = 17, cells = 500, balanced = TRUE)
DimHeatmap(object = obj, dims = 18, cells = 500, balanced = TRUE)
DimHeatmap(object = obj, dims = 19, cells = 500, balanced = TRUE)
DimHeatmap(object = obj, dims = 20, cells = 500, balanced = TRUE)

# Cluster the cells

#find neighbours and cluster the cells
obj <- FindNeighbors(object = obj, dims = 1:10)
obj <- FindClusters(object = obj, resolution = 0.5)

## Run Non-linear dimensional reduction (tSNE)

#run tsne
obj <- RunTSNE(object = obj, dims = 1:10)

#look how it looks
DimPlot(obj, reduction = 'tsne', label = TRUE, pt.size = 1.5)

#run umap
obj <- RunUMAP(obj, dims = 1:10)

#put the AB assay aside, let's continue with the RNA assay
#with the help of this command you can switch into the RNA assay
DefaultAssay(obj) <- 'RNA'

#what a beautiful UMAP!
DimPlot(obj, reduction = 'umap', label = TRUE, pt.size = 1.5, label.size = 12)

#look into some key markers for immune cells 
DefaultAssay(obj) <- 'RNA'
genes_rna <- c('CD3E', 'CD4', 'CD8A', 'CD14', 'CD19')
#CD3E - for T Cells
#CD4  - for CD4 T cells
#CD8A - for CD8A t cells
#CD14 - for myeloid cells
#CD19 - for B cells

#nice way to iterate in plotting
for (i in 1:length(genes_rna)){
    print(FeaturePlot(obj, reduction = 'umap', features = genes_rna[i], cols = c('grey','red'), pt.size = 1.5, label = TRUE))
}

#look for the markers in the AB assay
DefaultAssay(obj) <- 'ab'
rownames(obj)

#look for the markers in the AB assay
genes_ab <- c('CD3', 'CD4', 'CD8', 'CD19')
for (i in 1:length(genes_ab)){
    print(FeaturePlot(obj, reduction = 'umap', features = genes_ab[i], cols = c('grey','red'), pt.size = 1.5, min.cutoff = 'q25', label = TRUE))
}

# Cluster markers

#swith to the RNA assay
DefaultAssay(obj) <- 'RNA'

#look into the DE expression of each cluster
featuresobj <- rownames(obj)
markers.remove <- grep(pattern = "^TRAV|^TRBV|^TRGV|^TRDV|^RPL|^RPS", x = rownames(obj), value = TRUE) #remove TCR variable genes and ribosomal genes from the analysis
featuresobj <- featuresobj[!(featuresobj%in%markers.remove)]
obj.markers <- FindAllMarkers(object = obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, features = featuresobj)

DimPlot(obj, reduction = 'umap', label = TRUE, pt.size = 1.5, label.size = 12)

#visualise top 10 markers by avg Fold Change (log2)
top10 <- obj.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10

#look into the average heatmap
cluster.averages_obj <- AverageExpression(obj, assay = "RNA", return.seurat = TRUE) # , verbose = FALSE)

#look into the average heatmap
options(repr.plot.width=10, repr.plot.height=7)
DoHeatmap(cluster.averages_obj, features = top10$gene)

#look into the DotPlot
options(repr.plot.width=10, repr.plot.height=6)
top5 <- obj.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DotPlot(obj, features = unique(top5$gene), dot.scale = 8) + RotatedAxis()

#Lets find out a bit more about our T cells
#for Naiive T cells
options(repr.plot.width=8, repr.plot.height=6)
VlnPlot(obj, features = 'CCR7')
VlnPlot(obj, features = 'CD45RA', assay = 'ab')

#for cluster 5
VlnPlot(obj, features = 'S100A4') #memory subtype
VlnPlot(obj, features = 'DUSP2') #negative feedback loop
VlnPlot(obj, features = 'ZNF683') #differentiation factor

head(obj@meta.data)

#assign names for your clusters:
for(i in 1:nrow(obj@meta.data)){
    if(obj@meta.data$seurat_clusters[i] == 0){
        obj@meta.data$cluster_names[i] <- 'CD4 T-cells'
    }
    if(obj@meta.data$seurat_clusters[i] == 1){
        obj@meta.data$cluster_names[i] <- 'Myeloid'
    }
    if(obj@meta.data$seurat_clusters[i] == 2){
        obj@meta.data$cluster_names[i] <- 'B cells'
    }
    if(obj@meta.data$seurat_clusters[i] == 3){
        obj@meta.data$cluster_names[i] <- 'CD8 T-cells'
    }
    if(obj@meta.data$seurat_clusters[i] == 4){
        obj@meta.data$cluster_names[i] <- 'Myeloid-T doublets'
    }
    if(obj@meta.data$seurat_clusters[i] == 5){
        obj@meta.data$cluster_names[i] <- 'T memory'
    }
}

#check meta.data again
head(obj@meta.data)

#DimPlot new cluster names
DimPlot(obj, reduction = 'umap', label = TRUE, pt.size = 1.5, label.size = 10, group.by = 'cluster_names') #notice the group.by parameter

#you can assign idents to any column of meta.data you want!
Idents(obj) <- 'cluster_names'
levels(obj)
DimPlot(obj, reduction = 'umap', label = TRUE, pt.size = 1.5, label.size = 10)

# Visualize Ab signal

DefaultAssay(obj) <- 'ab'

#With feature umap
FeaturePlot(obj, reduction = 'umap', features = 'CD4', cols = c('grey','red'), pt.size = 1.5, label = TRUE, label.size = 10, min.cutoff = 'q25')
FeaturePlot(obj, reduction = 'umap', features = 'CD8', cols = c('grey','red'), pt.size = 1.5, label = TRUE, label.size = 10, min.cutoff = 'q25')

#With VlnPlot 
VlnPlot(obj, features = 'CD4')
VlnPlot(obj, features = 'CD8')
VlnPlot(obj, features = 'CD45RA')

#With RidgePlot 
RidgePlot(obj, features = 'CD4')
RidgePlot(obj, features = 'CD8')
RidgePlot(obj, features = 'CD45RA')

#With FeatureScatter
FeatureScatter(obj, feature1 = 'CD3', feature2  = 'CD4', pt.size = 1.5)
FeatureScatter(obj, feature1 = 'CD8', feature2  = 'CD4', pt.size = 1.5)

#gate selected cells!
plot <- FeatureScatter(obj, feature1 = 'CD8', feature2  = 'CD4', pt.size = 1.5)
select.cells <- CellSelector(plot = plot)

#look inside into the gating string
head(select.cells)

# Cluster with help of the RNA AND AB-Data

DefaultAssay(obj) <- 'ab'
#run pca on AB data
obj <- RunPCA(obj, features = rownames(obj), reduction.name = 'ab_pca')
obj

#Elbowplot to determine the number of the PCs
ElbowPlot(object = obj, ndims = 50, reduction = 'pca')
ElbowPlot(object = obj, ndims = 50, reduction = 'ab_pca')

#run the multimodal neighbours on BOTHs PCA reduction (rna and ab)
obj <- FindMultiModalNeighbors(
  obj, reduction.list = list("pca", "ab_pca"), 
  dims.list = list(1:10, 1:10), modality.weight.name = "RNA.weight"
)

#Find your clusters 
obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 3, resolution = 0.5, verbose = FALSE)

#run UMAP
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

#plot the results
DimPlot(obj, reduction = 'wnn.umap', label = TRUE, pt.size = 1.5, label.size = 10) #not the reduction parameter
DimPlot(obj, reduction = 'wnn.umap', label = TRUE, pt.size = 1.5, label.size = 10, group.by = 'cluster_names') #look on the previous clustering

#which cluster is clustererd with the help of RNA information the most?
VlnPlot(obj, features = "RNA.weight", group.by = 'seurat_clusters', sort = TRUE, pt.size = 1.5)

