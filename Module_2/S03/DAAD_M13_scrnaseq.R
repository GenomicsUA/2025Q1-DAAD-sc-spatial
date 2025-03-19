# Package installation
#remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)

# 
# remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)
# remotes::install_github("satijalab/azimuth", "seurat5", quiet = TRUE)
# remotes::install_github("satijalab/seurat-wrappers", "seurat5", quiet = TRUE)
# remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE)/
# load in the pbmc systematic comparative analysis dataset
# install dataset

# Alternatively: 
# download https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz 


pbmc.data <- Read10X(data.dir = "/Users/mkorshe/Documents/DAAD2025/scrna/filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

#https://datasets.cellxgene.cziscience.com/d60e4064-c784-436b-8e63-026d9242bfbc.rds
#Tabula Sapiens - Neural

data_2 <- readRDS('/Users/mkorshe/Documents/DAAD2025/scrna/data/d60e4064-c784-436b-8e63-026d9242bfbc.rds')

# Проаналізуйте метадану цього обєкту
str(data_2)
data_2$donor_assay
dim(data_2)

table(data_2$donor_id)
table(data_2$method)



# obj <- RunAzimuth(obj, reference = "pbmcref")

# currently, the object has two layers in the RNA assay: counts, and data

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# obj[["percent.mito"]] #<- PercentageFeatureSet(pbmc, pattern = "^MT-")

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA",'percent.mt'), ncol = 3)

###
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

all.genes <- rownames(pbmc)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pbmc <- SCTransform(pbmc)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
VizDimLoadings(pbmc, dims = 1:4, reduction = "pca")





# Testing integration methods
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)

DimPlot(pbmc, reduction = "umap")

#
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
#saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet",'other')
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# 
# 
# new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
#                      "NK", "DC", "Platelet")
# names(new.cluster.ids) <- levels(pbmc)
# pbmc <- RenameIdents(pbmc, new.cluster.ids)
# DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) #+ NoLegend()

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))



#### Other vis

features <- c("LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4")

RidgePlot(pbmc, features = features, ncol = 2)

DotPlot(pbmc, features = features) + RotatedAxis()

DoHeatmap(subset(pbmc, downsample = 100), features = features, size = 3)

FeaturePlot(pbmc, features = "MS4A1", min.cutoff = 1, max.cutoff = 3)

FeaturePlot(pbmc, features = c("MS4A1", "CD79A"), blend = TRUE)


# pbmc_data <- AggregateExpression(pbmc, assays = "RNA", return.seurat = F, group.by = c( "seurat_clusters"))

pbmc <- SCTransform(pbmc)

pbmc_data_SCT<- AggregateExpression(pbmc, assays = "RNA", return.seurat = F, group.by = c( "seurat_clusters"))
pbmc_data_SCT <- as.data.frame(pbmc_data_SCT)
pbmc3k.final.no.umap <- pbmc
pbmc3k.final.no.umap[["umap"]] <- NULL
DimPlot(pbmc3k.final.no.umap) + RotatedAxis()

# DoHeatmap now shows a grouping bar, splitting the heatmap into groups or clusters. This can
# be changed with the `group.by` parameter
DoHeatmap(pbmc, features = VariableFeatures(pbmc)[1:100], cells = 1:500, size = 4,
          angle = 90) + NoLegend() #+ DarkTheme()


plot <- FeaturePlot(pbmc, features = "MS4A1")

# Interactive 

HoverLocator(plot = plot, information = FetchData(pbmc, vars = c("ident", "PC_1", "nFeature_RNA")))

#

plot2 <- DimPlot(pbmc, reduction = "umap")
select.cells <- CellSelector(plot = plot2)
#head(select.cells)

plot3 <- FeatureScatter(pbmc, feature1 = "LYZ", feature2 = "CCL5")
B <- subset(x = pbmc, idents = "B")
FeatureScatter(B, feature1 = "LYZ", feature2 = "CCL5")
#levels(pbmc)

FeatureScatter(subset(x = pbmc, idents = "DC"), feature1 = "LYZ", feature2 = "CCL5")


bulk <- AggregateExpression(pbmc, group.by = "seurat_clusters", return.seurat = TRUE)
Cells(bulk)
gene_means_all_cells <- rowMeans(bulk@assays$RNA@layers$scale.data)

gene_means_all_cells <- rowMeans(B@assays$RNA@layers$scale.data)
hist(gene_means_all_cells)



#### CCC
remotes::install_github('saezlab/liana')
library(tidyverse)
library(magrittr)
library(liana)
show_resources()
show_methods()

testdata <- pbmc
liana_test <- liana_wrap(testdata)
liana_test %>% dplyr::glimpse()

liana_test <- liana_test %>%
  liana_aggregate()

dplyr::glimpse(liana_test)

liana_test %>%
  liana_dotplot(source_groups = c("B"),
                target_groups = c("NK", "CD8 T", "B"),
                ntop = 20)
###

liana_trunc <- liana_test %>%
  # only keep interactions concordant between methods
  filter(aggregate_rank <= 0.01) # note that these pvals are already corrected

heat_freq(liana_trunc)


### Obj2 ####


# InstallData("pbmcsca")
obj <- LoadData("pbmcsca")
obj <- subset(obj, nFeature_RNA > 1000) #Method

obj


# normalization  & DR

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

# Clustering

obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_clusters")


obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# visualize by batch and cell type annotation
# cell type annotations were previously added by Azimuth
DimPlot(obj, reduction = "umap.unintegrated", group.by = c("Method", "predicted.celltype.l2"))

DimPlot(obj, reduction = "umap.unintegrated")

obj_cel <- subset(obj,subset = Method == "CEL-Seq2")

###

InstallData("ifnb")
ifnb <- LoadData("ifnb")

ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
ifnb
# run standard anlaysis workflow
ifnb <- NormalizeData(ifnb)
ifnb <- FindVariableFeatures(ifnb)
ifnb <- ScaleData(ifnb)
ifnb <- RunPCA(ifnb)
ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
ifnb <- FindClusters(ifnb, resolution = 2, cluster.name = "unintegrated_clusters")
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(ifnb, reduction = "umap.unintegrated", group.by = c("stim", "seurat_clusters"))
###

ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                        verbose = FALSE) 
# check with 
#HarmonyIntegration
#JointPCAIntegration
#RPCAIntegration




# re-join layers after integration
ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])

ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 1)

ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca")
# Visualization
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))


##
# identify markers
Idents(ifnb) <- "seurat_annotations"
nk.markers <- FindConservedMarkers(ifnb, ident.1 = "NK", grouping.var = "stim", verbose = FALSE)
head(nk.markers)
# NEEDS TO BE FIXED AND SET ORDER CORRECTLY
Idents(ifnb) <- factor(Idents(ifnb), levels = c("pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono",
                                                "B Activated", "B", "CD8 T", "NK", "T activated", "CD4 Naive T", "CD4 Memory T"))

markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
                     "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
                     "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")
DotPlot(ifnb, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") +
  RotatedAxis()
###
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

aggregate_ifnb <- AggregateExpression(ifnb, group.by = c("seurat_annotations", "stim"), return.seurat = TRUE)
genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")

p1 <- CellScatter(aggregate_ifnb, "CD14 Mono_CTRL", "CD14 Mono_STIM", highlight = genes.to.label)
p2 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)

p3 <- CellScatter(aggregate_ifnb, "CD4 Naive T_CTRL", "CD4 Naive T_STIM", highlight = genes.to.label)
p4 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE)

p2 + p4

FeaturePlot(ifnb, features = c("CD3D", "GNLY", "IFI6"), split.by = "stim", max.cutoff = 3, cols = c("grey",
                                                                                                    "red"), reduction = "umap")

plots <- VlnPlot(ifnb, features = c("LYZ", "ISG15", "CXCL10"), split.by = "stim", group.by = "seurat_annotations",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

ifnb <- LoadData("ifnb")

# split datasets and process without integration
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
ifnb <- SCTransform(ifnb)
ifnb <- RunPCA(ifnb)
ifnb <- RunUMAP(ifnb, dims = 1:30)
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))

# integrate datasets
ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, normalization.method = "SCT", verbose = F)
ifnb <- FindNeighbors(ifnb, reduction = "integrated.dr", dims = 1:30)
ifnb <- FindClusters(ifnb, resolution = 0.6)
ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.dr")
DimPlot(ifnb, reduction = "umap", group.by = c("stim", "seurat_annotations"))


