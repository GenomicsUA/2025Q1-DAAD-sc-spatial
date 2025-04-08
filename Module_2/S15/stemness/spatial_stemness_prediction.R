# Title: Spatial stemness prediction script
# Author: Maycon Marção, Felipe Segato
# Contact: marcao.legatum@gmail.com



### Stemness prediction on Spatial Omics data (Visium example) ###
setwd("~/Desktop/spatial_course/stemness/")
library(Seurat)
library(DropletUtils)
library(hdf5r)
library(ggplot2)
library(glmGamPoi)

# 1. Load the spatial dataset

# Data downloaded from https://www.10xgenomics.com/datasets/human-brain-cancer-11-mm-capture-area-ffpe-2-standard
# It is a Human Brain Cancer, 11 mm Capture Area (FFPE)

untar("~/Desktop/spatial_course/stemness/CytAssist_11mm_FFPE_Human_Glioblastoma_spatial.tar.gz")

untar("~/Desktop/spatial_course/stemness/CytAssist_11mm_FFPE_Human_Glioblastoma_filtered_feature_bc_matrix.tar.gz")


# Create a .h5
feature_bc_matrix <- Read10X(data.dir = "~/Desktop/spatial_course/stemness/filtered_feature_bc_matrix")

write10xCounts("~/Desktop/spatial_course/stemness/filtered_feature_bc_matrix.h5", 
               feature_bc_matrix, 
               type = "HDF5",
               genome = "GRCh38",
               version = "3", overwrite = TRUE,
               gene.id = rownames(feature_bc_matrix),
               gene.symbol = rownames(feature_bc_matrix))



# Load the data into R
data.dir <- "/Users/elena/Desktop/spatial_course/stemness/"
spData <- Load10X_Spatial(
  data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL
)

plot2 <- SpatialFeaturePlot(spData, features = "nCount_Spatial", image.alpha = 0.2) + 
  theme(legend.position = "right")

print(plot2)



# 2. Transform the data before applying stemness model
# Transform the data

options(future.globals.maxSize = 3 * 1024^3)  # 3 GB

spData <- SCTransform(spData, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(spData, features = c("HPCA"))

# 3. Apply stemness model 
# Extract gene matrix 
gene_matrix <- as.matrix(spData@assays$SCT$data)
dim(gene_matrix)
gene_matrix[1:10, 1:10]

# Load Stemness model 
load("./model_RNA_MALTA.2018.Rda")
w = mm$w #extracting the stemness model weights (w) from the loaded object mm
w[1:5]
length(w) #12953 (number of genes on the model)
# Filter matrix expression by the genes on stemness model
matrix = gene_matrix
length(intersect(rownames(matrix), names(w))) #12259 

#filters the gene expression matrix to only include the genes that are present in the pretrained model (w).
predict.DATA = matrix[rownames(matrix) %in% names(w) ,]
length(rownames(predict.DATA)) # 12259
w = w[ rownames(predict.DATA) ]
length(intersect(names(w),rownames(predict.DATA))) # 12259 
w[1:5]
length(names(w)) # 12259
is.vector(w) #TRUE

# Score the Matrix `X` using Spearman correlation
s = apply( predict.DATA, 2, function(z) {cor( z, w, method="sp", use="complete.obs" )} )
s[1:5]

# Scale the scores to be between 0 and 1
s = s - min(s)
s = s / max(s)
s[1:5]
s = as.data.frame(s)
names(s) = "stemness"

# Check cells order before add stemness index into the data
cellbarcodes_seq_ori <- rownames(spData@meta.data)
cellbarcodes_seq_s <- rownames(s)

all(cellbarcodes_seq_ori %in% cellbarcodes_seq_s) #TRUE
identical(cellbarcodes_seq_ori, cellbarcodes_seq_s) #TRUE

# stemness scores to Seurat metadata
spData@meta.data$stemness <- s$stemness
SpatialFeaturePlot(spData, features = c("stemness"), image.alpha = 0)
