###########################. MERFISH  ######################

#https://pachterlab.github.io/voyager/articles/vig6_merfish.html 

BiocManager::install(c("Voyager", "SFEData", "SingleCellExperiment", "SpatialExperiment", 
                       "scater", "ggplot2", "patchwork", "stringr", "spdep", 
                       "BiocParallel", "BiocSingular", "gstat", "BiocNeighbors", 
                       "sf", "automap"), force = TRUE)


library(Voyager)
library(SFEData)
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater) 
library(ggplot2)
library(patchwork)
library(stringr)
library(spdep)
library(BiocParallel)
library(BiocSingular)
library(gstat)
library(BiocNeighbors)
library(sf)
library(automap)
theme_set(theme_bw())

sfe <- VizgenLiverData()

#count the number of cells in bins to better visualize cell density
plotCellBin2D(sfe, bins = 300, hex = TRUE)
names(colData(sfe))


#Negative controls

is_blank <- str_detect(rownames(sfe), "^Blank-")
sfe <- addPerCellQCMetrics(sfe, subset = list(blank = is_blank))
names(colData(sfe))

#the cell polygons are too small to see anyway, so plotting cell centroid points 
print(plotSpatialFeature(sfe, "nCounts", colGeometryName = "centroids",
                         scattermore = TRUE))

#inspect cell-cell contiguity
#nCounts seems to be related to cell size
#larger cells seem to have more total counts.
bbox_use <- c(xmin = 3000, xmax = 3500, ymin = 2500, ymax = 3000)
plotSpatialFeature(sfe, "nCounts", 
                   colGeometryName = "cellSeg", 
                   bbox = bbox_use)

rownames(sfe)
n_panel <- 347

get_neg_ctrl_outliers <- function(col, sfe, nmads = 3, log = FALSE) {
  inds <- colData(sfe)$nCounts > 0 & colData(sfe)[[col]] > 0
  df <- colData(sfe)[inds,]
  outlier_inds <- isOutlier(df[[col]], type = "higher", nmads = nmads, log = log)
  outliers <- rownames(df)[outlier_inds]
  col2 <- str_remove(col, "^subsets_")
  col2 <- str_remove(col2, "_percent$")
  new_colname <- paste("is", col2, "outlier", sep = "_")
  colData(sfe)[[new_colname]] <- colnames(sfe) %in% outliers
  sfe
}

sfe <- get_neg_ctrl_outliers("subsets_blank_percent", sfe, log = TRUE)

#What proportion of all cells are outliers?
mean(sfe$is_blank_outlier)

#What’s the cutoff for outlier?
min(sfe$subsets_blank_percent[sfe$is_blank_outlier])

#Remove the outliers and empty cells:
sfe <- sfe[, !sfe$is_blank_outlier & sfe$nCounts > 0]
sfe #check the N cells after outlirers removing


#spatial autocorrelation of QC metrics
colGraph(sfe, "knn5") <- findSpatialNeighbors(sfe, method = "knearneigh", 
                                              dist_type = "idw", k = 5, 
                                              style = "W")


#With the spatial neighborhood graph, compute Moran’s I for QC metrics
sfe <- colDataMoransI(sfe, c("nCounts", 
                             "nGenes", 
                             "volume"), 
                      colGraphName = "knn5")

#nCounts and nGenes have sizable negative Moran’s I’s, which is closer to 0 for volume.
colFeatureData(sfe)[c("nCounts", "nGenes", "volume"),]



sfe <- colDataUnivariate(sfe, type = "localmoran", 
                         features = c("nCounts", "nGenes", "volume"),
                         colGraphName = "knn5")

plotLocalResult(sfe, "localmoran", c("nCounts", "nGenes", "volume"),
                colGeometryName = "centroids", scattermore = TRUE,
                ncol = 2, divergent = TRUE, diverge_center = 0)
# niches around smaller blood vessels with positive local Moran’s I


#Moran’s I
sfe <- logNormCounts(sfe)

sfe <- runMoransI(sfe, BPPARAM = MulticoreParam(2))

#top genes with positive spatial autocorrelation
top_moran <- rownames(sfe)[order(rowData(sfe)$moran_sample01, decreasing = TRUE)[1:6]]
plotSpatialFeature(sfe, top_moran, colGeometryName = "centroids", scattermore = TRUE,
                   ncol = 2)

#the genes with the highest Moran’s I highlight different histological regions. 
#Some probably for zones in the hepatic lobule, and some for blood vessels. 
#interesting to compare spatial autocorrelation of marker genes among different tissues and cell types

#Negative Moran’s I means that nearby cells tend to be more dissimilar to each other


#PCA for larger datasets
sfe <- runPCA(sfe, ncomponents = 20, subset_row = !is_blank,
              exprs_values = "logcounts",
              scale = TRUE, BSPARAM = IrlbaParam())

plotDimLoadings(sfe)


#Many of these genes seem to be related to the endothelium
spatialReducedDim(sfe, "PCA", 4, colGeometryName = "centroids", scattermore = TRUE,
                  divergent = TRUE, diverge_center = 0)

#PC1 and PC4 highlight the major blood vessels, while PC2 and PC3 have less spatial structure.
