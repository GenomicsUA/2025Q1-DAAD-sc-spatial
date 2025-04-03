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

plotGeometry(sfe, colGeometryName = "cellSeg")

plotCellBin2D(sfe, bins = 300, hex = TRUE)

names(colData(sfe))

system.time(
  print(plotSpatialFeature(sfe, "nCounts", colGeometryName = "cellSeg"))
)

system.time({
  print(plotSpatialFeature(sfe, "nCounts", colGeometryName = "centroids",
                           scattermore = TRUE))
})

bbox_use <- c(xmin = 3000, xmax = 3500, ymin = 2500, ymax = 3000)


plotSpatialFeature(sfe, "nCounts", 
                   colGeometryName = "cellSeg", 
                   bbox = bbox_use)

rownames(sfe)


n_panel <- 347


colData(sfe)$nCounts_normed <- sfe$nCounts / n_panel
colData(sfe)$nGenes_normed <- sfe$nGenes / n_panel

plotColDataHistogram(sfe, c("nCounts_normed", "nGenes_normed"))

plotSpatialFeature(sfe, "nGenes", colGeometryName = "centroids", scattermore = TRUE)

plotSpatialFeature(sfe, "volume", colGeometryName = "centroids", scattermore = TRUE)

plotColData(sfe, x="nCounts", y="nGenes", bins = 100)

plotColData(sfe, x="volume", y="nCounts", bins = 100)

plotColData(sfe, x="volume", y="nGenes", bins = 100)

