library(data.table)
library(SingleCellExperiment)
library(hdf5r)
library(dplyr)
library(stringr)
library(purrr)
library(anndata)


# Read annotation measurements
annotDF <- fread("annotations.tsv", sep="\t")
head(annotDF)

# Read detection measurements
detectDF <- fread("detections.tsv", sep="\t")
head(detectDF)


dim(detectDF)  # Display total number of cells

# Extract columns with mean protein expression values
protDF <- detectDF[, .SD, .SDcols = grep("Cell: Mean", names(detectDF))]
head(protDF)

# Clean up column names
names(protDF) <-  map_chr(names(protDF), ~ str_split(.x, ":")[[1]][1])
head(protDF)

# Remove DAPI and Autofluorescence
# protDF <- protDF |>  
#   select(-`Channel 1`)
# print(class(protDF))

# Sort each column in descending order
protDF <- protDF |> 
  map(~ sort(.x, decreasing = TRUE)) |> 
  as.data.frame()

print(protDF)
head(protDF)

# Spatial columns
spatialCols <- c("Centroid X µm", "Centroid Y µm", "Cell: Area µm^2", "Image", "Object ID")
spatialNames <- c("spatial_X", "spatial_Y", "Area", "ImageID", "ObjectID")
spatDF <- detectDF %>% select(all_of(spatialCols))
setnames(spatDF, spatialCols, spatialNames)

# Clean up ImageID column
spatDF$ImageID <- str_replace(spatDF$ImageID, "\\.ome.tif$", "")
print(spatDF)

adata <- AnnData(
  X = as.matrix(protDF),
  obs = spatDF
)

write_h5ad(adata, "adata.h5ad")

computeTop20Btm10 <- function(ad) {
  top20btm10DF <- data.frame(ImageID = character(), 
                             Protein = character(), top20btm10 = numeric(), 
                             stringsAsFactors = FALSE)
  
  # For each sample
  for (sID in unique(ad$obs$ImageID)) {
    subAD <- ad[ad$obs$ImageID == sID, ]
    
    # For each protein in the sub-AnnData
    for (x in colnames(subAD$X)) {
      aX <- as.vector(subAD$X[, x])
      
      # Compute the 20 largest values in aX
      top20 <- sort(aX, decreasing = TRUE)[1:20]
      
      # Compute the mean of the bottom 10th percentile of aX
      btm10 <- sort(aX)[1:floor(length(aX) * 0.1)]
      
      # Compute top20/btm10 ratio
      top20btm10 <- mean(top20) / mean(btm10)
      
      # Append result to the data frame
      top20btm10DF <- rbind(top20btm10DF, data.frame(ImageID = sID, Protein = x, top20btm10 = top20btm10))
    }
  }
  
  # # Sort the data frame in ascending order of top20btm10
  # top20btm10DF <- top20btm10DF[order(top20btm10DF$top20btm10), ]
  
  # Round the top20btm10 values to 2 decimal places
  top20btm10DF$top20btm10 <- round(top20btm10DF$top20btm10, 2)
  
  return(top20btm10DF)
}

# Example usage (replace 'adata' with your actual AnnData object)
result <- computeTop20Btm10(adata)

# Print the result
print(result)
