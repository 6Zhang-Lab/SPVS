library(SPVS)

ExpressionTable <- readr::read_csv("YourPath/ExpressionMatrix.csv")
SpatialCoordinates <- readr::read_csv("YourPath/SpatialCoordinates.csv")

ExpressionMatrix <- as.matrix(ExpressionTable[, -1])
rownames(ExpressionMatrix) <- ExpressionTable[[1]]

SpatialData <- spvs_from_matrix(
  counts = ExpressionMatrix,
  coords = SpatialCoordinates,
  sample_id = "Sample01"
)
