library(SPVS)

SpatialData <- spvs_read_spatial(
  input = "YourPath/SpatialObject.rds",
  input_type = "seurat",
  sample_id = "Sample01"
)
