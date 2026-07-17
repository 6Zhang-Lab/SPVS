library(SPVS)

Sys.setenv(SPVS_PYTHON = "YourPath/python")

SpatialData <- spvs_read_spatial(
  input = "YourPath/SpatialData.h5ad",
  input_type = "h5ad",
  sample_id = "Sample01"
)

