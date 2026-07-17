library(SPVS)

SPVSResult <- run_spvs_on_spatial(
  input = "YourPath/SpatialData.h5ad",
  input_type = "h5ad",
  sample_id = "Sample01",
  out_dir = "YourPath/SPVS_results"
)

