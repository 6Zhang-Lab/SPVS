source("config/paths.R")

if (!file.exists(H5AD_FILE)) {
  stop("Study H5AD file not found: ", H5AD_FILE, call. = FALSE)
}
if (!requireNamespace("SPVS", quietly = TRUE)) {
  stop(
    "Install the bundled package first: ",
    "install.packages('SPVS_FINAL.tar.gz', repos = NULL, type = 'source')",
    call. = FALSE
  )
}

SPVSResult <- SPVS::run_spvs_on_spatial(
  input = H5AD_FILE,
  input_type = "h5ad",
  out_dir = file.path(REVIEW_DIR, "SPVS_application"),
  core_quantile = 0.80,
  k = 6,
  shell_steps = 5,
  boundary_shell = 2,
  boundary_gradient_quantile = 1
)

Prediction <- SPVSResult$spatial$data
Prediction$predicted_boundary <- as.logical(Prediction$PVI)
Prediction <- Prediction[, c(
  "sample", "spot_id", "x", "y",
  "SPVS_core_probability", "predicted_boundary"
)]

readr::write_csv(Prediction, PACKAGE_PREDICTION_FILE)
message("Saved model-derived PVI prediction: ", PACKAGE_PREDICTION_FILE)
