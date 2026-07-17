source("config/paths.R")

PathInventory <- data.frame(
  object = c(
    "H5AD_FILE",
    "MODEL_PREDICTION_FILE",
    "ANATOMICAL_REGION_FILE",
    "REFERENCE_BOUNDARY_FILE",
    "OUT_DIR"
  ),
  path = c(
    H5AD_FILE,
    MODEL_PREDICTION_FILE,
    ANATOMICAL_REGION_FILE,
    REFERENCE_BOUNDARY_FILE,
    OUT_DIR
  ),
  exists = c(
    file.exists(H5AD_FILE),
    file.exists(MODEL_PREDICTION_FILE),
    file.exists(ANATOMICAL_REGION_FILE),
    file.exists(REFERENCE_BOUNDARY_FILE),
    dir.exists(OUT_DIR)
  )
)

print(PathInventory, row.names = FALSE)
readr::write_csv(PathInventory, file.path(REVIEW_TABLE_DIR, "study_path_inventory.csv"))
