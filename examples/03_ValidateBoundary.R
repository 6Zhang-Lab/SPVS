source("config/paths.R")
source("R/Utilities.R")
source("R/AnatomicalReference.R")
source("R/BoundaryValidation.R")

# Manuscript prediction is preferred. The package rerun is a transparent
# fallback for reviewer inspection.
PredictionInput <- if (file.exists(MODEL_PREDICTION_FILE)) {
  MODEL_PREDICTION_FILE
} else {
  PACKAGE_PREDICTION_FILE
}

ReferenceInput <- if (file.exists(PREPARED_REFERENCE_FILE)) {
  PREPARED_REFERENCE_FILE
} else {
  REFERENCE_BOUNDARY_FILE
}

if (!file.exists(PredictionInput)) {
  stop("Model prediction file not found.", call. = FALSE)
}
if (!file.exists(ReferenceInput)) {
  stop("Independent anatomical reference file not found.", call. = FALSE)
}

Prediction <- prepare_model_prediction(PredictionInput)
Reference <- prepare_anatomical_reference(ReferenceInput)

ExactValidation <- evaluate_boundary(Prediction, Reference)
DistanceValidation <- evaluate_boundary_with_tolerance(
  Prediction,
  Reference,
  tolerance_spots = 0:3
)

readr::write_csv(
  ExactValidation$joined_spots,
  file.path(REVIEW_TABLE_DIR, "boundary_validation_exact_join_audit.csv")
)
readr::write_csv(
  ExactValidation$join_qc,
  file.path(REVIEW_TABLE_DIR, "boundary_validation_join_qc.csv")
)
readr::write_csv(
  ExactValidation$section_metrics,
  file.path(REVIEW_TABLE_DIR, "boundary_validation_by_section.csv")
)
readr::write_csv(
  ExactValidation$macro_metrics,
  file.path(REVIEW_TABLE_DIR, "boundary_validation_macro_average.csv")
)
readr::write_csv(
  DistanceValidation$section_metrics,
  file.path(REVIEW_TABLE_DIR, "boundary_distance_validation_by_section.csv")
)
readr::write_csv(
  DistanceValidation$macro_metrics,
  file.path(REVIEW_TABLE_DIR, "boundary_distance_validation_macro_average.csv")
)

print(ExactValidation$macro_metrics)
print(DistanceValidation$macro_metrics)
