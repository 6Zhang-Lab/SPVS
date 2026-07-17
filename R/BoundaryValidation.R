# Evaluate a predicted boundary against spot annotations or a reference boundary.
#
# Example input:
# Prediction <- readr::read_csv("YourPath/SPVSBoundary.csv")
# Reference <- readr::read_csv("YourPath/ReferenceBoundary.csv")
# Metrics <- evaluate_boundary(Prediction, Reference)

.binary_metrics <- function(predicted, reference) {
  predicted <- as.logical(predicted)
  reference <- as.logical(reference)
  tp <- sum(predicted & reference, na.rm = TRUE)
  fp <- sum(predicted & !reference, na.rm = TRUE)
  fn <- sum(!predicted & reference, na.rm = TRUE)
  precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  recall <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  f1 <- if (!is.finite(precision + recall) || (precision + recall) == 0) NA_real_ else 2 * precision * recall / (precision + recall)
  data.frame(precision = precision, recall = recall, F1 = f1)
}

#' Spot-level boundary evaluation
evaluate_boundary <- function(prediction, reference,
                              sample_col = "sample", spot_col = "spot_id",
                              prediction_col = "predicted_boundary",
                              reference_col = "reference_boundary") {
  keys <- c(sample_col, spot_col)
  if (!all(c(keys, prediction_col) %in% names(prediction)) ||
      !all(c(keys, reference_col) %in% names(reference))) {
    stop("Prediction and reference tables do not contain the requested columns.", call. = FALSE)
  }
  joined <- merge(
    prediction[, c(keys, prediction_col)],
    reference[, c(keys, reference_col)],
    by = keys,
    all = FALSE
  )
  if (!nrow(joined)) stop("No matched spatial spots were found.", call. = FALSE)
  metrics <- .binary_metrics(joined[[prediction_col]], joined[[reference_col]])
  metrics$n_matched_spots <- nrow(joined)
  metrics
}

#' Distance-tolerant boundary evaluation
#'
#' A predicted boundary spot is accepted when its nearest reference-boundary
#' spot lies within the user-defined spatial tolerance.
evaluate_boundary_with_tolerance <- function(prediction, reference, tolerance,
                                             x_col = "x", y_col = "y",
                                             prediction_col = "predicted_boundary",
                                             reference_col = "reference_boundary") {
  pred_xy <- as.matrix(prediction[as.logical(prediction[[prediction_col]]), c(x_col, y_col), drop = FALSE])
  ref_xy <- as.matrix(reference[as.logical(reference[[reference_col]]), c(x_col, y_col), drop = FALSE])
  if (!nrow(pred_xy) || !nrow(ref_xy)) stop("Both boundaries must contain spatial spots.", call. = FALSE)
  if (!requireNamespace("FNN", quietly = TRUE)) stop("Package 'FNN' is required.", call. = FALSE)
  pred_distance <- FNN::get.knnx(ref_xy, pred_xy, k = 1)$nn.dist[, 1]
  ref_distance <- FNN::get.knnx(pred_xy, ref_xy, k = 1)$nn.dist[, 1]
  pred_hit <- pred_distance <= tolerance
  ref_hit <- ref_distance <= tolerance
  precision <- mean(pred_hit)
  recall <- mean(ref_hit)
  data.frame(
    tolerance = tolerance,
    precision = precision,
    recall = recall,
    F1 = if ((precision + recall) == 0) NA_real_ else 2 * precision * recall / (precision + recall)
  )
}
