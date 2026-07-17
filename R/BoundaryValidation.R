# Boundary validation against an independently prepared anatomical reference.
# The reference is never used during SPVS prediction or model fitting.

.validation_pick_col <- function(data, candidates, required = TRUE, label = "column") {
  nms <- names(data)
  lower <- setNames(nms, tolower(nms))
  for (candidate in candidates) {
    key <- tolower(candidate)
    if (key %in% names(lower)) return(lower[[key]])
  }
  if (required) {
    stop("Cannot identify ", label, ". Candidates: ",
         paste(candidates, collapse = ", "), call. = FALSE)
  }
  NA_character_
}

.validation_flag <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(ifelse(is.na(x), NA, x != 0))
  value <- toupper(trimws(as.character(x)))
  out <- rep(NA, length(value))
  out[value %in% c("TRUE", "T", "1", "YES", "Y", "PVI", "BOUNDARY", "BND_IFTCORE")] <- TRUE
  out[value %in% c("FALSE", "F", "0", "NO", "N", "NON_PVI", "NON-PVI", "OUTSIDE", "CORE")] <- FALSE
  out
}

prepare_model_prediction <- function(input, prediction_col = NULL) {
  raw <- .read_spot_table(input)
  sample_col <- .validation_pick_col(
    raw,
    c("sample", "sample_id", "section", "section_id", "library_id", "orig.ident"),
    label = "sample column"
  )
  spot_col <- .validation_pick_col(
    raw,
    c("spot_id", "spot_barcode", "barcode", "_obs_index", "cell_id"),
    label = "spot/barcode column"
  )
  x_col <- .validation_pick_col(
    raw,
    c("x_plot", "x", "array_col", "spatial_x", "pxl_col_in_fullres"),
    label = "x-coordinate column"
  )
  y_col <- .validation_pick_col(
    raw,
    c("y_plot", "y", "array_row", "spatial_y", "pxl_row_in_fullres"),
    label = "y-coordinate column"
  )
  if (is.null(prediction_col)) {
    prediction_col <- .validation_pick_col(
      raw,
      c("predicted_boundary_for_PB_SDM", "predicted_boundary",
        "predicted_boundary_mask_initial", "PVI", "BND_IFTCore"),
      label = "model-predicted PVI column"
    )
  }

  prediction <- data.frame(
    sample = .norm_sample(raw[[sample_col]]),
    spot_id = as.character(raw[[spot_col]]),
    x = suppressWarnings(as.numeric(raw[[x_col]])),
    y = suppressWarnings(as.numeric(raw[[y_col]])),
    predicted_boundary = .validation_flag(raw[[prediction_col]]),
    prediction_source_column = prediction_col,
    stringsAsFactors = FALSE
  )

  key <- paste(prediction$sample, prediction$spot_id, sep = "|||")
  if (anyDuplicated(key)) {
    stop("Duplicated sample-plus-spot keys in model prediction.", call. = FALSE)
  }
  prediction
}

.binary_metrics <- function(predicted, reference) {
  predicted <- .validation_flag(predicted)
  reference <- .validation_flag(reference)
  keep <- !is.na(predicted) & !is.na(reference)
  predicted <- predicted[keep]
  reference <- reference[keep]

  tp <- sum(predicted & reference)
  fp <- sum(predicted & !reference)
  fn <- sum(!predicted & reference)
  tn <- sum(!predicted & !reference)
  precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  recall <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  specificity <- if ((tn + fp) == 0) NA_real_ else tn / (tn + fp)
  f1 <- if (!is.finite(precision + recall) || (precision + recall) == 0) {
    NA_real_
  } else {
    2 * precision * recall / (precision + recall)
  }
  iou <- if ((tp + fp + fn) == 0) NA_real_ else tp / (tp + fp + fn)

  data.frame(
    n = length(predicted), tp = tp, fp = fp, fn = fn, tn = tn,
    precision = precision, recall = recall, specificity = specificity,
    F1 = f1, IoU = iou
  )
}

.join_prediction_reference <- function(prediction, reference) {
  required_prediction <- c("sample", "spot_id", "x", "y", "predicted_boundary")
  required_reference <- c("sample", "spot_id", "reference_boundary")
  if (!all(required_prediction %in% names(prediction))) {
    stop("Prediction table is not in canonical format.", call. = FALSE)
  }
  if (!all(required_reference %in% names(reference))) {
    stop("Reference table is not in canonical format.", call. = FALSE)
  }

  joined <- merge(
    prediction[, required_prediction],
    reference[, required_reference],
    by = c("sample", "spot_id"),
    all = FALSE,
    sort = FALSE
  )
  if (!nrow(joined)) {
    stop("No exact sample-plus-spot matches were found.", call. = FALSE)
  }
  joined
}

evaluate_boundary <- function(prediction, reference) {
  joined <- .join_prediction_reference(prediction, reference)
  section_metrics <- do.call(rbind, lapply(split(joined, joined$sample), function(data) {
    out <- .binary_metrics(data$predicted_boundary, data$reference_boundary)
    out$sample <- data$sample[1]
    out
  }))
  rownames(section_metrics) <- NULL

  metric_cols <- c("precision", "recall", "specificity", "F1", "IoU")
  macro <- data.frame(
    n = sum(section_metrics$n),
    tp = sum(section_metrics$tp),
    fp = sum(section_metrics$fp),
    fn = sum(section_metrics$fn),
    tn = sum(section_metrics$tn),
    precision = mean(section_metrics$precision, na.rm = TRUE),
    recall = mean(section_metrics$recall, na.rm = TRUE),
    specificity = mean(section_metrics$specificity, na.rm = TRUE),
    F1 = mean(section_metrics$F1, na.rm = TRUE),
    IoU = mean(section_metrics$IoU, na.rm = TRUE),
    sample = "Macro_average"
  )
  for (column in metric_cols) {
    if (!is.finite(macro[[column]])) macro[[column]] <- NA_real_
  }

  list(
    joined_spots = joined,
    join_qc = data.frame(
      n_prediction_spots = nrow(prediction),
      n_reference_spots = nrow(reference),
      n_exactly_matched_spots = nrow(joined),
      prediction_match_fraction = nrow(joined) / nrow(prediction),
      reference_match_fraction = nrow(joined) / nrow(reference)
    ),
    section_metrics = section_metrics,
    macro_metrics = macro
  )
}

.section_spacing <- function(coordinates) {
  if (nrow(coordinates) < 2L) return(NA_real_)
  nearest <- FNN::get.knn(coordinates, k = 1)$nn.dist[, 1]
  stats::median(nearest[is.finite(nearest) & nearest > 0], na.rm = TRUE)
}

evaluate_boundary_with_tolerance <- function(prediction, reference,
                                             tolerance_spots = 0:3) {
  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Package 'FNN' is required for distance-aware validation.", call. = FALSE)
  }
  joined <- .join_prediction_reference(prediction, reference)

  per_section <- lapply(split(joined, joined$sample), function(data) {
    data <- data[
      is.finite(data$x) & is.finite(data$y),
      ,
      drop = FALSE
    ]
    if (nrow(data) < 2L) return(NULL)
    all_xy <- as.matrix(data[, c("x", "y")])
    spacing <- .section_spacing(all_xy)
    pred_xy <- as.matrix(data[.validation_flag(data$predicted_boundary) %in% TRUE, c("x", "y"), drop = FALSE])
    ref_xy <- as.matrix(data[.validation_flag(data$reference_boundary) %in% TRUE, c("x", "y"), drop = FALSE])
    if (!is.finite(spacing) || spacing <= 0 || !nrow(pred_xy) || !nrow(ref_xy)) {
      return(NULL)
    }

    pred_to_ref <- FNN::get.knnx(ref_xy, pred_xy, k = 1)$nn.dist[, 1] / spacing
    ref_to_pred <- FNN::get.knnx(pred_xy, ref_xy, k = 1)$nn.dist[, 1] / spacing
    symmetric <- c(pred_to_ref, ref_to_pred)

    do.call(rbind, lapply(tolerance_spots, function(tolerance) {
      precision <- mean(pred_to_ref <= tolerance)
      recall <- mean(ref_to_pred <= tolerance)
      data.frame(
        sample = data$sample[1],
        tolerance_spots = tolerance,
        precision = precision,
        recall = recall,
        F1 = if ((precision + recall) == 0) NA_real_ else
          2 * precision * recall / (precision + recall),
        ASSD_spots = mean(symmetric),
        HD95_spots = as.numeric(stats::quantile(symmetric, 0.95, na.rm = TRUE)),
        spot_spacing = spacing
      )
    }))
  })

  section_metrics <- do.call(rbind, per_section)
  if (is.null(section_metrics) || !nrow(section_metrics)) {
    stop("No evaluable sections contained both predicted and reference PVI.",
         call. = FALSE)
  }
  rownames(section_metrics) <- NULL

  macro_metrics <- stats::aggregate(
    section_metrics[, c("precision", "recall", "F1", "ASSD_spots", "HD95_spots")],
    by = list(tolerance_spots = section_metrics$tolerance_spots),
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  macro_metrics$sample <- "Macro_average"

  list(section_metrics = section_metrics, macro_metrics = macro_metrics)
}
