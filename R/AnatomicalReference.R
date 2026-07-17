# Morphology-anchored anatomical reference preparation.
#
# The reference is loaded from a spot-level file generated independently of
# SPVS prediction. A supplied binary interface column is used without
# re-estimation. If only anatomical compartments are available, a core-side
# interface shell can be reconstructed from spatial adjacency.

.reference_pick_col <- function(data, candidates, required = TRUE, label = "column") {
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

.read_spot_table <- function(input) {
  if (is.data.frame(input)) return(input)
  if (!is.character(input) || length(input) != 1L || !file.exists(input)) {
    stop("Input must be an existing spot-level file or data frame.", call. = FALSE)
  }
  ext <- tolower(tools::file_ext(input))
  switch(
    ext,
    csv = readr::read_csv(input, show_col_types = FALSE),
    tsv = readr::read_tsv(input, show_col_types = FALSE),
    txt = readr::read_tsv(input, show_col_types = FALSE),
    rds = readRDS(input),
    stop("Supported reference formats are CSV, TSV, TXT and RDS.", call. = FALSE)
  )
}

.reference_flag <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(ifelse(is.na(x), NA, x != 0))
  value <- toupper(trimws(as.character(x)))
  out <- rep(NA, length(value))
  out[value %in% c("TRUE", "T", "1", "YES", "Y", "PVI", "BND_IFTCORE")] <- TRUE
  out[value %in% c("FALSE", "F", "0", "NO", "N", "NON_PVI", "NON-PVI")] <- FALSE
  out
}

.canonical_anatomical_compartment <- function(x) {
  raw <- as.character(x)
  key <- tolower(gsub("[^a-z0-9]+", "", raw))
  out <- rep(NA_character_, length(key))

  out[grepl("advent", key)] <- "Adventitial compartment"
  out[grepl("media|^vsmc$|vsmc[1-4]|innermedia", key)] <- "Medial/VSMC compartment"
  out[grepl("intima|fibrous|atherosclerosis", key)] <-
    "Intimal-fibrous transition compartment"
  out[grepl("^plaque[0-9]*$|plaquecore|necroticcore|acellularmatrix", key)] <-
    "Plaque-core compartment"
  out[grepl("loosetissue|inaccessible|exclude", key)] <-
    "Excluded/technical"

  # Preserve already-standardized labels.
  exact <- raw %in% c(
    "Adventitial compartment",
    "Medial/VSMC compartment",
    "Intimal-fibrous transition compartment",
    "Plaque-core compartment",
    "Excluded/technical"
  )
  out[exact] <- raw[exact]
  out
}

.derive_core_side_shell <- function(data, shell_multiplier = 1.65) {
  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Package 'FNN' is required to derive the reference shell.", call. = FALSE)
  }

  out <- rep(NA, nrow(data))
  for (section in unique(data$sample)) {
    idx <- which(
      data$sample == section &
      is.finite(data$x) & is.finite(data$y) &
      !is.na(data$anatomical_compartment) &
      data$anatomical_compartment != "Excluded/technical"
    )
    if (length(idx) < 3L) next

    coordinates <- as.matrix(data[idx, c("x", "y")])
    nearest <- FNN::get.knn(coordinates, k = 1)$nn.dist[, 1]
    spacing <- stats::median(nearest[is.finite(nearest) & nearest > 0], na.rm = TRUE)
    if (!is.finite(spacing) || spacing <= 0) next

    compartment <- data$anatomical_compartment[idx]
    core_local <- which(compartment == "Plaque-core compartment")
    noncore_local <- which(compartment != "Plaque-core compartment")
    out[idx] <- FALSE
    if (!length(core_local) || !length(noncore_local)) next

    core_xy <- coordinates[core_local, , drop = FALSE]
    noncore_xy <- coordinates[noncore_local, , drop = FALSE]
    distance_to_noncore <- FNN::get.knnx(
      data = noncore_xy,
      query = core_xy,
      k = 1
    )$nn.dist[, 1]

    shell <- distance_to_noncore <= spacing * shell_multiplier
    out[idx[core_local[shell]]] <- TRUE
  }
  out
}

prepare_anatomical_reference <- function(input,
                                         reference_col = NULL,
                                         shell_multiplier = 1.65) {
  source_name <- if (is.character(input) && length(input) == 1L) {
    basename(input)
  } else {
    "in_memory_spot_table"
  }
  raw <- .read_spot_table(input)

  sample_col <- .reference_pick_col(
    raw,
    c("sample", "sample_id", "section", "section_id", "library_id", "orig.ident"),
    label = "sample column"
  )
  spot_col <- .reference_pick_col(
    raw,
    c("spot_id", "spot_barcode", "barcode", "_obs_index", "cell_id"),
    label = "spot/barcode column"
  )
  x_col <- .reference_pick_col(
    raw,
    c("x_plot", "x", "array_col", "spatial_x", "pxl_col_in_fullres"),
    label = "x-coordinate column"
  )
  y_col <- .reference_pick_col(
    raw,
    c("y_plot", "y", "array_row", "spatial_y", "pxl_row_in_fullres"),
    label = "y-coordinate column"
  )
  region_col <- .reference_pick_col(
    raw,
    c("anatomical_compartment", "region5_simple", "pathology_region5",
      "region4", "pathology_standard_class", "region_cluster_standardized",
      "region_cluster_raw", "region_cluster"),
    required = FALSE,
    label = "anatomical-compartment column"
  )

  if (is.null(reference_col)) {
    reference_col <- .reference_pick_col(
      raw,
      c("reference_boundary", "reference_boundary_iftcore", "BND_IFTCore",
        "transition_reference", "reference_transition"),
      required = FALSE,
      label = "reference-boundary column"
    )
  }

  reference <- data.frame(
    sample = .norm_sample(raw[[sample_col]]),
    spot_id = as.character(raw[[spot_col]]),
    x = suppressWarnings(as.numeric(raw[[x_col]])),
    y = suppressWarnings(as.numeric(raw[[y_col]])),
    stringsAsFactors = FALSE
  )
  reference$anatomical_compartment <- if (!is.na(region_col)) {
    .canonical_anatomical_compartment(raw[[region_col]])
  } else {
    NA_character_
  }

  key <- paste(reference$sample, reference$spot_id, sep = "|||")
  if (anyDuplicated(key)) {
    stop("Duplicated sample-plus-spot keys in anatomical reference: ",
         source_name, call. = FALSE)
  }

  if (!is.na(reference_col)) {
    reference$reference_boundary <- .reference_flag(raw[[reference_col]])
    definition <- paste0("Binary interface mask from column '", reference_col, "'")
  } else {
    if (all(is.na(reference$anatomical_compartment))) {
      stop("No reference-boundary or anatomical-compartment column was found.",
           call. = FALSE)
    }
    reference$reference_boundary <- .derive_core_side_shell(
      reference,
      shell_multiplier = shell_multiplier
    )
    definition <- paste0(
      "Core-side anatomical interface shell; spacing multiplier = ",
      shell_multiplier
    )
  }

  reference$reference_origin <- "morphology_anchored_anatomical_reference"
  reference$reference_definition <- definition
  reference$reference_source_file <- source_name
  reference
}
