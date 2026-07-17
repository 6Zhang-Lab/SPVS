
.summarize_by_compartment <- function(sp, programs, prefix, label_col) {
  .pkg(c("dplyr", "tidyr"))

  sp <- .score_programs(sp, programs, prefix, min_genes = 2)
  src <- attr(sp, paste0(prefix, "_sources"))
  cols <- src$output_column[src$source == "marker_program"]
  if (length(cols) == 0) stop("No programs scored for ", prefix, call. = FALSE)

  data <- sp$data
  long <- data %>%
    dplyr::select(sample, SPVS_compartment, dplyr::all_of(cols)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(cols), names_to = "output_column", values_to = "signal") %>%
    dplyr::left_join(src, by = "output_column") %>%
    dplyr::filter(is.finite(signal)) %>%
    dplyr::rename(!!label_col := program)

  by_sample <- long %>%
    dplyr::group_by(sample, SPVS_compartment, .data[[label_col]]) %>%
    dplyr::summarise(
      n_spots = dplyr::n(),
      mean_signal = mean(signal, na.rm = TRUE),
      median_signal = stats::median(signal, na.rm = TRUE),
      pct_positive = mean(signal > 0, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    dplyr::rename(compartment = SPVS_compartment)

  summary <- by_sample %>%
    dplyr::group_by(compartment, .data[[label_col]]) %>%
    dplyr::summarise(
      n_samples = dplyr::n_distinct(sample),
      mean_signal = mean(mean_signal, na.rm = TRUE),
      se_signal = stats::sd(mean_signal, na.rm = TRUE) / sqrt(dplyr::n_distinct(sample)),
      pct_positive = mean(pct_positive, na.rm = TRUE),
      .groups = "drop"
    )

  list(sp = sp, sources = src, long = long, by_sample = by_sample, summary = summary)
}

spvs_infer_boundary_ecology <- function(sp, cell_programs = spvs_default_cell_programs()) {
  .summarize_by_compartment(sp, cell_programs, "SPVS_cell", "cell_state")
}

spvs_infer_functional_programs <- function(sp, function_programs = spvs_default_function_programs()) {
  .summarize_by_compartment(sp, function_programs, "SPVS_function", "functional_program")
}
