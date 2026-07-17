
spvs_cell_composition <- function(ecology_result, temperature = 1.0,
                                  rank_bins = 120,
                                  rank_equal_sampling = FALSE,
                                  rank_sample_n = NULL,
                                  rank_seed = 2026) {
  .pkg(c("dplyr", "tidyr"))

  sp <- ecology_result$sp
  src <- ecology_result$sources
  cols <- src$output_column[src$source == "marker_program"]

  if (length(cols) == 0) stop("No cell-state programs available for composition.", call. = FALSE)

  data <- sp$data
  mat <- as.matrix(data[, cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  colnames(mat) <- src$program[match(cols, src$output_column)]

  # Marker-program activities are converted to relative cell-state composition by softmax.
  # This describes relative cell-state signal composition in each Visium spot, not literal cell counts.
  mat[!is.finite(mat)] <- 0
  mat <- mat / temperature
  mat <- mat - apply(mat, 1, max, na.rm = TRUE)
  exp_mat <- exp(mat)
  denom <- rowSums(exp_mat)
  denom[!is.finite(denom) | denom <= 0] <- 1
  prop <- sweep(exp_mat, 1, denom, "/", check.margin = FALSE) * 100

  prop_df <- as.data.frame(prop, check.names = FALSE)
  prop_df$sample <- data$sample
  prop_df$spot_id <- data$spot_id
  prop_df$compartment <- as.character(data$SPVS_compartment)
  prop_df$SPVS_core_probability <- data$SPVS_core_probability
  prop_df$SPVS_shell_to_PVI <- data$SPVS_shell_to_PVI
  prop_df$x <- data$x
  prop_df$y <- data$y

  long <- prop_df %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(colnames(prop)),
      names_to = "cell_state",
      values_to = "proportion"
    )

  base_rank <- data %>%
    dplyr::select(sample, spot_id, SPVS_compartment, SPVS_core_probability, SPVS_shell_to_PVI, x, y) %>%
    dplyr::mutate(compartment = as.character(SPVS_compartment)) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(
      center_x = stats::median(x, na.rm = TRUE),
      center_y = stats::median(y, na.rm = TRUE),
      angle = atan2(y - center_y, x - center_x)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      shell_safe = ifelse(is.finite(SPVS_shell_to_PVI), SPVS_shell_to_PVI, max(SPVS_shell_to_PVI, na.rm = TRUE)),
      rank_value = dplyr::case_when(
        compartment == "Outward non-core region" ~ -shell_safe + angle / 100,
        compartment == "PVI" ~ angle,
        compartment == "Inward plaque-core region" ~ angle + SPVS_core_probability / 100,
        TRUE ~ angle
      )
    )

  count_tbl <- base_rank %>%
    dplyr::count(sample, compartment, name = "n_available") %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(
      n_min_within_sample = min(n_available),
      n_target = if (isTRUE(rank_equal_sampling)) {
        if (is.null(rank_sample_n)) n_min_within_sample else pmin(as.integer(rank_sample_n), n_min_within_sample)
      } else {
        n_available
      },
      rank_mode = ifelse(isTRUE(rank_equal_sampling), "equal_random_sampling", "all_spots_ordered_binned")
    ) %>%
    dplyr::ungroup()

  base_rank <- base_rank %>%
    dplyr::left_join(count_tbl, by = c("sample", "compartment"))

  if (isTRUE(rank_equal_sampling)) {
    set.seed(rank_seed)
    keep <- rep(FALSE, nrow(base_rank))
    group_keys <- unique(base_rank[, c("sample", "compartment")])
    for (ii in seq_len(nrow(group_keys))) {
      idx <- which(base_rank$sample == group_keys$sample[ii] & base_rank$compartment == group_keys$compartment[ii])
      target <- unique(base_rank$n_target[idx])
      target <- target[is.finite(target)][1]
      if (length(target) == 0 || is.na(target) || target <= 0) next
      chosen <- sample(idx, size = min(length(idx), target), replace = FALSE)
      keep[chosen] <- TRUE
    }
    base_rank <- base_rank[keep, , drop = FALSE]
  }

  rank_base <- base_rank %>%
    dplyr::group_by(sample, compartment) %>%
    dplyr::arrange(rank_value, .by_group = TRUE) %>%
    dplyr::mutate(
      rank_index = dplyr::row_number(),
      n_rank_sampled = dplyr::n(),
      rank_percent = ifelse(n_rank_sampled > 1, (rank_index - 1) / (n_rank_sampled - 1) * 100, 50),
      rank_bin = pmin(rank_bins, pmax(1, floor(rank_percent / 100 * rank_bins) + 1)),
      rank_mid = (rank_bin - 0.5) * 100 / rank_bins
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(sample, spot_id, compartment, rank_value, rank_index, n_rank_sampled,
                  rank_percent, rank_bin, rank_mid, n_available, n_target, rank_mode)

  rank_long <- long %>%
    dplyr::inner_join(rank_base, by = c("sample", "spot_id", "compartment"))

  rank_binned <- rank_long %>%
    dplyr::group_by(sample, compartment, rank_bin, rank_mid, cell_state) %>%
    dplyr::summarise(
      proportion = mean(proportion, na.rm = TRUE),
      n_spots_bin = dplyr::n_distinct(spot_id),
      .groups = "drop"
    )

  # Renormalize each bin to exactly 100% for visually stable stacked bars.
  rank_binned <- rank_binned %>%
    dplyr::group_by(sample, compartment, rank_bin, rank_mid) %>%
    dplyr::mutate(
      bin_total = sum(proportion, na.rm = TRUE),
      proportion = ifelse(is.finite(bin_total) & bin_total > 0, proportion / bin_total * 100, proportion)
    ) %>%
    dplyr::ungroup()

  by_sample <- long %>%
    dplyr::group_by(sample, compartment, cell_state) %>%
    dplyr::summarise(
      n_spots = dplyr::n(),
      mean_proportion = mean(proportion, na.rm = TRUE),
      median_proportion = stats::median(proportion, na.rm = TRUE),
      sd_proportion = stats::sd(proportion, na.rm = TRUE),
      .groups = "drop"
    )

  summary <- by_sample %>%
    dplyr::group_by(sample, compartment, cell_state) %>%
    dplyr::summarise(
      n_spots = sum(n_spots),
      mean_proportion = mean(mean_proportion, na.rm = TRUE),
      median_proportion = mean(median_proportion, na.rm = TRUE),
      sd_proportion = mean(sd_proportion, na.rm = TRUE),
      .groups = "drop"
    )

  list(
    sp = sp,
    wide = prop_df,
    long = long,
    rank_long = rank_long,
    rank_binned = rank_binned,
    by_sample = by_sample,
    summary = summary,
    rank_sampling = count_tbl,
    rank_mode = unique(count_tbl$rank_mode)
  )
}
