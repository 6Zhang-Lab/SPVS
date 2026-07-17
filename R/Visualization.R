
.spvs_background_layer <- function(sp, alpha = 0.42) {
  bg <- attr(sp, "background_image")
  ext <- attr(sp, "background_extent")
  if (is.null(bg) || is.null(ext) || !file.exists(bg)) return(NULL)

  # v0.2.9: background is a downsampled matrix of hex colors written by Python.
  # This avoids installing the compiled R package 'png' on servers without C compilers.
  img_df <- tryCatch(
    readr::read_csv(bg, col_names = FALSE, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  if (is.null(img_df) || nrow(img_df) == 0 || ncol(img_df) == 0) return(NULL)

  img <- as.matrix(img_df)
  img <- matrix(as.character(img), nrow = nrow(img), ncol = ncol(img))
  img[is.na(img) | img == ""] <- "#FFFFFF"

  img_alpha <- grDevices::adjustcolor(as.vector(img), alpha.f = alpha)
  img_alpha <- matrix(img_alpha, nrow = nrow(img), ncol = ncol(img))
  raster_obj <- grDevices::as.raster(img_alpha)

  ggplot2::annotation_raster(
    raster = raster_obj,
    xmin = as.numeric(ext$xmin),
    xmax = as.numeric(ext$xmax),
    ymin = as.numeric(ext$ymin),
    ymax = as.numeric(ext$ymax),
    interpolate = TRUE
  )
}

plot_spvs_pvi <- function(sp, point_size = 1.35, show_histology = TRUE, histology_alpha = 0.20) {
  .pkg(c("ggplot2", "dplyr"))
  data <- sp$data

  data$SPVS_plot_class <- dplyr::case_when(
    as.character(data$SPVS_compartment) == "PVI" ~ "Boundary",
    as.character(data$SPVS_compartment) == "Inward plaque-core region" ~ "Core",
    TRUE ~ "Outside"
  )
  data$SPVS_plot_class <- factor(data$SPVS_plot_class, levels = c("Outside", "Boundary", "Core"))

  count_df <- data %>%
    dplyr::count(SPVS_plot_class, name = "n") %>%
    dplyr::mutate(SPVS_plot_label = paste0(as.character(SPVS_plot_class), " (n=", n, ")"))

  data <- data %>%
    dplyr::left_join(count_df, by = "SPVS_plot_class")

  label_levels <- count_df$SPVS_plot_label[match(c("Outside", "Boundary", "Core"), as.character(count_df$SPVS_plot_class))]
  label_levels <- label_levels[!is.na(label_levels)]
  data$SPVS_plot_label <- factor(data$SPVS_plot_label, levels = label_levels)

  label_cols <- c()
  if (any(grepl("^Outside", label_levels))) label_cols[label_levels[grepl("^Outside", label_levels)]] <- "#C8D2DE"
  if (any(grepl("^Boundary", label_levels))) label_cols[label_levels[grepl("^Boundary", label_levels)]] <- "#F2B85B"
  if (any(grepl("^Core", label_levels))) label_cols[label_levels[grepl("^Core", label_levels)]] <- "#E85D9E"

  ns <- unique(as.character(data$sample))
  title_text <- if (length(ns) == 1) {
    paste0("PB-SDM 3-class: ", ns)
  } else {
    "PB-SDM 3-class spatial map"
  }

  p <- ggplot2::ggplot(data, ggplot2::aes(x = x, y = -y))

  if (isTRUE(show_histology)) {
    bg_layer <- .spvs_background_layer(sp, alpha = histology_alpha)
    if (!is.null(bg_layer)) p <- p + bg_layer
  }

  p <- p +
    ggplot2::geom_point(
      ggplot2::aes(fill = SPVS_plot_label),
      shape = 21,
      size = point_size,
      stroke = 0.25,
      color = "#34404A",
      alpha = 0.96
    ) +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_manual(values = label_cols, drop = FALSE) +
    ggplot2::labs(title = title_text, x = NULL, y = NULL, fill = NULL) +
    ggplot2::theme_void(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 25, margin = ggplot2::margin(b = 6)),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = ggplot2::element_text(size = 13),
      legend.key.width = grid::unit(1.25, "lines"),
      legend.key.height = grid::unit(1.25, "lines"),
      plot.margin = ggplot2::margin(8, 10, 8, 10)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 5, shape = 21, color = "#34404A")))

  if (length(ns) > 1) {
    p <- p + ggplot2::facet_wrap(~sample)
  }
  p
}

plot_spvs_core_probability <- function(sp, point_size = 1.20, show_histology = TRUE, histology_alpha = 0.20) {
  .pkg(c("ggplot2", "scales"))
  data <- sp$data
  ns <- unique(as.character(data$sample))
  title_text <- if (length(ns) == 1) {
    paste0("Plaque-core latent probability: ", ns)
  } else {
    "Plaque-core latent probability"
  }

  p <- ggplot2::ggplot(data, ggplot2::aes(x = x, y = -y))

  if (isTRUE(show_histology)) {
    bg_layer <- .spvs_background_layer(sp, alpha = histology_alpha)
    if (!is.null(bg_layer)) p <- p + bg_layer
  }

  p <- p +
    ggplot2::geom_point(ggplot2::aes(color = SPVS_core_probability), size = point_size, alpha = 0.94) +
    ggplot2::coord_equal() +
    ggplot2::scale_color_gradientn(
      colors = c("#4A7494", "#E8EEF3", "#F2B85B", "#E85D9E"),
      limits = c(0, 1),
      oob = scales::squish
    ) +
    ggplot2::labs(
      title = title_text,
      subtitle = "Marker-informed spatial core-state inference",
      x = NULL,
      y = NULL,
      color = "Core\nprobability"
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 18),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )

  if (length(ns) > 1) {
    p <- p + ggplot2::facet_wrap(~sample)
  }
  p
}

plot_spvs_ecology <- function(ecology_result) {
  .pkg(c("ggplot2", "scales"))
  dd <- ecology_result$summary
  dd$compartment <- factor(dd$compartment, levels = .compartment_levels)

  ggplot2::ggplot(dd, ggplot2::aes(x = compartment, y = cell_state)) +
    ggplot2::geom_point(ggplot2::aes(size = pct_positive, color = mean_signal), alpha = 0.92) +
    ggplot2::scale_color_gradient2(
      low = "#4A7494", mid = "white", high = "#C51B7D",
      midpoint = 0, limits = c(-1.8, 1.8), oob = scales::squish
    ) +
    ggplot2::scale_size_continuous(range = c(1.4, 6.6), limits = c(0, 100)) +
    ggplot2::labs(
      title = "Boundary cell-state ecology",
      subtitle = "Marker-informed cell-state signals across spatial compartments",
      x = NULL, y = NULL, color = "Mean z", size = "% z > 0"
    ) +
    ggplot2::theme_classic(base_size = 10.5) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8.5),
      axis.text.x = ggplot2::element_text(face = "bold", angle = 15, hjust = 1),
      axis.text.y = ggplot2::element_text(face = "bold"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.32)
    )
}

plot_spvs_functions <- function(function_result) {
  .pkg(c("ggplot2", "scales"))
  dd <- function_result$summary
  dd$compartment <- factor(dd$compartment, levels = .compartment_levels)

  ggplot2::ggplot(dd, ggplot2::aes(x = compartment, y = functional_program)) +
    ggplot2::geom_point(ggplot2::aes(size = pct_positive, color = mean_signal), alpha = 0.92) +
    ggplot2::scale_color_gradient2(
      low = "#4A7494", mid = "white", high = "#C51B7D",
      midpoint = 0, limits = c(-1.8, 1.8), oob = scales::squish
    ) +
    ggplot2::scale_size_continuous(range = c(1.4, 6.6), limits = c(0, 100)) +
    ggplot2::labs(
      title = "Boundary functional programs",
      subtitle = "Functional vulnerability programs across spatial plaque compartments",
      x = NULL, y = NULL, color = "Mean z", size = "% z > 0"
    ) +
    ggplot2::theme_classic(base_size = 10.5) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8.5),
      axis.text.x = ggplot2::element_text(face = "bold", angle = 15, hjust = 1),
      axis.text.y = ggplot2::element_text(face = "bold"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.32)
    )
}

.spvs_cell_colors <- function() {
  c(
    "Other immune" = "#D95F02",
    "B / plasma" = "#B884A9",
    "T / NK" = "#56B4E9",
    "Endothelial" = "#009E73",
    "SMC / VSMC" = "#4A7494",
    "Fibroblast / stromal" = "#E69F00",
    "Myeloid / macrophage" = "#7F3C8D"
  )
}


plot_spvs_cell_composition_rankplot <- function(composition_result) {
  .pkg(c("ggplot2", "dplyr", "scales"))
  dd <- composition_result$rank_binned
  dd$compartment <- factor(dd$compartment, levels = .compartment_levels)

  cell_levels <- names(.spvs_cell_colors())
  dd$cell_state <- factor(dd$cell_state, levels = cell_levels)

  subtitle_text <- if (any(composition_result$rank_mode == "equal_random_sampling")) {
    "Equal numbers of spots are randomly sampled from each compartment, ordered and smoothed into 0–100% rank bins"
  } else {
    "All spots are spatially ordered within each compartment and smoothed into 0–100% rank bins"
  }

  bin_width <- 100 / max(dd$rank_bin, na.rm = TRUE)

  ggplot2::ggplot(dd, ggplot2::aes(x = rank_mid, y = proportion, fill = cell_state)) +
    ggplot2::geom_col(width = bin_width * 1.02, color = NA) +
    ggplot2::facet_grid(sample ~ compartment) +
    ggplot2::scale_fill_manual(values = .spvs_cell_colors(), drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    ggplot2::scale_x_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      labels = function(x) paste0(x, "%"),
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      title = "Aligned spot-level major cell-class composition",
      subtitle = subtitle_text,
      x = "Relative spot rank within compartment",
      y = "Relative signal (%)",
      fill = NULL
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 17),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      strip.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.45),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.45),
      panel.spacing = grid::unit(0.12, "lines"),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 10),
      axis.text = ggplot2::element_text(color = "grey20"),
      axis.title = ggplot2::element_text(face = "bold")
    )
}

plot_spvs_cell_proportion_statistics <- function(composition_result) {
  .pkg(c("ggplot2", "dplyr", "scales"))
  dd <- composition_result$summary
  dd$compartment <- factor(dd$compartment, levels = .compartment_levels)

  ggplot2::ggplot(dd, ggplot2::aes(x = compartment, y = mean_proportion, fill = cell_state)) +
    ggplot2::geom_col(width = 0.72, color = "white", linewidth = 0.20) +
    ggplot2::facet_wrap(~sample, nrow = 1) +
    ggplot2::scale_fill_manual(values = .spvs_cell_colors(), drop = FALSE) +
    ggplot2::scale_y_continuous(limits = c(0, 100), expand = c(0, 0), labels = function(x) paste0(x, "%")) +
    ggplot2::labs(
      title = "Cell-state proportion statistics across SPVS compartments",
      subtitle = "Stacked bars summarize marker-derived relative cell-state composition per compartment",
      x = NULL,
      y = "Mean relative signal (%)",
      fill = NULL
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 15),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9),
      axis.text.x = ggplot2::element_text(face = "bold", angle = 15, hjust = 1),
      axis.title.y = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.35),
      legend.position = "right"
    )
}
