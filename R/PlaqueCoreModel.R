
.score_programs <- function(sp, programs, prefix, min_genes = 2) {
  .pkg(c("Matrix", "dplyr"))

  data <- sp$data
  expr <- sp$expr
  rn <- rownames(expr)
  lib <- Matrix::colSums(expr)
  lib[!is.finite(lib) | lib <= 0] <- 1
  src <- list()

  for (nm in names(programs)) {
    genes <- .match_genes(programs[[nm]], rn)
    out_col <- paste0(prefix, "__", .sfile(nm))

    if (length(genes) >= min_genes) {
      mat <- expr[genes, , drop = FALSE]
      mat <- as.matrix(mat)
      storage.mode(mat) <- "numeric"
      lib_use <- lib[colnames(mat)]
      lib_use[!is.finite(lib_use) | lib_use <= 0] <- 1
      mat <- sweep(mat, 2, lib_use, "/", check.margin = FALSE) * 1e4
      mat <- log1p(mat)
      vals <- base::colMeans(mat)
      data[[out_col]] <- .gz(vals[match(data$spot_id, colnames(expr))], data$sample)
      source <- "marker_program"
    } else {
      data[[out_col]] <- NA_real_
      source <- "not_scored_insufficient_genes"
    }

    src[[length(src) + 1]] <- data.frame(
      program = nm,
      defined_genes = paste(programs[[nm]], collapse = ";"),
      matched_genes = paste(genes, collapse = ";"),
      n_matched = length(genes),
      output_column = out_col,
      source = source,
      stringsAsFactors = FALSE
    )
  }

  sp$data <- data
  attr(sp, paste0(prefix, "_sources")) <- dplyr::bind_rows(src)
  sp
}

.spvs_rescale01_by_sample <- function(x, sample) {
  out <- rep(NA_real_, length(x))
  for (s in unique(sample)) {
    idx <- which(sample == s)
    xx <- suppressWarnings(as.numeric(x[idx]))
    lo <- stats::quantile(xx, 0.02, na.rm = TRUE)
    hi <- stats::quantile(xx, 0.98, na.rm = TRUE)
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
      out[idx] <- 0.5
    } else {
      out[idx] <- pmin(1, pmax(0, (xx - lo) / (hi - lo)))
    }
  }
  out
}

.spvs_largest_components <- function(g, vertices, keep_n = 2, min_size = 10) {
  if (length(vertices) == 0) return(integer(0))
  sg <- igraph::induced_subgraph(g, vids = vertices)
  comp <- igraph::components(sg)$membership
  tab <- sort(table(comp), decreasing = TRUE)
  keep <- as.integer(names(tab)[seq_len(min(keep_n, length(tab)))])
  keep <- keep[tab[match(as.character(keep), names(tab))] >= min_size]
  if (length(keep) == 0) keep <- as.integer(names(tab)[1])
  vertices[comp %in% keep]
}

#' Delineate plaque-core-like latent state using the standalone SPVS model
#'
#' The standalone model uses fixed internal molecular programs and within-section
#' calibration. It does not require any user-provided reference boundary table.
#' @export
spvs_delineate_plaque_core <- function(sp,
                                       core_programs = spvs_default_core_programs(),
                                       core_quantile = 0.80,
                                       smooth_k = 6,
                                       min_core_spots = 10,
                                       keep_components = 2) {
  .pkg(c("dplyr", "FNN", "igraph"))

  # Score all default programs needed for the internal latent model.
  all_programs <- c(
    core_programs,
    list(
      smc_contractile_counter = c("ACTA2", "TAGLN", "MYH11", "CNN1", "MYL9"),
      stromal_counter = c("COL1A1", "COL3A1", "DCN", "LUM"),
      endothelial_counter = c("PECAM1", "VWF", "KDR", "CLDN5")
    )
  )
  sp <- .score_programs(sp, all_programs, "SPVS_core_program", min_genes = 2)
  data <- sp$data
  src <- attr(sp, "SPVS_core_program_sources")

  get_col <- function(nm) {
    hit <- src$output_column[src$program == nm & src$source == "marker_program"]
    if (length(hit) == 0) return(rep(0, nrow(data)))
    x <- suppressWarnings(as.numeric(data[[hit[1]]]))
    x[!is.finite(x)] <- 0
    x
  }

  # Fixed internal PB-SDM/SPVS latent model.
  # Positive terms capture lipid-necrotic, proteolytic, inflammatory and oxidative core-like vulnerability.
  # Counter terms reduce false core calls in contractile/stromal/endothelial wall regions.
  eta <- 0.34 * get_col("lipid_necrotic") +
    0.28 * get_col("proteolytic") +
    0.22 * get_col("inflammatory_core") +
    0.16 * get_col("hemorrhage_oxidative") -
    0.18 * get_col("smc_contractile_counter") -
    0.10 * get_col("stromal_counter") -
    0.06 * get_col("endothelial_counter")

  data$SPVS_core_latent_signal <- eta
  data$SPVS_core_probability_raw <- .sig(.gz(eta, data$sample))
  data$SPVS_core_probability <- .spvs_rescale01_by_sample(data$SPVS_core_probability_raw, data$sample)
  data$SPVS_PlaqueCore_raw <- FALSE
  data$SPVS_PlaqueCore <- FALSE
  data$SPVS_core_component <- NA_integer_

  for (s in unique(data$sample)) {
    idx <- which(data$sample == s)
    if (length(idx) <= smooth_k + 1) next

    coords <- as.matrix(data[idx, c("x", "y")])
    kk <- min(smooth_k, nrow(coords) - 1)
    nn <- FNN::get.knn(coords, k = kk)$nn.index

    # Local spatial smoothing stabilizes the model output.
    prob <- data$SPVS_core_probability[idx]
    smooth_prob <- vapply(seq_along(idx), function(i) {
      nbr <- unique(c(i, nn[i, ]))
      mean(prob[nbr], na.rm = TRUE)
    }, numeric(1))
    smooth_prob <- .spvs_rescale01_by_sample(smooth_prob, rep(s, length(smooth_prob)))
    data$SPVS_core_probability[idx] <- smooth_prob

    thr <- stats::quantile(smooth_prob, core_quantile, na.rm = TRUE)
    seed_thr <- stats::quantile(smooth_prob, max(0.90, core_quantile), na.rm = TRUE)
    raw <- smooth_prob >= thr
    seed <- smooth_prob >= seed_thr

    edges <- cbind(rep(seq_along(idx), each = kk), as.vector(nn))
    g <- igraph::as.undirected(igraph::graph_from_edgelist(edges, directed = FALSE), mode = "collapse")

    # Keep contiguous high-probability components connected to strongest seeds.
    cand <- which(raw)
    seed_vertices <- which(seed)
    keep_vertices <- .spvs_largest_components(g, cand, keep_n = keep_components, min_size = min_core_spots)
    if (length(seed_vertices) > 0 && length(keep_vertices) > 0) {
      # Prefer components that contain seed-like spots.
      sg <- igraph::induced_subgraph(g, vids = cand)
      comp <- igraph::components(sg)$membership
      cand_comp <- data.frame(vertex = cand, comp = comp)
      seed_comp <- unique(cand_comp$comp[cand_comp$vertex %in% seed_vertices])
      if (length(seed_comp) > 0) {
        keep_vertices <- cand_comp$vertex[cand_comp$comp %in% seed_comp]
        if (length(keep_vertices) < min_core_spots) {
          keep_vertices <- .spvs_largest_components(g, cand, keep_n = keep_components, min_size = min_core_spots)
        }
      }
    }

    data$SPVS_PlaqueCore_raw[idx[cand]] <- TRUE
    data$SPVS_PlaqueCore[idx[keep_vertices]] <- TRUE

    if (length(keep_vertices) > 0) {
      sg2 <- igraph::induced_subgraph(g, vids = keep_vertices)
      comp2 <- igraph::components(sg2)$membership
      data$SPVS_core_component[idx[keep_vertices]] <- comp2
    }
  }

  sp$data <- data
  attr(sp, "SPVS_model_mode") <- "standalone_fixed_internal_model"
  sp
}
