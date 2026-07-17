# Core model-derivation functions consolidated from the discovery analyses.
#
# The reviewer-facing example is provided in examples/05_BuildModel.R.

#' Standardize expression within each independent cohort
#'
#' @param expression A numeric gene-by-sample matrix.
#' @return A gene-by-sample matrix of within-cohort z values.
standardize_within_cohort <- function(expression) {
  expression <- as.matrix(expression)
  center <- rowMeans(expression, na.rm = TRUE)
  scale <- apply(expression, 1, stats::sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale == 0] <- 1
  sweep(sweep(expression, 1, center, "-"), 1, scale, "/")
}

#' Convert a two-sided association test to a signed normal deviate
signed_evidence_z <- function(effect, p_value) {
  p_value <- pmin(pmax(as.numeric(p_value), .Machine$double.xmin), 1)
  sign(as.numeric(effect)) * stats::qnorm(1 - p_value / 2)
}

#' Integrate direction-aware evidence across independent cohorts
#'
#' @param evidence Data frame containing gene, cohort, effect, p_value and
#'   optionally effective_n.
#' @return Gene-level meta-evidence with BH-adjusted P values.
combine_cohort_evidence <- function(evidence) {
  required <- c("gene", "cohort", "effect", "p_value")
  if (!all(required %in% names(evidence))) {
    stop("evidence must contain: ", paste(required, collapse = ", "), call. = FALSE)
  }
  if (!"effective_n" %in% names(evidence)) evidence$effective_n <- 1
  evidence$signed_z <- signed_evidence_z(evidence$effect, evidence$p_value)

  split_gene <- split(evidence, evidence$gene)
  out <- lapply(split_gene, function(x) {
    keep <- is.finite(x$signed_z) & is.finite(x$effective_n) & x$effective_n > 0
    x <- x[keep, , drop = FALSE]
    if (!nrow(x)) return(NULL)
    w <- sqrt(x$effective_n)
    z_meta <- sum(w * x$signed_z) / sqrt(sum(w^2))
    priority <- if ("axis_priority" %in% names(x)) stats::median(x$axis_priority, na.rm = TRUE) else 1
    if (!is.finite(priority)) priority <- 1
    data.frame(
      gene = x$gene[1],
      meta_z = z_meta,
      meta_p = 2 * stats::pnorm(-abs(z_meta)),
      direction_consistency = max(mean(x$effect >= 0), mean(x$effect < 0)),
      cohorts = nrow(x),
      axis_priority = priority,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out$meta_fdr <- stats::p.adjust(out$meta_p, method = "BH")
  rownames(out) <- NULL
  out
}

#' Derive a frozen, direction-aware SPVS model
#'
#' @param evidence Cohort-level association statistics.
#' @param stable_genes Genes representing smooth-muscle stability.
#' @param fdr Maximum meta-analysis FDR.
#' @param min_consistency Minimum directional consistency across cohorts.
#' @return A gene, direction and normalized-weight model table.
derive_spvs_model <- function(evidence, stable_genes = c("RGS5", "PLN"),
                              fdr = 0.05, min_consistency = 0.60) {
  meta <- combine_cohort_evidence(evidence)
  meta <- meta[meta$meta_fdr <= fdr & meta$direction_consistency >= min_consistency, , drop = FALSE]
  if (!nrow(meta)) stop("No genes passed the model criteria.", call. = FALSE)
  meta$direction <- ifelse(toupper(meta$gene) %in% toupper(stable_genes), -1, 1)
  raw_weight <- sqrt(pmin(-log10(pmax(meta$meta_fdr, .Machine$double.xmin)), 30)) *
    meta$direction_consistency * meta$axis_priority
  meta$weight <- raw_weight / stats::median(raw_weight, na.rm = TRUE)
  meta$weight[!is.finite(meta$weight)] <- 1
  meta$weight <- pmax(pmin(meta$weight, 3), 0.25)
  meta[, c("gene", "direction", "weight", "meta_z", "meta_fdr", "cohorts")]
}

#' Select the smallest cross-cohort-stable model size
#'
#' @param performance One row per model size with k, mean_auc, worst_auc and
#'   direction_success_rate.
#' @return The selected performance row.
select_stable_model_size <- function(performance, mean_margin = 0.02,
                                     worst_margin = 0.05,
                                     min_direction_success = 1) {
  required <- c("k", "mean_auc", "worst_auc", "direction_success_rate")
  if (!all(required %in% names(performance))) {
    stop("performance must contain: ", paste(required, collapse = ", "), call. = FALSE)
  }
  full <- performance[which.max(performance$k), , drop = FALSE]
  pass <- performance$mean_auc >= full$mean_auc - mean_margin &
    performance$worst_auc >= full$worst_auc - worst_margin &
    performance$direction_success_rate >= min_direction_success
  candidates <- performance[pass, , drop = FALSE]
  if (!nrow(candidates)) return(performance[which.max(performance$mean_auc), , drop = FALSE])
  candidates[which.min(candidates$k), , drop = FALSE]
}

#' Calculate SPVS model activity using frozen directions and weights
calculate_spvs_activity <- function(expression, model) {
  z <- standardize_within_cohort(expression)
  idx <- match(toupper(model$gene), toupper(rownames(z)))
  keep <- !is.na(idx) & is.finite(model$weight) & is.finite(model$direction)
  if (!any(keep)) stop("No model genes were found in the expression matrix.", call. = FALSE)
  coefficient <- model$direction[keep] * model$weight[keep]
  as.numeric(crossprod(coefficient, z[idx[keep], , drop = FALSE]) / sum(abs(model$weight[keep])))
}
