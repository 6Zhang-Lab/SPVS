
.pkg <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss) > 0) stop("Missing packages: ", paste(miss, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

.first <- function(x, candidates) {
  hit <- candidates[candidates %in% x]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

.norm_sample <- function(x) {
  x <- as.character(x)
  fw <- stringr::str_extract(x, "FW[0-9]+(?:_v2)?")
  ifelse(!is.na(fw) & fw != "", fw, gsub("\\s+", "", x))
}

.z <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 1e-8) return(ifelse(is.na(x), NA_real_, 0))
  (x - m) / s
}

.gz <- function(x, group) {
  x <- suppressWarnings(as.numeric(unlist(x, use.names = FALSE)))
  group <- as.character(unlist(group, use.names = FALSE))
  out <- rep(NA_real_, length(x))
  for (g in unique(group)) {
    idx <- which(group == g)
    out[idx] <- .z(x[idx])
  }
  out
}

.sig <- function(x) 1 / (1 + exp(-x))

.bool <- function(x) {
  if (is.logical(x)) return(ifelse(is.na(x), FALSE, x))
  if (is.numeric(x)) return(ifelse(is.na(x), FALSE, x > 0))
  toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1", "YES", "Y", "PVI", "BND_IFTCORE")
}

.sfile <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

.savep <- function(p, file, w = 7, h = 5) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  tryCatch(
    ggplot2::ggsave(file, p, width = w, height = h, device = grDevices::cairo_pdf),
    error = function(e) ggplot2::ggsave(file, p, width = w, height = h, device = "pdf")
  )
  invisible(file)
}

.match_genes <- function(genes, rn) {
  genes <- unique(toupper(genes))
  rn_u <- toupper(rn)
  out <- rn[match(genes, rn_u)]
  out[!is.na(out)]
}

.compartment_levels <- c("Outward non-core region", "PVI", "Inward plaque-core region")

.gene_union <- function(core_programs, cell_programs, function_programs) {
  unique(toupper(unlist(c(core_programs, cell_programs, function_programs), use.names = FALSE)))
}
