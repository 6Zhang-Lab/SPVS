
spvs_from_matrix <- function(counts, coords, sample_id = "sample1",
                             spot_col = NULL, x_col = NULL, y_col = NULL,
                             sample_col = NULL) {
  .pkg(c("Matrix", "dplyr"))
  counts <- as(counts, "dgCMatrix")
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("counts must have gene rownames and spot colnames.", call. = FALSE)
  }

  nms <- names(coords)
  if (is.null(spot_col)) spot_col <- .first(nms, c("spot_barcode", "barcode", "spot", "spot_id", "cell", "cell_id", "_obs_index"))
  if (is.null(x_col)) x_col <- .first(nms, c("x", "x_plot", "imagecol", "pxl_col_in_fullres", "spatial_x", "array_col"))
  if (is.null(y_col)) y_col <- .first(nms, c("y", "y_plot", "imagerow", "pxl_row_in_fullres", "spatial_y", "array_row"))
  if (is.na(spot_col) || is.na(x_col) || is.na(y_col)) {
    stop("coords must contain spot/x/y columns.", call. = FALSE)
  }

  meta <- coords %>%
    dplyr::mutate(
      spot_id = as.character(.data[[spot_col]]),
      x = suppressWarnings(as.numeric(.data[[x_col]])),
      y = suppressWarnings(as.numeric(.data[[y_col]])),
      sample = if (!is.null(sample_col) && sample_col %in% names(coords)) {
        .norm_sample(.data[[sample_col]])
      } else {
        sample_id
      }
    ) %>%
    dplyr::filter(spot_id %in% colnames(counts), is.finite(x), is.finite(y)) %>%
    dplyr::distinct(sample, spot_id, .keep_all = TRUE)

  counts <- counts[, meta$spot_id, drop = FALSE]
  list(data = meta, expr = counts, type = "SPVS_spatial")
}

# Spatial-data input functions.
#
# Example:
# SpatialData <- spvs_read_spatial(
#   input = "YourPath/SpatialData.h5ad",
#   input_type = "h5ad",
#   sample_id = "Sample01"
# )

.find_spvs_python <- function() {
  env <- Sys.getenv("SPVS_PYTHON", unset = "")
  candidates <- unique(c(
    env,
    Sys.which("python3"),
    Sys.which("python"),
    "/usr/bin/python3"
  ))
  candidates <- candidates[candidates != "" & file.exists(candidates)]

  for (py in candidates) {
    ok <- suppressWarnings(system2(
      py,
      c("-c", shQuote("import anndata, pandas, numpy, scipy; print('ok')")),
      stdout = TRUE,
      stderr = TRUE
    ))
    if (any(grepl("ok", ok))) return(py)
  }

  stop("No Python with anndata/pandas/numpy/scipy found. Set Sys.setenv(SPVS_PYTHON='YourPath/python').", call. = FALSE)
}

.spvs_write_h5ad_extract_script <- function(path) {
  py_lines <- c(
    "import argparse",
    "import numpy as np",
    "import pandas as pd",
    "from scipy import sparse",
    "import anndata as ad",
    "",
    "ap = argparse.ArgumentParser()",
    "ap.add_argument('--h5ad', required=True)",
    "ap.add_argument('--genes', required=True)",
    "ap.add_argument('--out', required=True)",
    "ap.add_argument('--inventory', required=True)",
    "ap.add_argument('--sample_id', default='')",
    "ap.add_argument('--bg', required=True)",
    "ap.add_argument('--bgmeta', required=True)",
    "args = ap.parse_args()",
    "",
    "genes = [g.strip().upper() for g in args.genes.split(',') if g.strip()]",
    "adata = ad.read_h5ad(args.h5ad)",
    "obs = adata.obs.copy()",
    "obs.insert(0, '_obs_index', obs.index.astype(str))",
    "",
    "if 'spatial' in adata.obsm.keys():",
    "    sp = np.asarray(adata.obsm['spatial'])",
    "    if sp.ndim == 2 and sp.shape[1] >= 2:",
    "        obs['spatial_x'] = sp[:, 0]",
    "        obs['spatial_y'] = sp[:, 1]",
    "",
    "sample_col = None",
    "for c in ['sample','sample_id','section','section_id','library_id','slide','slice','orig.ident','orig_ident']:",
    "    if c in obs.columns:",
    "        sample_col = c",
    "        break",
    "if sample_col is not None:",
    "    obs['SPVS_sample'] = obs[sample_col].astype(str)",
    "else:",
    "    ss = obs['_obs_index'].astype(str).str.extract(r'(FW[0-9]+(?:_v2)?)', expand=False)",
    "    obs['SPVS_sample'] = ss.fillna(args.sample_id if args.sample_id else 'sample1')",
    "",
    "choices = []",
    "if adata.raw is not None:",
    "    rv = list(map(str, adata.raw.var_names))",
    "    rn = [x.upper() for x in rv]",
    "    choices.append(('raw', sum(g in rn for g in genes), rv, rn, adata.raw.X))",
    "v = list(map(str, adata.var_names))",
    "vn = [x.upper() for x in v]",
    "choices.append(('X', sum(g in vn for g in genes), v, vn, adata.X))",
    "layer, nhit, var_names, var_norm, X = sorted(choices, key=lambda z: z[1], reverse=True)[0]",
    "",
    "idx = {}",
    "for i, g in enumerate(var_norm):",
    "    if g not in idx:",
    "        idx[g] = i",
    "",
    "found_idx = []",
    "found_genes = []",
    "inv = []",
    "for g in genes:",
    "    if g in idx:",
    "        found_idx.append(idx[g])",
    "        found_genes.append(g)",
    "        inv.append({'requested_gene': g, 'found': True, 'matched_var_name': var_names[idx[g]], 'layer_used': layer})",
    "    else:",
    "        inv.append({'requested_gene': g, 'found': False, 'matched_var_name': '', 'layer_used': layer})",
    "",
    "out = obs.copy()",
    "if len(found_idx) > 0:",
    "    mat = X[:, found_idx].toarray() if sparse.issparse(X) else np.asarray(X[:, found_idx])",
    "    for j, g in enumerate(found_genes):",
    "        out[g] = mat[:, j]",
    "",
    "# Optional H&E/background extraction from AnnData uns['spatial']",
    "# v0.2.9 writes a downsampled hex-color raster CSV, so R does not need the compiled png package.",
    "bg_rows = []",
    "try:",
    "    spatial = adata.uns.get('spatial', {}) if hasattr(adata, 'uns') else {}",
    "    lib_key = None",
    "    if isinstance(spatial, dict) and len(spatial) > 0:",
    "        keys = list(spatial.keys())",
    "        for k in keys:",
    "            if args.sample_id and str(args.sample_id) in str(k):",
    "                lib_key = k",
    "                break",
    "        if lib_key is None:",
    "            lib_key = keys[0]",
    "        lib = spatial.get(lib_key, {})",
    "        imgs = lib.get('images', {}) if isinstance(lib, dict) else {}",
    "        scalef = lib.get('scalefactors', {}) if isinstance(lib, dict) else {}",
    "        img_name = None",
    "        img = None",
    "        for nm in ['hires', 'lowres']:",
    "            if nm in imgs:",
    "                img_name = nm",
    "                img = np.asarray(imgs[nm])",
    "                break",
    "        if img is not None:",
    "            img_save = img.astype(float)",
    "            if img_save.max() <= 1.0:",
    "                img_save = np.clip(img_save * 255.0, 0, 255).astype(np.uint8)",
    "            else:",
    "                img_save = np.clip(img_save, 0, 255).astype(np.uint8)",
    "            if img_save.ndim == 2:",
    "                img_save = np.stack([img_save, img_save, img_save], axis=2)",
    "            if img_save.shape[2] > 3:",
    "                img_save = img_save[:, :, :3]",
    "            h, w = img_save.shape[0], img_save.shape[1]",
    "            max_w = 900",
    "            step = max(1, int(np.ceil(float(w) / max_w)))",
    "            img_small = img_save[::step, ::step, :]",
    "            hex_mat = np.empty((img_small.shape[0], img_small.shape[1]), dtype=object)",
    "            for rr in range(img_small.shape[0]):",
    "                row = img_small[rr]",
    "                hex_mat[rr, :] = ['#%02X%02X%02X' % (int(x[0]), int(x[1]), int(x[2])) for x in row]",
    "            pd.DataFrame(hex_mat).to_csv(args.bg, index=False, header=False)",
    "            scale_key = 'tissue_hires_scalef' if img_name == 'hires' else 'tissue_lowres_scalef'",
    "            scale = float(scalef.get(scale_key, 1.0)) if isinstance(scalef, dict) else 1.0",
    "            if not np.isfinite(scale) or scale <= 0:",
    "                scale = 1.0",
    "            full_w = float(w) / scale",
    "            full_h = float(h) / scale",
    "            bg_rows.append({'status':'written','image_file':args.bg,'library_id':lib_key,'image_name':img_name,'xmin':0,'xmax':full_w,'ymin':-full_h,'ymax':0,'image_width':w,'image_height':h,'raster_width':img_small.shape[1],'raster_height':img_small.shape[0],'scale':scale,'downsample_step':step})",
    "except Exception as e:",
    "    bg_rows.append({'status':'failed','image_file':'','library_id':'','image_name':'','xmin':0,'xmax':0,'ymin':0,'ymax':0,'image_width':0,'image_height':0,'raster_width':0,'raster_height':0,'scale':1,'downsample_step':1,'error':str(e)})",
    "if len(bg_rows) == 0:",
    "    bg_rows.append({'status':'not_found','image_file':'','library_id':'','image_name':'','xmin':0,'xmax':0,'ymin':0,'ymax':0,'image_width':0,'image_height':0,'raster_width':0,'raster_height':0,'scale':1,'downsample_step':1})",
    "pd.DataFrame(bg_rows).to_csv(args.bgmeta, index=False)",
    "",
    "pd.DataFrame(inv).to_csv(args.inventory, index=False)",
    "out.to_csv(args.out, index=False)",
    "print('SPVS_H5AD_LAYER', layer)",
    "print('SPVS_MARKERS_FOUND', len(found_idx), '/', len(genes))",
    "print('SPVS_N_OBS', out.shape[0])",
    "if len(bg_rows) > 0:",
    "    print('SPVS_BACKGROUND', bg_rows[0].get('status', 'unknown'))"
  )
  writeLines(py_lines, path)
}

.spvs_read_h5ad_python <- function(input, genes, sample_id = NULL,
                                   x_col = NULL, y_col = NULL,
                                   out_dir = tempdir(),
                                   use_histology = TRUE) {
  .pkg(c("readr", "dplyr", "Matrix", "stringr"))

  py <- .find_spvs_python()
  message("SPVS H5AD reader uses Python: ", py)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  py_script <- file.path(out_dir, "spvs_extract_h5ad.py")
  out_csv <- file.path(out_dir, "spvs_h5ad_marker_expression.csv")
  out_inv <- file.path(out_dir, "spvs_h5ad_marker_gene_inventory.csv")
  out_bg <- file.path(out_dir, "spvs_h5ad_background_hex.csv")
  out_bgmeta <- file.path(out_dir, "spvs_h5ad_background_meta.csv")
  .spvs_write_h5ad_extract_script(py_script)

  status <- system2(
    py,
    c(py_script,
      "--h5ad", input,
      "--genes", paste(unique(toupper(genes)), collapse = ","),
      "--out", out_csv,
      "--inventory", out_inv,
      "--sample_id", ifelse(is.null(sample_id), "", sample_id),
      "--bg", out_bg,
      "--bgmeta", out_bgmeta),
    stdout = TRUE,
    stderr = TRUE
  )
  message(paste(status, collapse = "\n"))

  if (!file.exists(out_csv) || !file.exists(out_inv)) {
    stop("Python H5AD extraction failed.", call. = FALSE)
  }

  dat <- readr::read_csv(out_csv, show_col_types = FALSE)
  inv <- readr::read_csv(out_inv, show_col_types = FALSE)

  if (!is.null(sample_id)) {
    keep <- .norm_sample(sample_id)
    dat$SPVS_sample_norm <- .norm_sample(dat$SPVS_sample)
    d2 <- dat[dat$SPVS_sample_norm == keep, , drop = FALSE]
    if (nrow(d2) == 0) {
      char_cols <- names(dat)[vapply(dat, is.character, logical(1))]
      hit <- Reduce(`|`, lapply(char_cols, function(cn) grepl(sample_id, dat[[cn]], fixed = TRUE)))
      d2 <- dat[hit, , drop = FALSE]
    }
    if (nrow(d2) > 0) dat <- d2
  }

  if (nrow(dat) == 0) stop("No spots retained from H5AD. Check sample_id.", call. = FALSE)

  xu <- if (!is.null(x_col) && x_col %in% names(dat)) {
    x_col
  } else {
    .first(names(dat), c("x", "x_plot", "spatial_x", "imagecol", "pxl_col_in_fullres", "array_col", "col"))
  }
  yu <- if (!is.null(y_col) && y_col %in% names(dat)) {
    y_col
  } else {
    .first(names(dat), c("y", "y_plot", "spatial_y", "imagerow", "pxl_row_in_fullres", "array_row", "row"))
  }
  if (is.na(xu) || is.na(yu)) stop("No spatial coordinates found.", call. = FALSE)

  found <- inv$requested_gene[inv$found]
  gene_cols <- intersect(found, names(dat))
  if (length(gene_cols) < 5) stop("Too few SPVS marker genes found.", call. = FALSE)

  mat <- as.matrix(t(dat[, gene_cols, drop = FALSE]))
  rownames(mat) <- gene_cols
  colnames(mat) <- as.character(dat$`_obs_index`)

  coords <- dat
  coords$spot_id <- as.character(dat$`_obs_index`)
  coords$sample <- .norm_sample(dat$SPVS_sample)
  coords$x <- suppressWarnings(as.numeric(dat[[xu]]))
  coords$y <- suppressWarnings(as.numeric(dat[[yu]]))

  sp <- spvs_from_matrix(
    counts = Matrix::Matrix(mat, sparse = TRUE),
    coords = coords,
    sample_id = ifelse(is.null(sample_id), "sample1", sample_id),
    spot_col = "spot_id",
    x_col = "x",
    y_col = "y",
    sample_col = "sample"
  )
  attr(sp, "h5ad_gene_inventory") <- inv

  if (isTRUE(use_histology) && file.exists(out_bg) && file.exists(out_bgmeta)) {
    bgmeta <- suppressWarnings(readr::read_csv(out_bgmeta, show_col_types = FALSE))
    if (nrow(bgmeta) > 0 && "status" %in% names(bgmeta) && bgmeta$status[1] == "written") {
      attr(sp, "background_image") <- out_bg
      attr(sp, "background_extent") <- as.list(bgmeta[1, c("xmin", "xmax", "ymin", "ymax")])
      attr(sp, "background_meta") <- bgmeta
    }
  }

  sp
}

spvs_read_spatial <- function(input, input_type = "auto", sample_id = NULL,
                              sample_col = NULL, spot_col = NULL,
                              x_col = NULL, y_col = NULL,
                              assay = NULL, image = NULL,
                              genes = NULL, temp_dir = tempdir(),
                              use_histology = TRUE) {
  .pkg(c("dplyr", "readr", "Matrix"))
  if (!file.exists(input)) stop("Input not found: ", input, call. = FALSE)

  if (input_type == "auto") {
    if (dir.exists(input)) {
      input_type <- "visium"
    } else {
      ext <- tolower(tools::file_ext(input))
      input_type <- if (ext == "rds") {
        "seurat"
      } else if (ext == "h5ad") {
        "h5ad"
      } else if (ext %in% c("csv", "tsv", "txt")) {
        "csv"
      } else {
        "csv"
      }
    }
  }

  if (is.null(genes)) {
    genes <- .gene_union(spvs_default_core_programs(), spvs_default_cell_programs(), spvs_default_function_programs())
  }

  if (input_type == "h5ad") {
    return(.spvs_read_h5ad_python(
      input = input,
      genes = genes,
      sample_id = sample_id,
      x_col = x_col,
      y_col = y_col,
      out_dir = temp_dir,
      use_histology = use_histology
    ))
  }

  if (input_type == "csv") {
    dat <- if (grepl("\\.tsv$|\\.txt$", input, ignore.case = TRUE)) {
      readr::read_tsv(input, show_col_types = FALSE)
    } else {
      readr::read_csv(input, show_col_types = FALSE)
    }
    nms <- names(dat)
    if (is.null(sample_col)) sample_col <- .first(nms, c("sample", "sample_id", "section", "library_id"))
    if (is.null(spot_col)) spot_col <- .first(nms, c("spot_barcode", "barcode", "spot", "spot_id", "cell", "cell_id"))
    if (is.null(x_col)) x_col <- .first(nms, c("x", "x_plot", "imagecol", "pxl_col_in_fullres", "spatial_x", "array_col"))
    if (is.null(y_col)) y_col <- .first(nms, c("y", "y_plot", "imagerow", "pxl_row_in_fullres", "spatial_y", "array_row"))
    if (is.na(spot_col) || is.na(x_col) || is.na(y_col)) stop("CSV must contain spot/x/y columns.", call. = FALSE)
    if (is.null(sample_id)) sample_id <- "sample1"

    meta <- dat %>%
      dplyr::mutate(
        sample = if (!is.na(sample_col)) .norm_sample(.data[[sample_col]]) else sample_id,
        spot_id = as.character(.data[[spot_col]]),
        x = suppressWarnings(as.numeric(.data[[x_col]])),
        y = suppressWarnings(as.numeric(.data[[y_col]]))
      ) %>%
      dplyr::filter(is.finite(x), is.finite(y))

    gene_cols <- names(dat)[toupper(names(dat)) %in% toupper(genes)]
    if (length(gene_cols) < 5) {
      meta_cols <- c(sample_col, spot_col, x_col, y_col, "sample", "spot_id", "x", "y")
      gene_cols <- setdiff(names(dat), unique(meta_cols[!is.na(meta_cols)]))
      gene_cols <- gene_cols[vapply(gene_cols, function(cc) is.numeric(dat[[cc]]) || is.integer(dat[[cc]]), logical(1))]
    }

    mat <- as.matrix(t(dat[, gene_cols, drop = FALSE]))
    colnames(mat) <- meta$spot_id
    rownames(mat) <- toupper(gene_cols)
    return(list(data = meta, expr = Matrix::Matrix(mat, sparse = TRUE), type = "SPVS_spatial"))
  }

  if (input_type %in% c("seurat", "visium")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat input requires the Seurat package.", call. = FALSE)
    obj <- if (input_type == "visium") {
      Seurat::Load10X_Spatial(input, assay = ifelse(is.null(assay), "Spatial", assay), slice = ifelse(is.null(sample_id), "slice1", sample_id))
    } else {
      readRDS(input)
    }
    if (is.null(assay)) assay <- Seurat::DefaultAssay(obj)
    cnt <- Seurat::GetAssayData(obj, assay = assay, slot = "counts")
    matched <- .match_genes(genes, rownames(cnt))
    if (length(matched) >= 5) cnt <- cnt[matched, , drop = FALSE]

    if (is.null(image)) image <- names(obj@images)[1]
    if (is.null(image) || !(image %in% names(obj@images))) stop("No spatial coordinates found in Seurat object.", call. = FALSE)
    coords <- Seurat::GetTissueCoordinates(obj, image = image)
    coords$spot_id <- rownames(coords)
    if (is.null(sample_id)) sample_id <- "sample1"

    return(spvs_from_matrix(
      counts = cnt,
      coords = coords,
      sample_id = sample_id,
      spot_col = "spot_id",
      x_col = .first(names(coords), c("imagecol", "x", "pxl_col_in_fullres", "col")),
      y_col = .first(names(coords), c("imagerow", "y", "pxl_row_in_fullres", "row"))
    ))
  }

  stop("Unsupported input_type: ", input_type, call. = FALSE)
}
