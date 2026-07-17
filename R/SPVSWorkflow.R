
# The workflow predicts SPVS/PVI without reading the anatomical reference.
# Reference annotations enter only through examples/03_ValidateBoundary.R.
run_spvs_on_spatial <- function(input, out_dir = "SPVS_results", input_type = "auto",
                                sample_id = NULL, sample_col = NULL, spot_col = NULL,
                                x_col = NULL, y_col = NULL, assay = NULL, image = NULL,
                                core_quantile = 0.80, k = 6, shell_steps = 5,
                                use_histology = TRUE,
                                histology_alpha = 0.20,
                                rank_bins = 120,
                                rank_equal_sampling = FALSE,
                                rank_sample_n = NULL,
                                rank_seed = 2026,
                                boundary_shell = 2,
                                boundary_gradient_quantile = 1,
                                cell_programs = spvs_default_cell_programs(),
                                function_programs = spvs_default_function_programs(),
                                core_programs = spvs_default_core_programs()) {
  .pkg(c("dplyr", "readr", "ggplot2"))

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  table_dir <- file.path(out_dir, "tables")
  fig_dir <- file.path(out_dir, "figures")
  temp_dir <- file.path(out_dir, "tmp")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

  genes <- .gene_union(core_programs, cell_programs, function_programs)

  message("SPVS: reading spatial sample")
  sp <- spvs_read_spatial(
    input = input,
    input_type = input_type,
    sample_id = sample_id,
    sample_col = sample_col,
    spot_col = spot_col,
    x_col = x_col,
    y_col = y_col,
    assay = assay,
    image = image,
    genes = genes,
    temp_dir = temp_dir,
    use_histology = use_histology
  )

  message("SPVS: standalone plaque-core-like latent-state inference")
  sp <- spvs_delineate_plaque_core(
    sp,
    core_programs = core_programs,
    core_quantile = core_quantile,
    smooth_k = k
  )
  readr::write_csv(sp$data, file.path(table_dir, "SPVS_00_core_probability_by_spot.csv"))

  message("SPVS: standalone BND_IFTCore / PVI definition")
  sp <- spvs_define_pvi(
    sp,
    k = k,
    shell_steps = shell_steps,
    boundary_shell = boundary_shell,
    boundary_gradient_quantile = boundary_gradient_quantile
  )
  readr::write_csv(sp$data, file.path(table_dir, "SPVS_02_BND_IFTCore_PVI_by_spot.csv"))

  if (!is.null(attr(sp, "h5ad_gene_inventory"))) {
    readr::write_csv(attr(sp, "h5ad_gene_inventory"), file.path(table_dir, "SPVS_h5ad_marker_gene_inventory.csv"))
  }
  if (!is.null(attr(sp, "background_meta"))) {
    readr::write_csv(attr(sp, "background_meta"), file.path(table_dir, "SPVS_h5ad_background_meta.csv"))
  }

  message("SPVS: boundary cell-state ecology")
  ecology <- spvs_infer_boundary_ecology(sp, cell_programs = cell_programs)
  readr::write_csv(ecology$sources, file.path(table_dir, "SPVS_03_cell_state_sources.csv"))
  readr::write_csv(ecology$by_sample, file.path(table_dir, "SPVS_03_boundary_cell_state_ecology_by_sample.csv"))
  readr::write_csv(ecology$summary, file.path(table_dir, "SPVS_03_boundary_cell_state_ecology_summary.csv"))

  message("SPVS: boundary functional programs")
  functions <- spvs_infer_functional_programs(sp, function_programs = function_programs)
  readr::write_csv(functions$sources, file.path(table_dir, "SPVS_04_functional_program_sources.csv"))
  readr::write_csv(functions$by_sample, file.path(table_dir, "SPVS_04_boundary_functional_programs_by_sample.csv"))
  readr::write_csv(functions$summary, file.path(table_dir, "SPVS_04_boundary_functional_programs_summary.csv"))

  message("SPVS: cell-state composition statistics")
  composition <- spvs_cell_composition(
    ecology,
    rank_bins = rank_bins,
    rank_equal_sampling = rank_equal_sampling,
    rank_sample_n = rank_sample_n,
    rank_seed = rank_seed
  )
  readr::write_csv(composition$wide, file.path(table_dir, "SPVS_05_cell_state_composition_by_spot.csv"))
  readr::write_csv(composition$rank_long, file.path(table_dir, "SPVS_05_cell_state_composition_rank_long.csv"))
  readr::write_csv(composition$rank_binned, file.path(table_dir, "SPVS_05_cell_state_composition_rank_binned.csv"))
  readr::write_csv(composition$summary, file.path(table_dir, "SPVS_06_cell_state_proportion_statistics.csv"))
  readr::write_csv(composition$rank_sampling, file.path(table_dir, "SPVS_05_rankplot_equal_sampling_inventory.csv"))

  message("SPVS: writing figures")
  .savep(plot_spvs_pvi(sp, show_histology = use_histology, histology_alpha = histology_alpha), file.path(fig_dir, "SPVS_01_BND_IFTCore_PVI_spatial_map.pdf"), 8.5, 7.2)
  .savep(plot_spvs_core_probability(sp, show_histology = use_histology, histology_alpha = histology_alpha), file.path(fig_dir, "SPVS_02_plaque_core_probability_map.pdf"), 8.5, 6.4)
  .savep(plot_spvs_ecology(ecology), file.path(fig_dir, "SPVS_03_boundary_cell_state_ecology_dotplot.pdf"), 6.8, 4.8)
  .savep(plot_spvs_functions(functions), file.path(fig_dir, "SPVS_04_boundary_functional_program_dotplot.pdf"), 6.8, 4.8)
  .savep(plot_spvs_cell_composition_rankplot(composition), file.path(fig_dir, "SPVS_05_aligned_cell_state_composition_rankplot.pdf"), 10.5, 2.5)
  .savep(plot_spvs_cell_proportion_statistics(composition), file.path(fig_dir, "SPVS_06_cell_state_proportion_statistics.pdf"), 7.8, 5.2)

  inventory <- data.frame(
    item = c("input", "input_type", "out_dir", "n_spots", "n_samples", "core_quantile", "k", "shell_steps", "boundary_shell", "boundary_gradient_quantile", "SPVS_definition"),
    value = c(
      input, input_type, out_dir, nrow(sp$data), length(unique(sp$data$sample)),
      core_quantile, k, shell_steps, boundary_shell, boundary_gradient_quantile,
      "SPVS is an interpretable standalone Spatial Plaque Vulnerability State model; BND_IFTCore/PVI is predicted from fixed molecular programs and graph-calibrated plaque-core transition inference."
    ),
    stringsAsFactors = FALSE
  )
  readr::write_csv(inventory, file.path(table_dir, "SPVS_run_inventory.csv"))

  message("SPVS completed. Main figures:")
  for (f in c(
    "SPVS_01_BND_IFTCore_PVI_spatial_map.pdf",
    "SPVS_02_plaque_core_probability_map.pdf",
    "SPVS_03_boundary_cell_state_ecology_dotplot.pdf",
    "SPVS_04_boundary_functional_program_dotplot.pdf",
    "SPVS_05_aligned_cell_state_composition_rankplot.pdf",
    "SPVS_06_cell_state_proportion_statistics.pdf"
  )) {
    message(" - ", file.path(fig_dir, f))
  }

  invisible(list(spatial = sp, ecology = ecology, functions = functions, composition = composition, out_dir = out_dir))
}

spvs_run_demo <- function(out_dir = "/work/zzh/ZKN/AS模型/PB_SDM/results/reviewer_code/SPVS_demo_FW104302") {
  h5ad <- "/data/public/scRNA/Aging/aging_artery/Visium_Slides_Novo_Atherosclerosis.h5ad"
  if (!file.exists(h5ad)) stop("Demo H5AD not found: ", h5ad, call. = FALSE)
  Sys.setenv(SPVS_PYTHON = Sys.getenv("SPVS_PYTHON", unset = "/usr/bin/python3"))
  run_spvs_on_spatial(
    input = h5ad,
    input_type = "h5ad",
    sample_id = "FW104302",
    out_dir = out_dir,
    core_quantile = 0.80,
    k = 6,
    shell_steps = 5
  )
}
