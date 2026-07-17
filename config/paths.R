# Fixed study paths used for the manuscript analysis.
# These paths intentionally reflect the analysis server rather than a portable
# end-user installation.

H5AD_FILE <- paste0(
  "/data/public/scRNA/Aging/aging_artery/",
  "Visium_Slides_Novo_Atherosclerosis.h5ad"
)

OUT_DIR <- "/work/zzh/ZKN/AS模型/PB_SDM/results/"

TABLE_DIR <- file.path(OUT_DIR, "tables")
FIGURE_DIR <- file.path(OUT_DIR, "figures")
REVIEW_DIR <- file.path(OUT_DIR, "reviewer_code")
REVIEW_TABLE_DIR <- file.path(REVIEW_DIR, "tables")
REVIEW_FIGURE_DIR <- file.path(REVIEW_DIR, "figures")

# Manuscript-analysis outputs used for reviewer-facing boundary validation.
MODEL_PREDICTION_FILE <- file.path(
  TABLE_DIR,
  "Step10F_PB_SDM_layers_by_spot.csv"
)

ANATOMICAL_REGION_FILE <- file.path(
  TABLE_DIR,
  "Step10G0_region_cluster_standardized_by_spot.csv"
)

REFERENCE_BOUNDARY_FILE <- file.path(
  TABLE_DIR,
  "Step10G0D_same_shell_boundary_by_spot.csv"
)

PREPARED_REFERENCE_FILE <- file.path(
  REVIEW_TABLE_DIR,
  "MorphologyAnchored_reference_PVI_by_spot.csv"
)

PACKAGE_PREDICTION_FILE <- file.path(
  REVIEW_TABLE_DIR,
  "SPVS_predicted_PVI_by_spot.csv"
)

dir.create(REVIEW_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REVIEW_FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
