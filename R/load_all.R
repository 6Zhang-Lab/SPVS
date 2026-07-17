# Source reviewer-facing functions in dependency order.

spvs_source_files <- c(
  "R/Utilities.R",
  "R/ModelDerivation.R",
  "R/SpatialDataInput.R",
  "R/MolecularPrograms.R",
  "R/PlaqueCoreModel.R",
  "R/PlaqueVulnerabilityInterface.R",
  "R/AnatomicalReference.R",
  "R/BoundaryValidation.R",
  "R/BoundaryEcology.R",
  "R/CellComposition.R",
  "R/Visualization.R",
  "R/SingleCellGrounding.R",
  "R/SPVSWorkflow.R"
)

invisible(lapply(spvs_source_files, source))
