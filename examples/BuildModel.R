source("R/ModelDerivation.R")

CohortEvidence <- readr::read_csv("YourPath/CohortEvidence.csv")

SPVSModel <- derive_spvs_model(
  evidence = CohortEvidence,
  stable_genes = c("RGS5", "PLN")
)

readr::write_csv(SPVSModel, "YourPath/SPVSModel.csv")

