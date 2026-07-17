source("R/BoundaryValidation.R")

Prediction <- readr::read_csv("YourPath/SPVSBoundary.csv")
Reference <- readr::read_csv("YourPath/ReferenceBoundary.csv")

BoundaryMetrics <- evaluate_boundary(Prediction, Reference)
print(BoundaryMetrics)
