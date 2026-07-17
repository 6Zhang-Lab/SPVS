source("R/ModelDerivation.R")
source("R/SingleCellGrounding.R")

SingleCell <- readr::read_rds("YourPath/SingleCellReference.rds.gz")
SPVSModel <- readr::read_csv("YourPath/SPVSModel.csv")

SingleCell <- project_spvs_to_seurat(SingleCell, SPVSModel)

