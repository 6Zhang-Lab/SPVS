# SPVS: an interpretable spatial model of plaque vulnerability

This repository provides the core code used to derive, apply and evaluate the SPVS model. The development scripts have been consolidated into a small set of function-oriented R files. The pre-built `SPVS_FINAL.tar.gz` package can be installed directly for spatial boundary prediction.

![SPVS workflow](assets/SPVS_workflow_overview.jpeg)

## Install the SPVS package

```r
install.packages("YourPath/SPVS_FINAL.tar.gz", repos = NULL, type = "source")
library(SPVS)
```

## Read spatial transcriptomics data

### Seurat object

```r
SpatialObject <- readr::read_rds("YourPath/SpatialObject.rds.gz")
```

### H5AD file

```r
Sys.setenv(SPVS_PYTHON = "YourPath/python")

SpatialData <- spvs_read_spatial(
  input = "YourPath/SpatialData.h5ad",
  input_type = "h5ad",
  sample_id = "Sample01"
)
```

### Expression matrix and coordinates

```r
ExpressionMatrix <- readr::read_csv("YourPath/ExpressionMatrix.csv")
SpatialCoordinates <- readr::read_csv("YourPath/SpatialCoordinates.csv")

SpatialData <- spvs_from_matrix(
  counts = as.matrix(ExpressionMatrix[, -1]),
  coords = SpatialCoordinates,
  sample_id = "Sample01"
)
```

The count matrix must contain genes in rows and spatial spots in columns. The coordinate table must contain a spot identifier and spatial `x` and `y` coordinates.

## Run the complete spatial workflow

```r
SPVSResult <- run_spvs_on_spatial(
  input = "YourPath/SpatialData.h5ad",
  input_type = "h5ad",
  sample_id = "Sample01",
  out_dir = "YourPath/SPVS_results"
)
```

The workflow applies the fixed SPVS molecular programs, infers plaque-core probability, delineates the plaque vulnerability interface (PVI), summarizes boundary-associated cell states and functional programs, and exports tables and publication-ready PDF figures.

## Core source files

| File | Purpose |
|---|---|
| `R/ModelDerivation.R` | Cross-cohort evidence integration, stable model-size selection and construction of the frozen SPVS model |
| `R/SpatialDataInput.R` | Input adapters for H5AD, Seurat, Visium and matrix-based spatial data |
| `R/PlaqueCoreModel.R` | Molecular-program activity and plaque-core probability |
| `R/PlaqueVulnerabilityInterface.R` | Spatial graph construction and PVI delineation |
| `R/SingleCellGrounding.R` | Projection of the frozen model to single-cell reference data |
| `R/BoundaryEcology.R` | Boundary-associated cell-state and functional-program analysis |
| `R/BoundaryValidation.R` | Spot-level and distance-aware boundary evaluation |
| `R/Visualization.R` | Spatial maps and summary figures |
| `R/SPVSWorkflow.R` | End-to-end analysis wrapper |

Minimal runnable examples are provided in `examples/`.
