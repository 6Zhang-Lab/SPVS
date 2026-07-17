# Apply the frozen SPVS model to a single-cell reference without gene reselection.
#
# Example input:
# SingleCell <- readr::read_rds("YourPath/SingleCellReference.rds.gz")
# Model <- readr::read_csv("YourPath/SPVSModel.csv")
# SingleCell <- project_spvs_to_seurat(SingleCell, Model)

#' Project SPVS model activity to a Seurat object
project_spvs_to_seurat <- function(object, model, assay = "RNA",
                                   layer = "data", name = "SPVS_activity") {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("Package 'SeuratObject' is required.", call. = FALSE)
  }
  expression <- tryCatch(
    SeuratObject::LayerData(object, assay = assay, layer = layer),
    error = function(e) SeuratObject::GetAssayData(object, assay = assay, slot = layer)
  )
  idx <- match(toupper(model$gene), toupper(rownames(expression)))
  present <- unique(idx[!is.na(idx)])
  if (!length(present)) stop("No model genes were found in the Seurat object.", call. = FALSE)
  activity <- calculate_spvs_activity(as.matrix(expression[present, , drop = FALSE]), model)
  names(activity) <- colnames(expression)
  object[[name]] <- activity[colnames(object)]
  object
}
