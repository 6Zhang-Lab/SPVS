
spvs_default_core_programs <- function() {
  list(
    lipid_necrotic = c("APOE", "LPL", "TREM2", "LGALS3", "SPP1"),
    proteolytic = c("MMP9", "MMP12", "CTSK", "CTSL", "MMP14"),
    inflammatory_core = c("IL1B", "TNF", "CCL2", "CXCL8", "NFKBIA"),
    hemorrhage_oxidative = c("HMOX1", "FTL", "FTH1", "CD163")
  )
}

spvs_default_cell_programs <- function() {
  list(
    "Fibroblast / stromal" = c("COL1A1", "COL3A1", "DCN", "LUM", "FN1", "COL8A1", "MMP2", "MMP14", "PI16", "DPT"),
    "Myeloid / macrophage" = c("APOE", "LPL", "TREM2", "LGALS3", "SPP1", "MMP9", "CTSB", "LST1", "FCN1"),
    "SMC / VSMC" = c("ACTA2", "TAGLN", "MYH11", "MYL9", "CNN1", "TPM2"),
    "Endothelial" = c("PECAM1", "VWF", "KDR", "CLDN5", "ENG"),
    "T / NK" = c("CD3D", "CD3E", "TRAC", "NKG7", "GNLY"),
    "B / plasma" = c("MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN"),
    "Other immune" = c("HLA-DRA", "FCER1A", "CLEC10A", "TPSAB1", "KIT")
  )
}

spvs_default_function_programs <- function() {
  list(
    inflammation = c("IL1B", "TNF", "CCL2", "CXCL8", "NFKBIA", "SOCS3"),
    lipid_foamy_macrophage = c("APOE", "LPL", "TREM2", "LGALS3", "SPP1"),
    matrix_degradation = c("MMP9", "MMP12", "CTSK", "CTSL", "MMP14"),
    stromal_matrix_remodeling = c("COL1A1", "COL3A1", "DCN", "LUM", "FN1", "COL8A1", "MMP2"),
    SMC_inflammatory_transition = c("ACTA2", "TAGLN", "MYL9", "VCAM1", "CCL2", "IL6"),
    hemorrhage_oxidative = c("HMOX1", "FTL", "FTH1", "CD163")
  )
}
