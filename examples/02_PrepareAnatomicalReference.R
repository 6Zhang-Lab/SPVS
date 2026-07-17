source("config/paths.R")
source("R/Utilities.R")
source("R/AnatomicalReference.R")

# Prefer the study's prepared morphology-anchored interface file. If it is not
# present, reconstruct the reference shell from the standardized anatomical
# compartment table. Neither file contains SPVS prediction as an input.
ReferenceInput <- if (file.exists(REFERENCE_BOUNDARY_FILE)) {
  REFERENCE_BOUNDARY_FILE
} else {
  ANATOMICAL_REGION_FILE
}

if (!file.exists(ReferenceInput)) {
  stop("No anatomical reference input was found.", call. = FALSE)
}

AnatomicalReference <- prepare_anatomical_reference(
  input = ReferenceInput,
  shell_multiplier = 1.65
)

readr::write_csv(AnatomicalReference, PREPARED_REFERENCE_FILE)
message("Saved independent anatomical reference: ", PREPARED_REFERENCE_FILE)
