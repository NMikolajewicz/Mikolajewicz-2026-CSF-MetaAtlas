# Run all publication analyses in order.

helper <- if (file.exists("scripts/00_functions.R")) {
  "scripts/00_functions.R"
} else {
  "00_functions.R"
}
source(helper)
root <- axis_repo_root()

scripts <- c(
  "01_matrix_correction.R",
  "02_assayable_universe.R",
  "03_core_proteome.R",
  "04_detection_model.R",
  "05_core_classes.R",
  "06_meta_programs.R",
  "07_injury_homeostasis_axis.R",
  "08_enrichment.R",
  "09_axis_prediction.R"
)

for (script in scripts) {
  message("\nRunning ", script)
  sys.source(
    file.path(root, "scripts", script),
    envir = new.env(parent = globalenv())
  )
}

message("\nVerifying headline results")
sys.source(
  file.path(root, "scripts", "verify_results.R"),
  envir = new.env(parent = globalenv())
)

message("\nAll analyses completed and verified.")
