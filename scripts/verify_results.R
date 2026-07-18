# Check the principal numerical results produced by scripts/run_all.R.

helper <- if (file.exists("scripts/00_functions.R")) {
  "scripts/00_functions.R"
} else {
  "00_functions.R"
}
source(helper)
root <- axis_repo_root()

read_result <- function(directory, filename) {
  path <- file.path(root, "analysis", directory, filename)
  if (!file.exists(path)) {
    stop("Result file not found: ", path, ". Run scripts/run_all.R first.")
  }
  read_tsv(path)
}

expect_equal <- function(observed, expected, label, tolerance = 1e-10) {
  same <- isTRUE(all.equal(
    unname(observed),
    unname(expected),
    tolerance = tolerance,
    check.attributes = FALSE
  ))
  if (!same) {
    stop(
      label,
      " differs from the release value. Observed: ",
      paste(observed, collapse = ", "),
      "; expected: ",
      paste(expected, collapse = ", ")
    )
  }
}

correction <- read_result(
  "01_matrix_correction",
  "correction_method_scores.tsv"
)
limma <- correction[correction$method_id == "limma", , drop = FALSE]
expect_equal(nrow(limma), 1L, "Number of limma benchmark rows")
expect_equal(limma$composite_score, 0.674434012339625, "Limma composite score")
expect_equal(limma$composite_rank, 1L, "Limma benchmark rank")

tiers <- read_result("03_core_proteome", "core_tier_counts.tsv")
expect_equal(tiers$tier, c("A", "B", "C", "D", "E"), "Core tier labels")
expect_equal(
  tiers$n_proteins,
  c(337L, 410L, 497L, 1518L, 9844L),
  "Core tier counts"
)
core <- read_result("03_core_proteome", "core_proteins.tsv")
expect_equal(nrow(core), 747L, "Core-proteome size")

detection <- read_result(
  "04_detection_model",
  "detection_model_performance.tsv"
)
source_head <- detection[detection$model == "source_head", , drop = FALSE]
expect_equal(source_head$n, 2762L, "Detection-model protein count")
expect_equal(source_head$spearman_rho, 0.540232352131842, "Detection-model rho")
expect_equal(source_head$auroc, 0.832664653651828, "Detection-model AUROC")
expect_equal(source_head$auprc, 0.690502356175033, "Detection-model AUPRC")

classes <- read_result("05_core_classes", "core_class_counts.tsv")
expect_equal(
  classes$core_class,
  paste0("CC", seq_len(6L)),
  "Core Class labels"
)
expect_equal(
  classes$n_proteins,
  c(190L, 185L, 92L, 193L, 42L, 45L),
  "Core Class counts"
)

funnel <- read_result("06_meta_programs", "program_filtering_funnel.tsv")
expect_equal(
  funnel$n_programs,
  c(28795L, 22271L, 10805L, 634L, 247L, 13L),
  "Program-filtering series"
)
programs <- read_result("06_meta_programs", "meta_program_proteins.tsv")
expect_equal(length(unique(programs$meta_program)), 13L, "Meta-program count")
expect_equal(length(unique(programs$gene)), 524L, "Meta-program protein count")

axis <- read_result("07_injury_homeostasis_axis", "axis_summary.tsv")
axis_value <- function(statistic) {
  value <- axis$value[axis$statistic == statistic]
  if (length(value) != 1L) stop("Expected one axis result for ", statistic)
  value
}
expect_equal(axis_value("PC1_variance_fraction"), 0.555557802027216, "Axis PC1")
expect_equal(axis_value("injury_proteins"), 71L, "Injury-associated proteins")
expect_equal(axis_value("homeostasis_proteins"), 236L, "Homeostasis-associated proteins")
expect_equal(axis_value("C5_spearman"), 0.679050495344441, "C5 axis correlation")
expect_equal(axis_value("NRCAM_spearman"), -0.706793618055807, "NRCAM axis correlation")

prediction <- read_result(
  "09_axis_prediction",
  "all_study_recomputed_meta_analysis.tsv"
)
pooled_r2 <- prediction$estimate[prediction$metric == "R2"]
expect_equal(pooled_r2, 0.812413668622287, "All-study pooled R-squared")

message("Release verification passed.")
