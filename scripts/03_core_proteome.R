# Core CSF proteome
#
# Input:  data/03_core_proteome/protein_detection_meta_analysis.tsv.gz
# Output: analysis/03_core_proteome/
# Panels: Figure 2a-e and associated Extended Data panels

helper <- if (file.exists("scripts/00_functions.R")) {
  "scripts/00_functions.R"
} else {
  "../scripts/00_functions.R"
}
source(helper)

root <- axis_repo_root()
out_dir <- analysis_dir(root, "03_core_proteome")
core <- read_tsv(file.path(
  root, "data", "03_core_proteome", "protein_detection_meta_analysis.tsv.gz"
))

require_columns(
  core,
  c(
    "protein_id", "coverage_aware_detection", "coverage_aware_lower95",
    "coverage_naive_detection", "n_studies", "loso_class_demoted", "tier"
  )
)

tier_levels <- c("A", "B", "C", "D", "E")
core$tier <- factor(core$tier, levels = tier_levels)
tier_counts <- as.data.frame(table(core$tier), stringsAsFactors = FALSE)
names(tier_counts) <- c("tier", "n_proteins")
tier_counts$core <- tier_counts$tier %in% c("A", "B")

core_proteins <- core[core$tier %in% c("A", "B"), ]
core_proteins <- core_proteins[order(core_proteins$tier, -core_proteins$coverage_aware_detection), ]

write_tsv(tier_counts, file.path(out_dir, "core_tier_counts.tsv"))
write_tsv(core_proteins, file.path(out_dir, "core_proteins.tsv"))

with_pdf(file.path(out_dir, "core_proteome.pdf"), width = 9, height = 5.5, {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  layout(matrix(c(1, 2), nrow = 1), widths = c(0.8, 1.3))

  par(mar = c(4, 4, 2, 1))
  bar_cols <- c(axis_colours["red"], axis_colours["gold"], axis_colours["teal"],
                axis_colours["blue"], axis_colours["grey"])
  barplot(
    tier_counts$n_proteins,
    names.arg = paste0("Tier ", tier_counts$tier),
    col = bar_cols,
    border = NA,
    ylab = "Proteins",
    log = "y"
  )
  title("Core-proteome tiers", adj = 0)

  par(mar = c(4, 4, 2, 1))
  point_cols <- bar_cols[as.integer(core$tier)]
  plot(
    core$coverage_naive_detection,
    core$coverage_aware_detection,
    pch = 16,
    cex = 0.45,
    col = grDevices::adjustcolor(point_cols, alpha.f = 0.45),
    xlab = "Coverage-naive detection",
    ylab = "Coverage-aware detection"
  )
  abline(0, 1, lty = 2, col = axis_colours["dark_grey"])
  title("Effect of assayability correction", adj = 0)
})

message(
  "Core proteins (Tier A + B): ", nrow(core_proteins),
  "; proteins evaluated: ", nrow(core)
)
message_done("Core-proteome analysis", out_dir)
