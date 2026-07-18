# C5 and NRCAM prediction of the injury-homeostasis axis
#
# Input:  data/09_axis_prediction/
# Output: analysis/09_axis_prediction/
# Panels: Figure 4m-o

helper <- if (file.exists("scripts/00_functions.R")) {
  "scripts/00_functions.R"
} else {
  "../scripts/00_functions.R"
}
source(helper)

root <- axis_repo_root()
out_dir <- analysis_dir(root, "09_axis_prediction")

ladder <- read_tsv(file.path(root, "data", "09_axis_prediction", "model_ladder.tsv"))
forest <- read_tsv(file.path(
  root, "data", "09_axis_prediction", "all_study_models.tsv.gz"
))
reported_meta <- read_tsv(file.path(
  root, "data", "09_axis_prediction", "all_study_meta_analysis.tsv"
))

if (!all(forest$axis_omits_c5_nrcam, na.rm = TRUE)) {
  stop("C5 or NRCAM is present in an axis outcome used for prediction.")
}

eligible <- forest$evaluable & is.finite(forest$fisher_z) & is.finite(forest$se_z)
pooled <- random_effects_dl(forest$fisher_z[eligible], forest$se_z[eligible])
pooled_r <- tanh(pooled$estimate)
pooled_low_r <- tanh(pooled$lower95)
pooled_high_r <- tanh(pooled$upper95)

recomputed_meta <- data.frame(
  metric = "R2",
  studies = pooled$k,
  estimate = pooled_r^2,
  lower95 = pooled_low_r^2,
  upper95 = pooled_high_r^2,
  tau2_fisher_z = pooled$tau2,
  i2 = pooled$i2
)

reported_r2 <- reported_meta$estimate[reported_meta$metric_id == "r2"]
if (abs(recomputed_meta$estimate - reported_r2) > 1e-6) {
  stop("The repooled all-study R-squared does not match the supplied meta-analysis.")
}

write_tsv(ladder, file.path(out_dir, "cohort_model_ladder.tsv"))
write_tsv(recomputed_meta, file.path(out_dir, "all_study_recomputed_meta_analysis.tsv"))

with_pdf(file.path(out_dir, "axis_prediction.pdf"), width = 12, height = 6.5, {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 1, 1.2))

  bader <- ladder[ladder$panel_id == "A", ]
  bader_label <- paste(
    sub(" validation$", "", sub(" train CV$", "", bader$cohort_label)),
    bader$model_label,
    sep = " | "
  )
  par(mar = c(4, 10, 2, 1))
  barplot(
    bader$estimate,
    names.arg = bader_label,
    horiz = TRUE,
    las = 1,
    col = rep(c(axis_colours["blue"], axis_colours["teal"]), each = 5),
    border = NA,
    xlab = expression(R^2),
    xlim = c(0, 0.95),
    cex.names = 0.52
  )
  title("Bader platform analyses", adj = 0)

  riviere <- ladder[ladder$panel_id == "B", ]
  par(mar = c(4, 8, 2, 1))
  barplot(
    riviere$estimate,
    names.arg = riviere$model_label,
    horiz = TRUE,
    las = 1,
    col = ifelse(grepl("NRCAM", riviere$model_label), axis_colours["red"], axis_colours["navy"]),
    border = NA,
    xlab = expression(R^2),
    xlim = c(0, 0.95),
    cex.names = 0.58
  )
  title("Riviere-Cazaux analysis", adj = 0)

  forest_plot <- forest[eligible, ]
  forest_plot <- forest_plot[order(forest_plot$r2), ]
  y <- seq_len(nrow(forest_plot))
  par(mar = c(4, 8, 2, 1))
  plot(
    forest_plot$r2,
    y,
    pch = 16,
    col = axis_colours["blue"],
    xlim = c(0, 1),
    yaxt = "n",
    xlab = expression(Study~R^2),
    ylab = ""
  )
  segments(forest_plot$ci_low_r2, y, forest_plot$ci_high_r2, y, col = axis_colours["blue"])
  axis(2, at = y, labels = forest_plot$study_label, las = 2, cex.axis = 0.48)
  abline(v = recomputed_meta$estimate, lty = 2, col = axis_colours["red"])
  title("All-study C5 + NRCAM models", adj = 0)
})

message_done("Axis-prediction analysis", out_dir)
