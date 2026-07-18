# Protein-detection model
#
# Input:  data/04_detection_model/
# Output: analysis/04_detection_model/
# Panels: Figure 2f-i

helper <- if (file.exists("scripts/00_functions.R")) {
  "scripts/00_functions.R"
} else {
  "../scripts/00_functions.R"
}
source(helper)

root <- axis_repo_root()
out_dir <- analysis_dir(root, "04_detection_model")

pred <- read_tsv(file.path(
  root, "data", "04_detection_model", "cross_validated_predictions.tsv.gz"
))
coef <- read_tsv(file.path(
  root, "data", "04_detection_model", "model_coefficients.tsv.gz"
))

require_columns(
  pred,
  c(
    "observed", "core_protein", "predicted_source_head", "predicted_flat",
    "predicted_mean"
  )
)

model_columns <- c(
  source_head = "predicted_source_head",
  flat_elastic_net = "predicted_flat",
  mean_baseline = "predicted_mean"
)

performance <- do.call(rbind, lapply(names(model_columns), function(model) {
  score <- pred[[model_columns[[model]]]]
  spearman <- if (stats::sd(score, na.rm = TRUE) == 0) {
    NA_real_
  } else {
    cor(pred$observed, score, method = "spearman", use = "complete.obs")
  }
  data.frame(
    model = model,
    n = sum(is.finite(score) & is.finite(pred$observed)),
    spearman_rho = spearman,
    rmse = sqrt(mean((pred$observed - score)^2, na.rm = TRUE)),
    auroc = roc_auc(pred$core_protein, score),
    auprc = average_precision(pred$core_protein, score)
  )
}))

write_tsv(performance, file.path(out_dir, "detection_model_performance.tsv"))

top_coef <- coef[order(-coef$abs_importance), ]
top_coef <- top_coef[!duplicated(top_coef$feature), ]
top_coef <- head(top_coef, 25)
write_tsv(top_coef, file.path(out_dir, "top_model_coefficients.tsv"))

with_pdf(file.path(out_dir, "detection_model.pdf"), width = 9, height = 5.5, {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  layout(matrix(c(1, 2), nrow = 1))

  par(mar = c(4, 4, 2, 1))
  plot(
    pred$observed,
    pred$predicted_source_head,
    pch = 16,
    cex = 0.55,
    col = grDevices::adjustcolor(axis_colours["blue"], alpha.f = 0.35),
    xlab = "Observed cross-study detection",
    ylab = "Cross-validated prediction",
    xlim = c(0, 1),
    ylim = c(0, 1)
  )
  abline(0, 1, lty = 2, col = axis_colours["dark_grey"])
  title("Source-head elastic net", adj = 0)

  curve <- roc_points(pred$core_protein, pred$predicted_source_head)
  par(mar = c(4, 4, 2, 1))
  plot(
    curve$false_positive_rate,
    curve$true_positive_rate,
    type = "l",
    lwd = 2.5,
    col = axis_colours["red"],
    xlab = "False-positive rate",
    ylab = "True-positive rate",
    xlim = c(0, 1),
    ylim = c(0, 1),
    asp = 1
  )
  abline(0, 1, lty = 2, col = axis_colours["grey"])
  legend(
    "bottomright",
    legend = sprintf(
      "AUROC %.3f\nAUPRC %.3f",
      performance$auroc[performance$model == "source_head"],
      performance$auprc[performance$model == "source_head"]
    ),
    bty = "n"
  )
  title("Tier A+B discrimination", adj = 0)
})

message_done("Detection-model analysis", out_dir)
