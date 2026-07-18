# Matrix correction benchmark
#
# Input:  data/01_matrix_correction/benchmark_metrics.tsv.gz
# Output: analysis/01_matrix_correction/
# Panels: Extended Data Figure 4

helper <- if (file.exists("scripts/00_functions.R")) {
  "scripts/00_functions.R"
} else {
  "../scripts/00_functions.R"
}
source(helper)

root <- axis_repo_root()
out_dir <- analysis_dir(root, "01_matrix_correction")
benchmark <- read_tsv(file.path(
  root, "data", "01_matrix_correction", "benchmark_metrics.tsv.gz"
))

require_columns(
  benchmark,
  c("method", "batch_score", "biology_score", "reported_composite")
)

benchmark$composite_score <-
  0.4 * benchmark$batch_score + 0.6 * benchmark$biology_score
benchmark$composite_rank <- rank(-benchmark$composite_score, ties.method = "min")
benchmark <- benchmark[order(benchmark$composite_rank), ]

if (max(abs(benchmark$composite_score - benchmark$reported_composite)) > 1e-12) {
  stop("The recomputed composite score differs from the supplied benchmark table.")
}

write_tsv(benchmark, file.path(out_dir, "correction_method_scores.tsv"))

with_pdf(file.path(out_dir, "correction_benchmark.pdf"), width = 9, height = 5.5, {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.25, 1))

  ordered <- benchmark[order(benchmark$composite_score), ]
  cols <- ifelse(ordered$method_id == "limma", axis_colours["red"], axis_colours["blue"])
  par(mar = c(4, 8, 2, 1))
  barplot(
    ordered$composite_score,
    names.arg = ordered$method,
    horiz = TRUE,
    las = 1,
    col = cols,
    border = NA,
    xlim = c(0, 0.75),
    xlab = "Composite score"
  )
  title("Correction benchmark", adj = 0)

  par(mar = c(4, 1, 2, 1))
  y <- rev(seq_len(nrow(benchmark)))
  plot(
    c(benchmark$batch_score, benchmark$biology_score),
    rep(y, 2),
    type = "n",
    yaxt = "n",
    xlab = "Component score",
    ylab = "",
    xlim = c(0, 1)
  )
  segments(benchmark$batch_score, y, benchmark$biology_score, y, col = "#D9D9D9")
  points(benchmark$batch_score, y, pch = 16, col = axis_colours["navy"])
  points(benchmark$biology_score, y, pch = 16, col = axis_colours["gold"])
  legend(
    "bottomright",
    legend = c("Batch correction", "Biology preservation"),
    pch = 16,
    col = c(axis_colours["navy"], axis_colours["gold"]),
    bty = "n",
    cex = 0.8
  )
  title("Component scores", adj = 0)
})

message_done("Matrix correction benchmark", out_dir)
