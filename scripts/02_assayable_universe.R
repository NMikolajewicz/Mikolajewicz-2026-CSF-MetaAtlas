# Rolling local-support assayable universe
#
# Input:  data/02_assayable_universe/
# Output: analysis/02_assayable_universe/
# Panels: Extended Data Figure 1a-b

helper <- if (file.exists("scripts/00_functions.R")) {
  "scripts/00_functions.R"
} else {
  "../scripts/00_functions.R"
}
source(helper)

root <- axis_repo_root()
out_dir <- analysis_dir(root, "02_assayable_universe")

windows <- read_tsv(file.path(
  root, "data", "02_assayable_universe", "study_windows.tsv"
))
assay <- read_tsv(file.path(
  root, "data", "02_assayable_universe", "protein_study_assayability.tsv.gz"
))

require_columns(
  assay,
  c(
    "protein_id", "study_id", "study_depth_n", "final_assayable",
    "forced_into_universe", "peer_support_n", "window_size", "n_window_peers"
  )
)

unsupported <- assay$final_assayable &
  !(assay$peer_support_n >= 2 | assay$forced_into_universe)
if (any(unsupported, na.rm = TRUE)) {
  stop("A final assayability call lacks two-peer support or held-out detection.")
}

split_assay <- split(assay, assay$study_id)
study_summary <- do.call(rbind, lapply(split_assay, function(x) {
  data.frame(
    study_id = x$study_id[1],
    study_depth_n = x$study_depth_n[1],
    proteins_with_two_peer_support = sum(x$peer_support_n >= 2, na.rm = TRUE),
    heldout_detected_additions = sum(x$forced_into_universe, na.rm = TRUE),
    final_assayable_n = sum(x$final_assayable, na.rm = TRUE)
  )
}))
row.names(study_summary) <- NULL

study_summary <- merge(
  windows,
  study_summary,
  by = c("study_id", "study_depth_n"),
  all.x = TRUE,
  sort = FALSE
)
study_summary <- study_summary[order(study_summary$study_rank), ]

if (any(study_summary$final_universe_n != study_summary$final_assayable_n)) {
  stop("The study-level universe counts do not agree with the protein-level flags.")
}

write_tsv(study_summary, file.path(out_dir, "study_assayable_universes.tsv"))

with_pdf(file.path(out_dir, "assayable_universe.pdf"), width = 9, height = 5.5, {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  par(mar = c(7, 4, 2, 1))
  x <- seq_len(nrow(study_summary))
  plot(
    x,
    study_summary$final_universe_n,
    type = "o",
    pch = 16,
    col = axis_colours["red"],
    ylim = range(c(study_summary$study_depth_n, study_summary$final_universe_n)),
    xaxt = "n",
    xlab = "",
    ylab = "Number of proteins"
  )
  lines(x, study_summary$study_depth_n, type = "o", pch = 1, col = axis_colours["navy"])
  axis(1, at = x, labels = study_summary$study_id, las = 2, cex.axis = 0.65)
  legend(
    "topleft",
    legend = c("Rolling local-support universe", "Observed study depth"),
    col = c(axis_colours["red"], axis_colours["navy"]),
    pch = c(16, 1),
    lty = 1,
    bty = "n"
  )
  title("Coverage-aware assayable universe", adj = 0)
})

message_done("Assayable-universe analysis", out_dir)
