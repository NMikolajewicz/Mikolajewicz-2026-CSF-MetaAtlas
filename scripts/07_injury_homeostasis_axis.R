# Injury-homeostasis axis
#
# Input:  data/07_injury_homeostasis_axis/
# Output: analysis/07_injury_homeostasis_axis/
# Panels: Figure 4a-c,f and Extended Data Figures 8-10

helper <- if (file.exists("scripts/00_functions.R")) {
  "scripts/00_functions.R"
} else {
  "../scripts/00_functions.R"
}
source(helper)

root <- axis_repo_root()
out_dir <- analysis_dir(root, "07_injury_homeostasis_axis")

corr_long <- read_tsv(file.path(
  root, "data", "07_injury_homeostasis_axis", "mp_correlations.tsv"
))
pca_var <- read_tsv(file.path(
  root, "data", "07_injury_homeostasis_axis", "pca_variance.tsv"
))
loadings <- read_tsv(file.path(
  root, "data", "07_injury_homeostasis_axis", "pc1_loadings.tsv"
))
permutation <- read_tsv(file.path(
  root, "data", "07_injury_homeostasis_axis", "permutation_test.tsv"
))
protein_axis <- read_tsv(file.path(
  root, "data", "07_injury_homeostasis_axis", "protein_axis_correlations.tsv.gz"
))
sensitivity <- read_tsv(file.path(
  root, "data", "07_injury_homeostasis_axis", "method_sensitivity.tsv.gz"
))

mps <- sprintf("MP%02d", 1:13)
corr_matrix <- matrix(NA_real_, 13, 13, dimnames = list(mps, mps))
row_index <- match(corr_long$mp_row, mps)
col_index <- match(corr_long$mp_col, mps)
corr_matrix[cbind(row_index, col_index)] <- corr_long$pooled_rho
if (any(!is.finite(corr_matrix))) stop("The meta-program correlation matrix is incomplete.")

cluster_order <- hclust(as.dist(1 - corr_matrix), method = "ward.D2")$order
ordered_mps <- mps[cluster_order]

protein_axis$axis_set <- "not_selected"
protein_axis$axis_set[
  protein_axis$pooled_rho >= 0.4 & protein_axis$fdr_value < 0.05
] <- "injury"
protein_axis$axis_set[
  protein_axis$pooled_rho <= -0.4 & protein_axis$fdr_value < 0.05
] <- "homeostasis"
signature <- protein_axis[protein_axis$axis_set != "not_selected", ]
signature <- signature[order(signature$axis_set, -abs(signature$pooled_rho)), ]

wide_sensitivity <- merge(
  sensitivity[sensitivity$method == "pearson", c("protein_id", "pooled_rho")],
  sensitivity[sensitivity$method == "spearman", c("protein_id", "pooled_rho")],
  by = "protein_id",
  suffixes = c("_pearson", "_spearman")
)
sensitivity_congruence <- cor(
  wide_sensitivity$pooled_rho_pearson,
  wide_sensitivity$pooled_rho_spearman,
  use = "complete.obs"
)

marker_value <- function(marker) {
  value <- protein_axis$pooled_rho[protein_axis$protein_id == marker]
  if (length(value) != 1L) stop("Expected one row for ", marker)
  value
}

axis_summary <- data.frame(
  statistic = c(
    "PC1_variance_fraction", "permutation_p_value", "permutations",
    "injury_proteins", "homeostasis_proteins", "C5_spearman",
    "NRCAM_spearman", "pearson_spearman_congruence"
  ),
  value = c(
    pca_var$variance_explained[pca_var$PC == "PC1"],
    permutation$p_value,
    permutation$n_permutations,
    sum(signature$axis_set == "injury"),
    sum(signature$axis_set == "homeostasis"),
    marker_value("C5"),
    marker_value("NRCAM"),
    sensitivity_congruence
  )
)

write_tsv(axis_summary, file.path(out_dir, "axis_summary.tsv"))
write_tsv(signature, file.path(out_dir, "axis_protein_sets.tsv"))
write_tsv(loadings, file.path(out_dir, "pc1_loadings.tsv"))

with_pdf(file.path(out_dir, "injury_homeostasis_axis.pdf"), width = 12, height = 5.5, {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1.1, 0.9, 1.1))

  par(mar = c(7, 7, 2, 1))
  ordered_corr <- corr_matrix[ordered_mps, ordered_mps]
  image(
    seq_len(13),
    seq_len(13),
    t(ordered_corr[13:1, ]),
    col = grDevices::colorRampPalette(c("#2B6CB0", "white", "#C53030"))(101),
    zlim = c(-1, 1),
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  axis(1, at = seq_len(13), labels = ordered_mps, las = 2, cex.axis = 0.65)
  axis(2, at = seq_len(13), labels = rev(ordered_mps), las = 2, cex.axis = 0.65)
  box()
  title("Meta-program correlations", adj = 0)

  par(mar = c(7, 4, 2, 1))
  loading_order <- order(loadings$loading_observed)
  barplot(
    loadings$loading_observed[loading_order],
    names.arg = loadings$meta_program[loading_order],
    las = 2,
    col = ifelse(
      loadings$loading_observed[loading_order] > 0,
      axis_colours["red"],
      axis_colours["navy"]
    ),
    border = NA,
    ylab = "PC1 loading"
  )
  abline(h = 0, col = axis_colours["dark_grey"])
  title("Injury-homeostasis loadings", adj = 0)

  par(mar = c(4, 4, 2, 1))
  point_colours <- c(
    injury = axis_colours["red"],
    homeostasis = axis_colours["navy"],
    not_selected = axis_colours["grey"]
  )
  plot(
    protein_axis$pooled_rho,
    -log10(pmax(protein_axis$fdr_value, .Machine$double.xmin)),
    pch = 16,
    cex = 0.55,
    col = grDevices::adjustcolor(point_colours[protein_axis$axis_set], alpha.f = 0.7),
    xlab = "Pooled protein-axis correlation",
    ylab = expression(-log[10](FDR))
  )
  abline(v = c(-0.4, 0.4), lty = 2, col = axis_colours["dark_grey"])
  abline(h = -log10(0.05), lty = 2, col = axis_colours["dark_grey"])
  markers <- protein_axis$protein_id %in% c("C5", "NRCAM")
  marker_rows <- protein_axis[markers, ]
  text(
    marker_rows$pooled_rho,
    -log10(marker_rows$fdr_value),
    labels = marker_rows$protein_id,
    pos = ifelse(marker_rows$pooled_rho > 0, 2, 4),
    cex = 0.8,
    offset = 0.4
  )
  title("Axis-associated proteins", adj = 0)
})

message_done("Injury-homeostasis axis analysis", out_dir)
