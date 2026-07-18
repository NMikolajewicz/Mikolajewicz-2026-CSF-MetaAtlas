# Gene-set enrichment
#
# Input:  data/08_enrichment/
# Output: analysis/08_enrichment/
# Panels: Figure 2k, Figure 3e and Figure 4d-e

helper <- if (file.exists("scripts/00_functions.R")) {
  "scripts/00_functions.R"
} else {
  "../scripts/00_functions.R"
}
source(helper)

root <- axis_repo_root()
out_dir <- analysis_dir(root, "08_enrichment")

class_enrichment <- read_tsv(file.path(
  root, "data", "08_enrichment", "core_class_enrichment.tsv"
))
mp_enrichment <- read_tsv(file.path(
  root, "data", "08_enrichment", "meta_program_enrichment.tsv.gz"
))
axis_gsea <- read_tsv(file.path(
  root, "data", "08_enrichment", "axis_gsea.tsv.gz"
))
mp_network <- read_tsv(file.path(
  root, "data", "08_enrichment", "meta_program_network.tsv"
))

class_sig <- class_enrichment[class_enrichment$fdr < 0.05, ]
mp_sig <- mp_enrichment[mp_enrichment$fdr_global < 0.05, ]
axis_sig <- axis_gsea[axis_gsea$padj < 0.05, ]

class_sig <- class_sig[order(class_sig$fdr, -abs(log2(class_sig$odds_ratio_cc))), ]
mp_sig <- mp_sig[order(mp_sig$fdr_global, -mp_sig$odds_ratio), ]
axis_sig <- axis_sig[order(axis_sig$padj, -abs(axis_sig$NES)), ]

enrichment_summary <- data.frame(
  analysis = c("core_classes", "meta_programs", "axis_ranked"),
  significant_terms = c(nrow(class_sig), nrow(mp_sig), nrow(axis_sig)),
  recorded_universe = c(
    747,
    paste(sort(unique(mp_enrichment$universe_size)), collapse = ";"),
    917
  )
)

write_tsv(class_sig, file.path(out_dir, "core_class_enrichment_significant.tsv"))
write_tsv(mp_sig, file.path(out_dir, "meta_program_enrichment_significant.tsv"))
write_tsv(axis_sig, file.path(out_dir, "axis_gsea_significant.tsv"))
write_tsv(mp_network, file.path(out_dir, "meta_program_enrichment_network.tsv"))
write_tsv(enrichment_summary, file.path(out_dir, "enrichment_summary.tsv"))

top_class <- head(class_sig, 15)
top_mp <- head(mp_sig, 15)
top_axis <- head(axis_sig[order(abs(axis_sig$NES), decreasing = TRUE), ], 15)
short_label <- function(x, width = 36) {
  x <- gsub("_", " ", x)
  ifelse(nchar(x) > width, paste0(substr(x, 1, width - 3), "..."), x)
}

with_pdf(file.path(out_dir, "gene_set_enrichment.pdf"), width = 12, height = 6, {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  layout(matrix(c(1, 2, 3), nrow = 1))

  par(mar = c(4, 8, 2, 1))
  class_values <- log2(top_class$odds_ratio_cc)
  barplot(
    rev(class_values),
    names.arg = rev(short_label(paste(top_class$core_class, top_class$class_value, sep = ": "))),
    horiz = TRUE,
    las = 1,
    col = ifelse(rev(class_values) > 0, axis_colours["red"], axis_colours["navy"]),
    border = NA,
    xlab = expression(log[2](odds~ratio))
  )
  title("Core classes", adj = 0)

  par(mar = c(4, 8, 2, 1))
  barplot(
    rev(-log10(top_mp$fdr_global)),
    names.arg = rev(short_label(paste(top_mp$meta_program, sub(".*::", "", top_mp$pathway), sep = ": "))),
    horiz = TRUE,
    las = 1,
    col = axis_colours["teal"],
    border = NA,
    xlab = expression(-log[10](global~FDR)),
    cex.names = 0.58
  )
  title("Meta-programs", adj = 0)

  par(mar = c(4, 8, 2, 1))
  barplot(
    rev(top_axis$NES),
    names.arg = rev(short_label(top_axis$pathway_label)),
    horiz = TRUE,
    las = 1,
    col = ifelse(rev(top_axis$NES) > 0, axis_colours["red"], axis_colours["navy"]),
    border = NA,
    xlab = "Normalized enrichment score",
    cex.names = 0.65
  )
  title("Injury-homeostasis axis", adj = 0)
})

message_done("Enrichment analysis", out_dir)
