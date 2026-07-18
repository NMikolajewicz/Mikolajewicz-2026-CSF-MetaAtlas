# Consensus proteomic meta-programs
#
# Input:  data/06_meta_programs/
# Output: analysis/06_meta_programs/
# Panels: Figure 3 and Extended Data Figures 5-6

helper <- if (file.exists("scripts/00_functions.R")) {
  "scripts/00_functions.R"
} else {
  "../scripts/00_functions.R"
}
source(helper)

root <- axis_repo_root()
out_dir <- analysis_dir(root, "06_meta_programs")

funnel <- read_tsv(file.path(root, "data", "06_meta_programs", "filter_progression.tsv"))
membership <- read_tsv(file.path(
  root, "data", "06_meta_programs", "component_membership.tsv.gz"
))
programs <- read_tsv(file.path(
  root, "data", "06_meta_programs", "meta_program_proteins.tsv.gz"
))
loso <- read_tsv(file.path(root, "data", "06_meta_programs", "loso_results.tsv.gz"))
pca_ceiling <- read_tsv(file.path(root, "data", "06_meta_programs", "pca_ceiling.tsv"))
null_benchmark <- read_tsv(file.path(root, "data", "06_meta_programs", "null_benchmark.tsv"))

funnel_full <- rbind(
  funnel[, c("stage", "n_programs")],
  data.frame(
    stage = c("assigned_components", "meta_programs"),
    n_programs = c(nrow(membership), length(unique(membership$meta_program)))
  )
)

loso_by_program <- do.call(rbind, lapply(split(loso, loso$meta_program), function(x) {
  data.frame(
    meta_program = x$meta_program[1],
    n_studies = nrow(x),
    median_spearman = median(x$spearman_rho, na.rm = TRUE),
    mean_spearman = mean(x$spearman_rho, na.rm = TRUE),
    recovered_fraction = mean(x$recovered, na.rm = TRUE)
  )
}))
row.names(loso_by_program) <- NULL
loso_by_program <- loso_by_program[order(loso_by_program$meta_program), ]

program_summary <- data.frame(
  statistic = c(
    "candidate_components", "nonredundant_components", "assigned_components",
    "meta_programs", "unique_meta_program_proteins", "median_program_spearman",
    "median_program_recovery"
  ),
  value = c(
    funnel_full$n_programs[funnel_full$stage == "all_candidates"],
    funnel_full$n_programs[funnel_full$stage == "nonredundant"],
    nrow(membership),
    length(unique(membership$meta_program)),
    length(unique(programs$gene)),
    median(loso_by_program$median_spearman),
    median(loso_by_program$recovered_fraction)
  )
)

write_tsv(funnel_full, file.path(out_dir, "program_filtering_funnel.tsv"))
write_tsv(loso_by_program, file.path(out_dir, "loso_by_meta_program.tsv"))
write_tsv(program_summary, file.path(out_dir, "meta_program_summary.tsv"))
write_tsv(programs, file.path(out_dir, "meta_program_proteins.tsv"))

with_pdf(file.path(out_dir, "meta_programs.pdf"), width = 11, height = 5.5, {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(0.9, 1.2, 1))

  par(mar = c(8, 4, 2, 1))
  barplot(
    funnel_full$n_programs,
    names.arg = gsub("_", " ", funnel_full$stage),
    las = 2,
    log = "y",
    col = axis_colours["blue"],
    border = NA,
    ylab = "Programs or components"
  )
  title("Program filtering", adj = 0)

  par(mar = c(7, 4, 2, 1))
  loso$meta_program <- factor(loso$meta_program, levels = sprintf("MP%02d", 1:13))
  boxplot(
    spearman_rho ~ meta_program,
    data = loso,
    las = 2,
    outline = FALSE,
    col = axis_colours["teal"],
    border = "white",
    ylab = "Held-out Spearman correlation"
  )
  abline(h = 0, lty = 2, col = axis_colours["grey"])
  title("Leave-one-study-out transfer", adj = 0)

  par(mar = c(4, 4, 2, 1))
  plot(
    pca_ceiling$fve_pca_k_insample,
    pca_ceiling$fve_pca_k_cv,
    pch = 21,
    bg = axis_colours["gold"],
    col = "white",
    cex = 1.1,
    xlab = "In-sample PCA ceiling",
    ylab = "Cross-validated PCA ceiling",
    xlim = c(0, 1),
    ylim = c(0, 1)
  )
  abline(0, 1, lty = 2, col = axis_colours["grey"])
  title("Study-wise variance ceiling", adj = 0)
})

message_done("Meta-program analysis", out_dir)
