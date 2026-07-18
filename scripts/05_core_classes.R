# Core protein classes
#
# Input:  data/05_core_classes/
# Output: analysis/05_core_classes/
# Panels: Figure 2j-k and Extended Data Figure 2

helper <- if (file.exists("scripts/00_functions.R")) {
  "scripts/00_functions.R"
} else {
  "../scripts/00_functions.R"
}
source(helper)

root <- axis_repo_root()
out_dir <- analysis_dir(root, "05_core_classes")

features <- read_tsv(file.path(root, "data", "05_core_classes", "core_features.tsv.gz"))
assignments <- read_tsv(file.path(root, "data", "05_core_classes", "class_assignments.tsv"))
feature_meta <- read_tsv(file.path(root, "data", "05_core_classes", "feature_metadata.tsv"))

require_columns(features, "protein_id")
require_columns(assignments, c("protein_id", "core_class"))
if (anyDuplicated(assignments$protein_id)) stop("Core-class assignments are not unique.")

features <- features[match(assignments$protein_id, features$protein_id), ]
if (any(is.na(features$protein_id))) stop("A class-assigned protein is absent from the feature table.")

gower_distance <- function(x, metadata) {
  n <- nrow(x)
  numerator <- matrix(0, n, n)
  denominator <- matrix(0, n, n)

  for (feature in metadata$feature) {
    value <- x[[feature]]
    if (metadata$clustering_encoding[metadata$feature == feature][1] == "continuous") {
      value <- suppressWarnings(as.numeric(value))
      observed <- is.finite(value)
      feature_range <- diff(range(value[observed], na.rm = TRUE))
      if (!is.finite(feature_range) || feature_range == 0) next
      scaled <- value / feature_range
      difference <- abs(outer(scaled, scaled, "-"))
    } else {
      value <- as.character(value)
      observed <- !is.na(value) & nzchar(value)
      difference <- outer(value, value, "!=") * 1
    }
    valid <- outer(observed, observed, "&")
    difference[!valid] <- 0
    numerator <- numerator + difference
    denominator <- denominator + valid
  }

  distance <- numerator / denominator
  distance[!is.finite(distance)] <- 1
  diag(distance) <- 0
  stats::as.dist(distance)
}

distance <- gower_distance(features, feature_meta)
coordinates <- cmdscale(distance, k = 2, add = TRUE, eig = TRUE)
pcoa <- data.frame(
  protein_id = assignments$protein_id,
  core_class = assignments$core_class,
  PCo1 = coordinates$points[, 1],
  PCo2 = coordinates$points[, 2]
)

class_counts <- as.data.frame(table(assignments$core_class), stringsAsFactors = FALSE)
names(class_counts) <- c("core_class", "n_proteins")
class_counts <- class_counts[order(class_counts$core_class), ]

numeric_features <- feature_meta$feature[
  feature_meta$clustering_encoding == "continuous"
]
numeric_matrix <- sapply(numeric_features, function(feature) {
  value <- suppressWarnings(as.numeric(features[[feature]]))
  value[!is.finite(value)] <- median(value[is.finite(value)], na.rm = TRUE)
  as.numeric(scale(value))
})
profiles <- aggregate(numeric_matrix, list(core_class = assignments$core_class), mean)

write_tsv(class_counts, file.path(out_dir, "core_class_counts.tsv"))
write_tsv(pcoa, file.path(out_dir, "core_class_pcoa.tsv"))
write_tsv(profiles, file.path(out_dir, "core_class_feature_profiles.tsv"))

class_colours <- setNames(
  c("#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1", "#76B7B2"),
  paste0("CC", 1:6)
)

with_pdf(file.path(out_dir, "core_classes.pdf"), width = 9, height = 5.5, {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  layout(matrix(c(1, 2), nrow = 1), widths = c(1.2, 0.8))

  par(mar = c(4, 4, 2, 1))
  plot(
    pcoa$PCo1,
    pcoa$PCo2,
    pch = 21,
    bg = class_colours[pcoa$core_class],
    col = "white",
    cex = 0.8,
    xlab = "Gower PCoA 1",
    ylab = "Gower PCoA 2"
  )
  legend(
    "topright",
    legend = names(class_colours),
    pt.bg = class_colours,
    pch = 21,
    bty = "n",
    cex = 0.8
  )
  title("Core protein classes", adj = 0)

  par(mar = c(4, 4, 2, 1))
  barplot(
    class_counts$n_proteins,
    names.arg = class_counts$core_class,
    col = class_colours[class_counts$core_class],
    border = NA,
    ylab = "Proteins"
  )
  title("Class sizes", adj = 0)
})

message_done("Core-class analysis", out_dir)
