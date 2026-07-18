# Shared functions used by the analysis scripts

axis_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  for (i in seq_len(8L)) {
    if (dir.exists(file.path(current, "data")) &&
        dir.exists(file.path(current, "scripts"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) break
    current <- parent
  }
  stop("Run the script from the repository root or one of its subdirectories.")
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else path
  on.exit(if (inherits(con, "connection")) close(con), add = TRUE)
  read.delim(
    con,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", "")
  )
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(
    x,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
  invisible(path)
}

analysis_dir <- function(root, name) {
  path <- file.path(root, "analysis", name)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

require_columns <- function(x, columns, label = deparse(substitute(x))) {
  missing <- setdiff(columns, names(x))
  if (length(missing)) {
    stop(label, " is missing: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

with_pdf <- function(path, width = 8, height = 6, expr) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  eval.parent(substitute(expr))
  invisible(path)
}

axis_colours <- c(
  navy = "#1F4E79",
  blue = "#4E79A7",
  teal = "#59A89C",
  gold = "#E0A02B",
  red = "#C44E52",
  grey = "#B7B7B7",
  dark_grey = "#555555"
)

roc_auc <- function(y, score) {
  keep <- is.finite(score) & !is.na(y)
  y <- as.logical(y[keep])
  score <- score[keep]
  n_pos <- sum(y)
  n_neg <- sum(!y)
  if (!n_pos || !n_neg) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

average_precision <- function(y, score) {
  keep <- is.finite(score) & !is.na(y)
  y <- as.logical(y[keep])
  score <- score[keep]
  if (!sum(y)) return(NA_real_)
  ord <- order(score, decreasing = TRUE)
  y <- y[ord]
  precision <- cumsum(y) / seq_along(y)
  sum(precision[y]) / sum(y)
}

roc_points <- function(y, score) {
  keep <- is.finite(score) & !is.na(y)
  y <- as.logical(y[keep])
  score <- score[keep]
  ord <- order(score, decreasing = TRUE)
  y <- y[ord]
  data.frame(
    false_positive_rate = c(0, cumsum(!y) / sum(!y)),
    true_positive_rate = c(0, cumsum(y) / sum(y))
  )
}

random_effects_dl <- function(estimate, standard_error) {
  keep <- is.finite(estimate) & is.finite(standard_error) & standard_error > 0
  yi <- estimate[keep]
  vi <- standard_error[keep]^2
  k <- length(yi)
  if (k < 2L) stop("At least two studies are required for random-effects pooling.")
  w <- 1 / vi
  fixed <- sum(w * yi) / sum(w)
  q <- sum(w * (yi - fixed)^2)
  c_value <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max(0, (q - (k - 1)) / c_value)
  wr <- 1 / (vi + tau2)
  pooled <- sum(wr * yi) / sum(wr)
  se <- sqrt(1 / sum(wr))
  list(
    k = k,
    estimate = pooled,
    standard_error = se,
    lower95 = pooled - stats::qnorm(0.975) * se,
    upper95 = pooled + stats::qnorm(0.975) * se,
    tau2 = tau2,
    q = q,
    i2 = max(0, (q - (k - 1)) / q)
  )
}

message_done <- function(name, out_dir) {
  message(name, " complete. Results: ", out_dir)
}
