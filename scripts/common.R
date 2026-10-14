# Run scripts from the repository root. No data are bundled with this repository.
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1) args[1] else "data"
out_dir <- if (length(args) >= 2) args[2] else "results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
samples <- read.delim("config/samples.tsv", check.names = FALSE)
comparisons <- read.delim("config/comparisons.tsv", check.names = FALSE)
require_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
}
save_session <- function(label) {
  writeLines(capture.output(sessionInfo()), file.path(out_dir, paste0(label, "_sessionInfo.txt")))
}
read_counts <- function(ids) {
  meta <- samples[match(ids, samples$sample_id), , drop = FALSE]
  if (anyNA(meta$sample_id) || anyDuplicated(ids)) stop("Unknown or duplicate sample IDs")
  tables <- lapply(meta$file, function(f) {
    path <- file.path(data_dir, f)
    if (!file.exists(path)) stop("Missing input: ", path, ". See README for GEO/input instructions.")
    x <- read.delim(path, comment.char = "#", check.names = FALSE)
    if (!all(c("Geneid", "Length") %in% names(x))) stop("Expected featureCounts Geneid and Length: ", path)
    cols <- setdiff(names(x), c("Geneid", "Chr", "Start", "End", "Strand", "Length"))
    if (length(cols) != 1) stop("Expected one count column per sample: ", path)
    if (anyDuplicated(x$Geneid)) stop("Duplicate gene IDs in ", path)
    v <- as.numeric(x[[cols]])
    if (any(!is.finite(v) | v < 0 | v != floor(v))) stop("Non-integer/raw-count input: ", path)
    list(gene = x$Geneid, count = v, length = as.numeric(x$Length))
  })
  genes <- tables[[1]]$gene
  for (x in tables) {
    if (!setequal(genes, x$gene)) stop("Gene sets differ between samples; align/annotate inputs first")
    if (!isTRUE(all.equal(tables[[1]]$length, x$length[match(genes, x$gene)])))
      stop("Gene lengths differ; use a consistent annotation")
  }
  mat <- do.call(cbind, lapply(tables, function(x) x$count[match(genes, x$gene)]))
  dimnames(mat) <- list(genes, ids)
  list(counts = mat, lengths = tables[[1]]$length, meta = meta)
}
as_tpm <- function(x) {
  if (any(!is.finite(x$lengths) | x$lengths <= 0)) stop("Invalid gene lengths")
  rpk <- sweep(x$counts, 1, x$lengths / 1000, "/")
  if (any(colSums(rpk) <= 0)) stop("Empty library")
  sweep(rpk, 2, colSums(rpk) / 1e6, "/")
}
comparison_samples <- function(row) c(strsplit(row$control, ";", fixed = TRUE)[[1]],
                                      strsplit(row$treatment, ";", fixed = TRUE)[[1]])
classify_degs <- function(res, padj_cutoff, lfc_cutoff) {
  status <- rep("Not significant", nrow(res))
  sig <- !is.na(res$padj) & !is.na(res$log2FoldChange) & res$padj < padj_cutoff
  status[sig & res$log2FoldChange > lfc_cutoff] <- "Up"
  status[sig & res$log2FoldChange < -lfc_cutoff] <- "Down"
  factor(status, levels = c("Not significant", "Down", "Up"))
}
read_degs <- function(name) {
  p <- file.path(out_dir, "DEG", paste0(name, "_all_genes.csv"))
  if (!file.exists(p)) stop("Run scripts/01_differential_expression.R first: ", p)
  read.csv(p, check.names = FALSE)
}
