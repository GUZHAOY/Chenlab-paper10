source("scripts/common.R")
require_packages(c("ComplexHeatmap", "circlize"))
genes <- read.delim("config/genes.tsv")
panels <- read.delim("config/heatmaps.tsv")
dir.create(file.path(out_dir, "heatmaps"), showWarnings = FALSE)
for (i in seq_len(nrow(panels))) {
  cfg <- panels[i, ]
  ids <- strsplit(cfg$samples, ";", fixed = TRUE)[[1]]
  x <- read_counts(ids)
  wanted <- genes[genes$panel == cfg$gene_panel, , drop = FALSE]
  missing <- wanted$gene[!wanted$gene %in% rownames(x$counts)]
  write.csv(data.frame(gene = missing), file.path(out_dir, "heatmaps", paste0(cfg$panel, "_missing_genes.csv")), row.names = FALSE)
  wanted <- wanted[wanted$gene %in% rownames(x$counts), , drop = FALSE]
  if (!nrow(wanted)) stop("No matching genes for ", cfg$panel)
  # Historical Fig2A workflow normalized counts per million without gene lengths.
  # Preserve that calculation and label it as CPM; other panels use real TPM.
  expr <- if (cfg$normalization == "CPM") {
    if (any(colSums(x$counts) == 0)) stop("Empty library")
    sweep(x$counts, 2, colSums(x$counts) / 1e6, "/")
  } else as_tpm(x)
  mat <- log2(expr[wanted$gene, , drop = FALSE] + 1)
  write.csv(mat, file.path(out_dir, "heatmaps", paste0(cfg$panel, "_log2_expression.csv")))
  splits <- factor(wanted$category, levels = unique(wanted$category))
  for (mode in c("log2", "zscore")) {
    z <- mat
    legend <- paste0("log2(", cfg$normalization, "+1)")
    if (mode == "zscore") {
      z <- t(scale(t(mat)))
      z[!is.finite(z)] <- 0 # A constant-expression row has zero deviation.
      legend <- "Row Z-score"
      breaks <- c(-2, 0, 2)
    } else {
      limits <- range(z)
      if (diff(limits) == 0) limits <- limits + c(-0.5, 0.5)
      breaks <- c(limits[1], mean(limits), limits[2])
    }
    write.csv(z, file.path(out_dir, "heatmaps", paste0(cfg$panel, "_", mode, "_plotted_matrix.csv")))
    pdf(file.path(out_dir, "heatmaps", paste0(cfg$panel, "_", mode, ".pdf")),
        width = max(7, length(ids) * 0.42 + 3), height = max(5, nrow(z) * 0.16 + 2))
    ht <- ComplexHeatmap::Heatmap(z, name = legend, row_split = splits,
      cluster_rows = FALSE, cluster_columns = FALSE, cluster_row_slices = FALSE,
      column_labels = ids, column_names_rot = 55,
      row_names_gp = grid::gpar(fontsize = 8), column_names_gp = grid::gpar(fontsize = 8),
      col = circlize::colorRamp2(breaks, c("#38629E", "#F3F4B7", "#BF4431")))
    ComplexHeatmap::draw(ht)
    dev.off()
  }
}
save_session("03_heatmaps")
