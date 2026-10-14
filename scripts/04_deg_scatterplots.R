source("scripts/common.R")
require_packages("ggplot2")
dir.create(file.path(out_dir, "scatterplots"), showWarnings = FALSE)
for (i in seq_len(nrow(comparisons))) {
  cfg <- comparisons[i, ]
  x <- read_counts(comparison_samples(cfg))
  res <- read_degs(cfg$comparison)
  # Historical DEG panels plot mean ln(TPM+1), rather than a volcano plot.
  expr <- log(as_tpm(x) + 1)
  n_control <- length(strsplit(cfg$control, ";", fixed = TRUE)[[1]])
  plot_data <- data.frame(gene = rownames(expr),
    control = rowMeans(expr[, seq_len(n_control), drop = FALSE]),
    treatment = rowMeans(expr[, seq.int(n_control + 1, ncol(expr)), drop = FALSE]))
  plot_data$status <- res$status[match(plot_data$gene, res$gene)]
  if (anyNA(plot_data$status)) stop("DEG/input genes differ")
  plot_data$status <- factor(plot_data$status, levels = c("Not significant", "Down", "Up"))
  plot_data <- plot_data[order(plot_data$status), ]
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(control, treatment, color = status)) +
    ggplot2::geom_point(size = 0.7, alpha = 0.6) + ggplot2::theme_classic(base_size = 12) +
    ggplot2::scale_color_manual(values = c("Not significant" = "#D9D9B2", "Down" = "#5A7CAF", "Up" = "#C96455"), drop = FALSE) +
    ggplot2::labs(title = cfg$comparison, x = "Control mean ln(TPM+1)",
                 y = "Treatment mean ln(TPM+1)", color = NULL,
                 subtitle = sprintf("padj < %g; |log2FC| > %g", cfg$padj_cutoff, cfg$lfc_cutoff))
  ggplot2::ggsave(file.path(out_dir, "scatterplots", paste0(cfg$comparison, ".pdf")), p, width = 7, height = 6)
  write.csv(plot_data, file.path(out_dir, "scatterplots", paste0(cfg$comparison, "_plotted_data.csv")), row.names = FALSE)
}
save_session("04_scatterplots")
