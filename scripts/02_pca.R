source("scripts/common.R")
require_packages(c("DESeq2", "ggplot2"))
ids <- c("CGS223", "CGS224", "CGS243", "CGS244", "CGS247", "CGS248",
         "CGS266", "CGS267", "CGS517", "CGS518", "CGS532", "CGS533")
x <- read_counts(ids)
counts <- x$counts[rowSums(x$counts > 1) > 0, , drop = FALSE]
meta <- data.frame(group = factor(x$meta$group), row.names = ids)
dds <- DESeq2::DESeqDataSetFromMatrix(counts, meta, design = ~group)
dds <- DESeq2::DESeq(dds)
vsd <- DESeq2::vst(dds, blind = FALSE)
# plotPCA uses the 500 most variable genes by default. Labels use its own variance.
pca <- DESeq2::plotPCA(vsd, intgroup = "group", ntop = 500, returnData = TRUE)
percent <- 100 * attr(pca, "percentVar")
p <- ggplot2::ggplot(pca, ggplot2::aes(PC1, PC2, color = group)) +
  ggplot2::geom_point(size = 4, alpha = 0.8) + ggplot2::theme_classic(base_size = 12) +
  ggplot2::labs(x = sprintf("PC1 (%.1f%%)", percent[1]), y = sprintf("PC2 (%.1f%%)", percent[2]), color = NULL)
ggplot2::ggsave(file.path(out_dir, "PCA_Fig2E.pdf"), p, width = 7, height = 5)
write.csv(pca, file.path(out_dir, "PCA_scores.csv"), row.names = FALSE)
save_session("02_PCA")
