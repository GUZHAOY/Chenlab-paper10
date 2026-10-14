source("scripts/common.R")
require_packages("DESeq2")
dir.create(file.path(out_dir, "DEG"), showWarnings = FALSE)
summary_rows <- list()
for (i in seq_len(nrow(comparisons))) {
  cfg <- comparisons[i, ]
  x <- read_counts(comparison_samples(cfg))
  n_control <- length(strsplit(cfg$control, ";", fixed = TRUE)[[1]])
  condition <- factor(c(rep("control", n_control), rep("treatment", ncol(x$counts) - n_control)),
                      levels = c("control", "treatment"))
  # Raw integer counts only; TPM/log-transformed values must not enter DESeq2.
  dds <- DESeq2::DESeqDataSetFromMatrix(x$counts,
    data.frame(condition, row.names = colnames(x$counts)), design = ~condition)
  dds <- DESeq2::DESeq(dds)
  res <- as.data.frame(DESeq2::results(dds, contrast = c("condition", "treatment", "control"),
                                     alpha = cfg$padj_cutoff))
  res$gene <- rownames(res)
  res$status <- classify_degs(res, cfg$padj_cutoff, cfg$lfc_cutoff)
  res <- res[, c("gene", setdiff(names(res), "gene"))]
  prefix <- file.path(out_dir, "DEG", cfg$comparison)
  write.csv(res, paste0(prefix, "_all_genes.csv"), row.names = FALSE)
  for (direction in c("Up", "Down"))
    write.csv(res[res$status == direction, ], paste0(prefix, "_", tolower(direction), ".csv"), row.names = FALSE)
  summary_rows[[i]] <- data.frame(comparison = cfg$comparison, up = sum(res$status == "Up"),
    down = sum(res$status == "Down"), padj_cutoff = cfg$padj_cutoff, abs_log2fc_cutoff = cfg$lfc_cutoff)
}
write.csv(do.call(rbind, summary_rows), file.path(out_dir, "DEG", "summary.csv"), row.names = FALSE)
save_session("01_DEG")
