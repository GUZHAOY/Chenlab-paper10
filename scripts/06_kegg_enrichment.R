source("scripts/common.R")
require_packages(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2"))
dir.create(file.path(out_dir, "KEGG"), showWarnings = FALSE)
# Both receiving-cell profiles are retained because the local archive has two
# sorting experiments; choose the final panel's samples in comparisons.tsv.
cfgs <- comparisons[grepl("receiving", comparisons$comparison), , drop = FALSE]
for (i in seq_len(nrow(cfgs))) {
  cfg <- cfgs[i, ]
  res <- read_degs(cfg$comparison)
  mapped <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = res$gene,
    column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  ids <- unique(stats::na.omit(unname(mapped[res$gene[res$status == "Up"]])))
  prefix <- file.path(out_dir, "KEGG", cfg$comparison)
  write.csv(data.frame(ENTREZID = ids), paste0(prefix, "_up_input.csv"), row.names = FALSE)
  if (!length(ids)) next
  # Uses the current KEGG service; keep output/session date for reproducibility.
  ke <- clusterProfiler::enrichKEGG(ids, organism = "hsa", pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  tab <- as.data.frame(ke)
  write.csv(tab, paste0(prefix, ".csv"), row.names = FALSE)
  if (nrow(tab)) {
    p <- enrichplot::dotplot(ke, showCategory = 10, title = paste(cfg$comparison, "upregulated KEGG"))
    ggplot2::ggsave(paste0(prefix, ".pdf"), p, width = 10, height = 6)
  }
}
save_session("06_KEGG")
