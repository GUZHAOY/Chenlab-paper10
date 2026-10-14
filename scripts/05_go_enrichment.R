source("scripts/common.R")
require_packages(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2"))
dir.create(file.path(out_dir, "GO"), showWarnings = FALSE)
for (i in seq_len(nrow(comparisons))) {
  cfg <- comparisons[i, ]
  res <- read_degs(cfg$comparison)
  # SYMBOL-based featureCounts annotation. Supply ENTREZ IDs only once per gene.
  mapped <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = res$gene,
    column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  for (direction in c("Up", "Down")) {
    ids <- unique(stats::na.omit(unname(mapped[res$gene[res$status == direction]])))
    prefix <- file.path(out_dir, "GO", paste(cfg$comparison, direction, sep = "_"))
    write.csv(data.frame(ENTREZID = ids), paste0(prefix, "_input.csv"), row.names = FALSE)
    if (!length(ids)) next
    for (ont in c("BP", "MF", "CC")) {
      # Historical workflow used the OrgDb annotation background (no universe).
      go <- clusterProfiler::enrichGO(ids, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "ENTREZID", ont = ont, pAdjustMethod = "BH", pvalueCutoff = 0.05,
        qvalueCutoff = 0.2, readable = TRUE)
      tab <- as.data.frame(go)
      write.csv(tab, paste0(prefix, "_", ont, ".csv"), row.names = FALSE)
      if (nrow(tab)) {
        p <- enrichplot::dotplot(go, showCategory = 10, title = paste(cfg$comparison, direction, ont))
        ggplot2::ggsave(paste0(prefix, "_", ont, ".pdf"), p, width = 10, height = 6)
      }
    }
  }
}
save_session("05_GO")
