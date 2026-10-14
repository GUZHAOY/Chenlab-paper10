# Release validation

Validation date: 2026-10-14. Local raw gene counts were used only for validation;
neither those counts nor generated figures/results are distributed in this
code-only repository.

All eight R files parse successfully. The six analysis entry points were run on
real project inputs. The input audit covered 30 unique samples, each with 60,606
unique gene symbols and no duplicate IDs. Gene sets and lengths were checked at
runtime before combining samples.

## Differential expression

The default `padj < 0.05`, `abs(log2FC) > 2` thresholds produced:

| Comparison | Up | Down |
| --- | ---: | ---: |
| RAC1-KD versus control | 250 | 96 |
| DIAPH1-KD versus control | 166 | 48 |
| ACTA receiving, initial | 6 | 1 |
| ACTA receiving, later sorted | 88 | 745 |

The first two comparisons match the DEG counts in the supplied V51 manuscript.
This does not establish exact reproduction of all plotted values, enrichment
terms or every panel's sample selection.

## Plots and enrichment

PCA, all five configured heatmap profiles (absolute expression and row Z-score),
all four DEG scatter plots, GO BP/MF/CC workflows and both receiving-cell KEGG
profiles executed. Empty enrichment results are valid outputs and do not imply
an execution error. Current KEGG gave two significant pathways for the initial
receiving profile and none for the later sorted profile under these cutoffs;
the manuscript's TGF-beta pathway enrichment was not reproduced under this
uniform significant-DEG workflow. GO results also differ from the manuscript's
historical terms. The old enrichment scripts often selected genes from the full
DEG table by fold change alone, without reapplying `padj`; this release uses the
same significant DEG lists consistently. Database changes additionally affect
term membership and significance. Consult the historical source-data analysis
before using this release to replace the published enrichment panels.

Representative PCA, lineage heatmap, ACTA-OE heatmap, RAC1 scatter, RAC1 GO-BP,
DIAPH1 GO-BP and receiving-cell KEGG plots were rendered and visually inspected
for labels, legends and layout. Generated PDFs were checked for readability and
page presence. This is representative visual QA, not inspection of every PDF.

## Validation environment

- R 4.5.0
- DESeq2 1.48.2
- ggplot2 4.0.2
- ComplexHeatmap 2.24.1
- circlize 0.4.17
- clusterProfiler 4.16.0
- org.Hs.eg.db 3.21.0
- enrichplot 1.28.4

These differ from the historical manuscript versions. `sessionInfo()` is written
by each analysis entry point, so future executions retain their actual software
environment. The installation helper was syntax-checked; no package installation
was necessary in the validation environment.
