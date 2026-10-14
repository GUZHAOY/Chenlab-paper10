# Chenlab paper 10: cytoneme RNA-seq analysis

Representative R workflows for the cytoneme/hESC study: differential expression,
PCA, expression heatmaps, DEG expression scatter plots, GO enrichment and KEGG
enrichment. The scripts consolidate the project's historical analysis scripts
into one workflow per plot type. They cover downstream bulk RNA-seq analysis;
FlowJo, confocal image processing and GraphPad Prism analyses are outside this
repository.

## Data

Sequencing data are deposited under **GSE287522**:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE287522

This repository contains **code and configuration only**. Obtain the relevant
processed counts from GEO, or generate gene-level featureCounts files using
hg38/GENCODE v41 as described in the manuscript. Put one raw-count file per
sample in `data/`, with the filenames in `config/samples.tsv`. GEO filenames may
need to be renamed or the configuration edited to match the local inputs.
No automatic GEO-to-CGS mapping is assumed.

Inputs must have the featureCounts columns `Geneid`, `Chr`, `Start`, `End`,
`Strand`, `Length` and one integer count column. A leading `#` command line is
supported. Gene IDs are human gene symbols (`featureCounts -g gene_name`).
Lengths must come from the same annotation across samples. Missing inputs,
duplicate genes, non-integer counts, unequal gene sets and inconsistent lengths
raise errors rather than silently altering the analysis.

## Run

Use R with Bioconductor and run all commands from the repository root:

```sh
Rscript scripts/00_install_packages.R
Rscript scripts/01_differential_expression.R data results
Rscript scripts/02_pca.R data results
Rscript scripts/03_expression_heatmaps.R data results
Rscript scripts/04_deg_scatterplots.R data results
Rscript scripts/05_go_enrichment.R data results
Rscript scripts/06_kegg_enrichment.R data results
```

The input/output arguments are optional and default to `data` and `results`.
PDF plots, plotted matrices, DEG tables, enrichment tables and `sessionInfo()`
records are written to the output directory. The installation script installs
compatible available packages; it does not recreate a frozen historical library.
KEGG requires internet access and uses the current KEGG annotation service.

## Workflow and figure correspondence

| Analysis | Script | Manuscript panels / use |
| --- | --- | --- |
| DESeq2 | `01_differential_expression.R` | Statistical input for DEG, GO and KEGG plots |
| PCA | `02_pca.R` | Figure 2E |
| Expression heatmaps | `03_expression_heatmaps.R` | Figure 2A; Figures S2D, S4E and S4H |
| DEG expression scatter | `04_deg_scatterplots.R` | Figures S2E and S2G |
| GO dot plots | `05_go_enrichment.R` | Figures S2F, S2G and S4I |
| KEGG dot plots | `06_kegg_enrichment.R` | Figure 4F |

All comparison directions are explicitly **treatment versus control**. Defaults
follow the supplied manuscript: BH-adjusted `padj < 0.05` and
`abs(log2FoldChange) > 2` (a four-fold change, not a two-fold change).
The same significant-gene lists feed scatter plots and enrichment.
DESeq2 takes raw integer counts. PCA uses DESeq2 VST and the 500 most variable
genes; axis variance percentages are computed from that same PCA.

Scatter plots use the historical mean **ln(TPM+1)** expression axes. Heatmaps
provide absolute log2 expression and row Z-score variants. TPM uses gene lengths
and per-sample scaling. A historical lineage heatmap calculated CPM without
length correction; that calculation is retained with an accurate **CPM** label.
GO runs BP, MF and CC separately, using the historical OrgDb background and BH
adjustment. Enrichment tables retain term descriptions and gene lists.

## Configuration and historical interpretation

- `config/samples.tsv`: sample IDs, filenames and the group labels from the
  selected source workflows.
- `config/comparisons.tsv`: control/treatment IDs and explicit cutoffs.
- `config/heatmaps.tsv` and `config/genes.tsv`: sample order and gene panels.
- `docs/source_workflows.tsv`: provenance of the consolidated code.

The local archive contains multiple generations of knockdown and sorting
experiments. This release selects the CGS517/518 RAC1-KD and CGS532/533 DIAPH1-KD
workflows for the main hESC comparison. It includes both the initial receiving
cell experiment (CGS625/626 versus CGS623/624) and a later sorted receiving-cell
profile (CGS677/678 versus CGS623/624). These are separate profiles, not pooled
replicates. The final Figure 4F/S4H/S4I experiment should be selected against the
final source-data sheet before claiming exact panel reproduction.

Two label inconsistencies exist in the historical lineage workflows:
CGS243/244 and CGS247/248 have opposite endoderm/mesoderm labels between the PCA
and lineage heatmap scripts; CGS185/186 were called `Formative_hPGCLC` in one
script while the manuscript describes formative hESCs. PCA retains its source
workflow's labels; the lineage heatmap displays CGS sample IDs. Confirm these
labels against the original sample submission sheet when making a final figure.
The historical `STRIT1` candidate is retained as written; absent symbols are
reported in missing-gene tables.

These are cleaned, representative analysis workflows, not a byte-for-byte
archive of every published panel. Historical cutoffs, package versions,
normalizations and sample selections differ across old scripts. No automated
change of cutoff or sample exclusion is performed to match a published DEG
count. See `docs/validation.md` for checks run on this release.

## Software

The manuscript reports R 4.3.1, DESeq2 1.40.2, ggplot2 3.4.4,
ComplexHeatmap 2.15.4 and clusterProfiler 4.9.0.2. This release additionally
uses circlize, org.Hs.eg.db, AnnotationDbi and enrichplot. Actual validation
versions are recorded separately; outputs always record the running session.
