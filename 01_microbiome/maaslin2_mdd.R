#!/usr/bin/env Rscript
suppressMessages({
  if (!requireNamespace("Maaslin2", quietly=TRUE)) {
    if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
    BiocManager::install("Maaslin2", update=FALSE, ask=FALSE)
  }
  library(Maaslin2)
})
setwd("/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/merged_genus")

## ---- raw counts (features x samples) + metadata ----
raw <- read.delim("merged_genus_raw.tsv", skip=1, header=TRUE, row.names=1,
                  check.names=FALSE, comment.char="")
meta <- read.delim("merged_genus_metadata.tsv", header=TRUE, check.names=FALSE)
rownames(meta) <- meta$`sample-id`
raw <- raw[, rownames(meta)]

## MaAsLin2 wants samples as ROWS in the feature table (samples x features)
feat <- as.data.frame(t(raw))            # samples x genera
meta$disease <- factor(meta$disease, levels=c("Control","MDD"))  # Control = reference
meta$study   <- factor(meta$study)

## ---- run MaAsLin2 ----
## TSS normalization + LOG transform (MaAsLin2 default, compositional-aware);
## disease is the effect, study is a fixed covariate; 10% prevalence to match ANCOM-BC2.
fit <- Maaslin2(
  input_data     = feat,
  input_metadata = meta,
  output         = "maaslin2_output",
  fixed_effects  = c("study","disease"),
  reference      = c("disease,Control"),
  normalization  = "TSS",
  transform      = "LOG",
  min_prevalence = 0.10,
  max_significance = 0.05,
  standardize    = FALSE,
  plot_heatmap   = FALSE, plot_scatter = FALSE,
  cores = 4)

## ---- extract the disease results ----
res <- fit$results
dis <- res[res$metadata=="disease" & res$value=="MDD", ]
dis <- dis[order(dis$qval), ]
sig <- dis[dis$qval < 0.05, ]
write.table(dis, "maaslin2_MDD_vs_Control.tsv", sep="\t", quote=FALSE, row.names=FALSE)

cat("\n===== MaAsLin2: disease (MDD vs Control), study-adjusted =====\n")
cat("total tested:", nrow(dis), "  significant (q<0.05):", nrow(sig), "\n")
if (nrow(sig)>0) print(sig[, c("feature","coef","qval")], row.names=FALSE)
cat("\nWROTE: maaslin2_MDD_vs_Control.tsv, maaslin2_output/\n")
