#!/usr/bin/env Rscript
suppressMessages(library(Maaslin2))
setwd("/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA880407/maaslin2")

raw <- read.delim("genus_raw.tsv", skip=1, header=TRUE, row.names=1,
                  check.names=FALSE, comment.char="")
meta <- read.delim("../metadata/metadata_qiime2.tsv", header=TRUE, check.names=FALSE)
rownames(meta) <- meta$`sample-id`
raw <- raw[, rownames(meta)]

feat <- as.data.frame(t(raw))
meta$diagnosis <- factor(meta$diagnosis, levels=c("HC","SZ"))  # HC = reference (matches ANCOM-BC)

fit <- Maaslin2(
  input_data     = feat,
  input_metadata = meta,
  output         = "maaslin2_output",
  fixed_effects  = c("diagnosis"),
  reference      = c("diagnosis,HC"),
  normalization  = "TSS",
  transform      = "LOG",
  min_prevalence = 0.10,
  max_significance = 0.05,
  standardize    = FALSE,
  plot_heatmap   = FALSE, plot_scatter = FALSE,
  cores = 4)

res <- fit$results
dis <- res[res$metadata=="diagnosis" & res$value=="SZ", ]
dis <- dis[order(dis$qval), ]
sig <- dis[dis$qval < 0.05, ]
write.table(dis, "maaslin2_SZ_vs_HC.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat("\n===== MaAsLin2: diagnosis (SZ vs HC) =====\n")
cat("total tested:", nrow(dis), "  significant (q<0.05):", nrow(sig), "\n")
if (nrow(sig)>0) print(sig[, c("feature","coef","qval")], row.names=FALSE)
