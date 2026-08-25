library(vegan)
setwd("/rds/projects/e/elhamsak-dissertation-autism/Microbiome/pooled")
raw  <- read.delim("merged_export/feature-table.tsv", skip=1, row.names=1, check.names=FALSE)
raw  <- raw[!grepl("^Unassigned", rownames(raw)), , drop=FALSE]
corr <- read.delim("genus_batch_corrected.tsv", row.names=1, check.names=FALSE)
meta <- read.delim("metadata_pooled.tsv", row.names=1, check.names=FALSE)
meta <- meta[rownames(meta) != "#q2:types", , drop=FALSE]
meta <- meta[colnames(corr), , drop=FALSE]
meta$diagnosis <- factor(meta$diagnosis, levels=c("Control","ASD"))
d_bray <- vegdist(t(corr), method="bray")
bd <- betadisper(d_bray, meta$diagnosis)
cat("Mean distance-to-centroid (higher = more dispersed):\n")
print(tapply(bd$distances, meta$diagnosis, mean))
