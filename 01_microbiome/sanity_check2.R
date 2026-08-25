setwd("/rds/projects/e/elhamsak-dissertation-autism/Microbiome/pooled")
meta <- read.delim("metadata_pooled.tsv", check.names=FALSE)
meta <- meta[meta[[1]]!="#q2:types",]
ids <- meta[[1]]
cat("=========== SANITY CHECK 2: cross-layer sample consistency ===========\n\n")

# pathway table (Stage 6)
p <- read.delim("pooled_pathways.tsv", row.names=1, check.names=FALSE)
cat("PATHWAY TABLE (Stage 6):\n")
cat("  samples:", ncol(p), "\n")
cat("  identical sample set to metadata?", setequal(colnames(p), ids), "(expect TRUE)\n")
cat("  any pathway-table samples missing from metadata:", sum(!colnames(p) %in% ids), "\n\n")

# UniFrac distance matrix (rarefied, so may be < 661 if any dropped)
u <- read.delim("unifrac/unweighted_export/distance-matrix.tsv", row.names=1, check.names=FALSE)
cat("UNIFRAC MATRIX (rarefied 4030):\n")
cat("  samples:", ncol(u), "\n")
cat("  all UniFrac samples in metadata?", all(colnames(u) %in% ids), "(expect TRUE)\n")
cat("  (n<661 is OK here only if rarefaction dropped samples - you reported 0 dropped)\n\n")

# pathway meta diagnosis split matches
pm <- read.delim("pooled_pathways_meta.tsv", check.names=FALSE)
pm <- pm[pm[[1]]!="#q2:types",]
cat("PATHWAY METADATA diagnosis split:\n")
print(table(pm$diagnosis))
cat("  (expect ASD 396, Control 265)\n")
