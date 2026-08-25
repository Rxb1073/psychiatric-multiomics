library(MMUPHin); library(vegan)
setwd("/rds/projects/e/elhamsak-dissertation-autism/Microbiome/pooled")

otu <- read.delim("merged_export/feature-table.tsv", skip=1, row.names=1, check.names=FALSE)
meta <- read.delim("metadata_pooled.tsv", row.names=1, check.names=FALSE)
meta <- meta[rownames(meta) != "#q2:types", , drop=FALSE]
meta <- meta[colnames(otu), , drop=FALSE]
stopifnot(all(colnames(otu) == rownames(meta)))
cat("samples:", ncol(otu), " genera:", nrow(otu), "\n")

# drop Unassigned bin if present (not a real genus)
otu <- otu[!grepl("^Unassigned", rownames(otu)), , drop=FALSE]
cat("genera after removing Unassigned:", nrow(otu), "\n")

# 10% prevalence filter (same threshold as MDD)
keep <- rowSums(otu > 0) >= 0.10 * ncol(otu)
otu_f <- otu[keep, ]
cat("genera after 10% prevalence filter:", nrow(otu_f), "of", nrow(otu), "\n")

# PERMANOVA BEFORE correction
set.seed(1)
d_before <- vegdist(t(otu_f), method="bray")
a_study_before   <- adonis2(d_before ~ study,     data=meta, permutations=999)
a_disease_before <- adonis2(d_before ~ diagnosis, data=meta, permutations=999)
cat("\nBEFORE correction:\n")
cat("  study   R2 =", round(a_study_before$R2[1],4),   " p =", a_study_before$`Pr(>F)`[1], "\n")
cat("  disease R2 =", round(a_disease_before$R2[1],4), " p =", a_disease_before$`Pr(>F)`[1], "\n")

# MMUPHin adjust_batch
fit <- adjust_batch(feature_abd = otu_f, batch = "study", covariates = "diagnosis", data = meta)
otu_corr <- fit$feature_abd_adj

# PERMANOVA AFTER correction
set.seed(1)
d_after <- vegdist(t(otu_corr), method="bray")
a_study_after   <- adonis2(d_after ~ study,     data=meta, permutations=999)
a_disease_after <- adonis2(d_after ~ diagnosis, data=meta, permutations=999)
cat("\nAFTER correction:\n")
cat("  study   R2 =", round(a_study_after$R2[1],4),   " p =", a_study_after$`Pr(>F)`[1], "\n")
cat("  disease R2 =", round(a_disease_after$R2[1],4), " p =", a_disease_after$`Pr(>F)`[1], "\n")

write.table(otu_corr, "genus_batch_corrected.tsv", sep="\t", quote=FALSE, col.names=NA)
cat("\nSaved corrected table -> pooled/genus_batch_corrected.tsv\n")
