#!/usr/bin/env Rscript
setwd("/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/merged_genus")
anc <- read.delim("ancombc2_MDD_vs_Control.tsv", check.names=FALSE)
maa <- read.delim("maaslin2_MDD_vs_Control.tsv", check.names=FALSE)
anc_sig <- anc[anc$significant %in% c(TRUE,"TRUE"), ]
maa_sig <- maa[maa$qval < 0.05, ]
last_tok <- function(x){ p <- unlist(strsplit(x, "[;.]")); p <- p[p!="" & p!="__"]; tail(p,1) }
anc_sig$key <- sapply(anc_sig$genus,   last_tok)
maa_sig$key <- sapply(maa_sig$feature, last_tok)
both <- intersect(anc_sig$key, maa_sig$key)
both <- both[!both %in% c("g__","__","")]
cat("ANCOM-BC2 significant:", nrow(anc_sig), "\n")
cat("MaAsLin2 significant:", nrow(maa_sig), "\n")
cat("Overlap (both methods):", length(both), "\n\n")
ov <- data.frame(genus=both)
ov$ancombc_log2FC <- anc_sig$log2FC_MDD_vs_Control[match(ov$genus, anc_sig$key)]
ov$ancombc_q      <- anc_sig$q_value[match(ov$genus, anc_sig$key)]
ov$maaslin_coef   <- maa_sig$coef[match(ov$genus, maa_sig$key)]
ov$maaslin_q      <- maa_sig$qval[match(ov$genus, maa_sig$key)]
ov$same_direction <- sign(ov$ancombc_log2FC) == sign(ov$maaslin_coef)
ov <- ov[order(ov$ancombc_q), ]
write.table(ov, "robust_genera_both_methods.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat("===== ROBUST GENERA (both methods) =====\n")
print(ov, row.names=FALSE)
cat("\ndirection concordant:", sum(ov$same_direction, na.rm=TRUE), "/", nrow(ov), "\n")
