#!/usr/bin/env Rscript
suppressMessages({
  personal <- Sys.getenv("R_LIBS_USER"); if (nchar(personal)==0) personal <- path.expand("~/R/mmuphin_lib")
  .libPaths(c(personal, .libPaths()))
  library(MMUPHin); library(vegan)
})
setwd("/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/merged_genus")
HARM <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/harmonization"

raw <- read.delim("merged_genus_raw.tsv", skip=1, header=TRUE, row.names=1,
                  check.names=FALSE, comment.char="")
cat("raw genus table:", nrow(raw), "genera x", ncol(raw), "samples\n")

meta_all <- read.delim(file.path(HARM,"master_sample_sheet.tsv"), header=TRUE, check.names=FALSE)
meta <- meta_all[meta_all$layer=="microbiome", ]
rownames(meta) <- meta$`sample-id`
common <- intersect(colnames(raw), rownames(meta))
cat("table:", ncol(raw), " metadata:", nrow(meta), " common:", length(common), "\n")
stopifnot(length(common) == 158)
raw  <- raw[, common]
meta <- meta[common, ]
meta$study   <- factor(meta$study)
meta$disease <- factor(meta$disease, levels=c("Control","MDD"))
cat("study x disease:\n"); print(table(meta$study, meta$disease))

prev <- rowSums(raw > 0) / ncol(raw)
keep <- prev >= 0.10
cat("genera before/after prevalence >=10%:", nrow(raw), "->", sum(keep), "\n")
feat <- as.matrix(raw[keep, ])

fit <- adjust_batch(feature_abd=feat, batch="study", covariates="disease",
                    data=meta, control=list(verbose=FALSE))
corrected <- fit$feature_abd_adj
cat("corrected:", nrow(corrected), "genera x", ncol(corrected), "samples\n")

bc_before <- vegdist(t(feat),      method="bray")
bc_after  <- vegdist(t(corrected), method="bray")
a_before <- adonis2(bc_before ~ study + disease, data=meta, permutations=999, by="margin")
a_after  <- adonis2(bc_after  ~ study + disease, data=meta, permutations=999, by="margin")
cat("\n--- PERMANOVA BEFORE ---\n"); print(a_before)
cat("\n--- PERMANOVA AFTER ---\n");  print(a_after)

write.table(data.frame(genus=rownames(corrected), corrected, check.names=FALSE),
            "merged_genus_corrected.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(data.frame(`sample-id`=rownames(meta),
            meta[,c("disease","study","tissue_region")], check.names=FALSE),
            "merged_genus_metadata.tsv", sep="\t", quote=FALSE, row.names=FALSE)
saveRDS(list(fit=fit, meta=meta, before=a_before, after=a_after), "mmuphin_result.rds")
cat("\nWROTE: merged_genus_corrected.tsv, merged_genus_metadata.tsv, mmuphin_result.rds\n")
