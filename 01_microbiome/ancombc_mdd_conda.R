#!/usr/bin/env Rscript
suppressMessages({ library(ANCOMBC); library(TreeSummarizedExperiment) })
setwd("/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/merged_genus")

raw <- read.delim("merged_genus_raw.tsv", skip=1, header=TRUE, row.names=1,
                  check.names=FALSE, comment.char="")
meta <- read.delim("merged_genus_metadata.tsv", header=TRUE, check.names=FALSE)
rownames(meta) <- meta$`sample-id`
raw <- raw[, rownames(meta)]
meta$disease <- factor(meta$disease, levels=c("Control","MDD"))
meta$study   <- factor(meta$study)

tse <- TreeSummarizedExperiment(assays=list(counts=as.matrix(raw)), colData=meta)
cat("TSE:", nrow(tse), "genera x", ncol(tse), "samples\n")

set.seed(1)
res <- ancombc2(
  data=tse, assay_name="counts", tax_level=NULL,
  fix_formula="study + disease", group="disease",
  struc_zero=TRUE, neg_lb=TRUE, p_adj_method="BH",
  prv_cut=0.10, lib_cut=0, alpha=0.05,
  global=FALSE, pairwise=FALSE, n_cl=4, verbose=TRUE)

out <- res$res
lfc_col  <- grep("^lfc_disease",  names(out), value=TRUE)[1]
q_col    <- grep("^q_disease",    names(out), value=TRUE)[1]
diff_col <- grep("^diff_disease", names(out), value=TRUE)[1]
cat("using columns:", lfc_col, q_col, diff_col, "\n")

tab <- data.frame(genus=out$taxon,
  log2FC_MDD_vs_Control=out[[lfc_col]],
  q_value=out[[q_col]], significant=out[[diff_col]], check.names=FALSE)
tab <- tab[order(tab$q_value), ]
write.table(tab, "ancombc2_MDD_vs_Control.tsv", sep="\t", quote=FALSE, row.names=FALSE)

sig <- tab[which(tab$significant), ]
cat("\n===== ANCOM-BC2: significant genera (q<0.05, study-adjusted) =====\n")
cat("total significant:", nrow(sig), "\n")
if (nrow(sig)>0) print(sig[, c("genus","log2FC_MDD_vs_Control","q_value")], row.names=FALSE)
saveRDS(res, "ancombc2_result.rds")
cat("\nWROTE: ancombc2_MDD_vs_Control.tsv, ancombc2_result.rds\n")
