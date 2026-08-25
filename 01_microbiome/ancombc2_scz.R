suppressMessages({ library(ANCOMBC); library(TreeSummarizedExperiment) })
setwd("/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA880407")
raw  <- read.delim("genus_raw.tsv", skip=1, header=TRUE, row.names=1, check.names=FALSE, comment.char="")
meta <- read.delim("metadata/metadata_qiime2.tsv", header=TRUE, check.names=FALSE)
meta <- meta[meta[[1]] != "#q2:types", ]
rownames(meta) <- meta[[1]]
stopifnot(all(rownames(meta) %in% colnames(raw)))
raw <- raw[, rownames(meta)]
meta$diagnosis <- factor(meta$diagnosis, levels=c("HC","SZ"))
cat("group sizes:\n"); print(table(meta$diagnosis))
tse <- TreeSummarizedExperiment(assays=list(counts=as.matrix(raw)), colData=meta)
cat("TSE:", nrow(tse), "genera x", ncol(tse), "samples\n")
set.seed(1)
res <- ancombc2(data=tse, assay_name="counts", tax_level=NULL,
  fix_formula="diagnosis", group="diagnosis",
  struc_zero=TRUE, neg_lb=TRUE, p_adj_method="BH",
  prv_cut=0.10, lib_cut=0, alpha=0.05,
  global=FALSE, pairwise=FALSE, n_cl=4, verbose=TRUE)
out <- res$res
lfc <- grep("^lfc_diagnosis", names(out), value=TRUE)[1]
qv  <- grep("^q_diagnosis",   names(out), value=TRUE)[1]
df  <- grep("^diff_diagnosis",names(out), value=TRUE)[1]
cat("columns:", lfc, qv, df, "\n")
tab <- data.frame(genus=out$taxon, log2FC_SZ_vs_HC=out[[lfc]], q_value=out[[qv]], significant=out[[df]])
tab <- tab[order(tab$q_value), ]
write.table(tab, "ancombc2_SZ_vs_HC.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat("\nsignificant (q<0.05):", sum(tab$significant, na.rm=TRUE), "  [ANCOM-BC v1/Holm gave 29]\n")
print(head(tab[which(tab$significant), ], 15), row.names=FALSE)
saveRDS(res, "ancombc2_result_scz.rds")
