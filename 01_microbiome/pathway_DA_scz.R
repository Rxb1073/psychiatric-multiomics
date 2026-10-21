suppressMessages(library(Maaslin2))
PICR <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA880407/picrust2"
OUT  <- file.path(PICR, "pathway_DA")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

p <- read.delim(gzfile(file.path(PICR,"picrust2_output/pathways_out/path_abun_unstrat.tsv.gz")),
                row.names=1, check.names=FALSE)
cat("pathways x samples:", nrow(p), "x", ncol(p), "\n")

meta <- read.delim("/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA880407/metadata/metadata_qiime2.tsv",
                   check.names=FALSE)
rownames(meta) <- meta$`sample-id`
p <- p[, rownames(meta)]
meta$diagnosis <- factor(meta$diagnosis, levels=c("HC","SZ"))

feat <- as.data.frame(t(p))
fit <- Maaslin2(input_data=feat, input_metadata=meta,
  output=file.path(OUT,"maaslin2_pathways"),
  fixed_effects=c("diagnosis"), reference=c("diagnosis,HC"),
  normalization="TSS", transform="LOG", min_prevalence=0.10,
  standardize=FALSE, plot_heatmap=FALSE, plot_scatter=FALSE)

res <- fit$results; res <- res[res$metadata=="diagnosis" & res$value=="SZ", ]
res <- res[order(res$qval), ]
write.table(res, file.path(OUT,"pathway_SZ_vs_HC.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
nsig <- sum(res$qval<0.05, na.rm=TRUE)
cat(sprintf("\n=== significant pathways (q<0.05): %d ===\n", nsig))
if (nsig>0){
  sig <- res[res$qval<0.05, c("feature","coef","pval","qval")]
  sig$direction <- ifelse(sig$coef<0,"down_in_SZ","up_in_SZ")
  print(head(sig,40), row.names=FALSE)
  cat(sprintf("\n  down in SZ: %d  |  up in SZ: %d\n", sum(sig$coef<0), sum(sig$coef>0)))
}
cat("\n===== BUTYRATE / SCFA PATHWAYS =====\n")
but <- res[grepl("butanoate|butyrate|PWY-5022|CENTFERM|PWY-5676|PWY-6590|P162-PWY|propano|fermentation", res$feature, ignore.case=TRUE), ]
if (nrow(but)){ but$direction <- ifelse(but$coef<0,"down_in_SZ","up_in_SZ")
  print(but[,c("feature","coef","pval","qval","direction")], row.names=FALSE)
} else cat("  none matched\n")
saveRDS(list(res=res, meta=meta), file.path(OUT,"pathway_DA_results.rds"))
cat("\nWROTE pathway_SZ_vs_HC.tsv\n")
