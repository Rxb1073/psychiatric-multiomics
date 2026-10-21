suppressMessages({ library(Maaslin2) })
MB   <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome"
PICR <- file.path(MB, "picrust2")
HARM <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/harmonization"
OUT  <- file.path(PICR, "pathway_DA")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
read_path <- function(d) {
  f <- file.path(PICR, d, "picrust2_output/pathways_out/path_abun_unstrat.tsv.gz")
  read.delim(gzfile(f), row.names=1, check.names=FALSE)
}
p1 <- read_path("PRJNA591924"); p2 <- read_path("PRJNA776170")
cat("PRJNA591924:", nrow(p1), "pathways x", ncol(p1), "samples\n")
cat("PRJNA776170:", nrow(p2), "pathways x", ncol(p2), "samples (pre-filter)\n")
read_harm <- function(name) {
  h <- read.delim(file.path(HARM, paste0(name, "_harmonized.tsv")), check.names=FALSE)
  h <- h[!grepl("^#", h[[1]]), ]
  setNames(h$disease, h$`sample-id`)
}
lab1 <- read_harm("PRJNA591924"); lab2 <- read_harm("PRJNA776170")
k1 <- intersect(colnames(p1), names(lab1)); p1 <- p1[, k1]
k2 <- intersect(colnames(p2), names(lab2)); p2 <- p2[, k2]
cat("after filter: PRJNA591924", ncol(p1), " PRJNA776170", ncol(p2), "\n")
common <- intersect(rownames(p1), rownames(p2))
cat("common MetaCyc pathways:", length(common), "\n")
merged <- cbind(p1[common, ], p2[common, ])
meta <- data.frame(sample=colnames(merged),
  disease=c(lab1[k1], lab2[k2]),
  study=c(rep("PRJNA591924", ncol(p1)), rep("PRJNA776170", ncol(p2))),
  row.names=colnames(merged), stringsAsFactors=FALSE)
meta$disease <- factor(meta$disease, levels=c("Control","MDD"))
cat("\npooled:", nrow(merged), "pathways x", ncol(merged), "samples\n")
print(table(meta$study, meta$disease))
feat <- as.data.frame(t(merged))
fit <- Maaslin2(input_data=feat, input_metadata=meta,
  output=file.path(OUT, "maaslin2_pathways"),
  fixed_effects=c("study","disease"), reference=c("disease,Control"),
  normalization="TSS", transform="LOG", min_prevalence=0.10,
  standardize=FALSE, plot_heatmap=FALSE, plot_scatter=FALSE)
res <- fit$results; res <- res[res$metadata=="disease", ]; res <- res[order(res$qval), ]
write.table(res, file.path(OUT, "pathway_MDD_vs_Control.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
nsig <- sum(res$qval < 0.05, na.rm=TRUE)
cat(sprintf("\n=== significant pathways (q<0.05): %d ===\n", nsig))
if (nsig>0) {
  sig <- res[res$qval<0.05, c("feature","coef","pval","qval")]
  sig$direction <- ifelse(sig$coef<0, "down_in_MDD", "up_in_MDD")
  print(head(sig, 30), row.names=FALSE)
  cat(sprintf("\n  down in MDD: %d  |  up in MDD: %d\n", sum(sig$coef<0), sum(sig$coef>0)))
}
cat("\n===== BUTYRATE / SCFA PATHWAYS =====\n")
but <- res[grepl("butanoate|butyrate|PWY-5022|CENTFERM|PWY-5676|PWY-6590|P162-PWY|propano|fermentation", res$feature, ignore.case=TRUE), ]
if (nrow(but)) { but$direction <- ifelse(but$coef<0,"down_in_MDD","up_in_MDD")
  print(but[, c("feature","coef","pval","qval","direction")], row.names=FALSE)
} else cat("  none matched by name/ID\n")
saveRDS(list(res=res, meta=meta, merged=merged), file.path(OUT,"pathway_DA_results.rds"))
cat("\nWROTE pathway_MDD_vs_Control.tsv\n")
