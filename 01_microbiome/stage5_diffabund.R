suppressMessages({library(ANCOMBC); library(TreeSummarizedExperiment); library(Maaslin2)})
setwd("/rds/projects/e/elhamsak-dissertation-autism/Microbiome/pooled")

## ---- load RAW counts (integers), drop Unassigned, align metadata ----
otu <- read.delim("merged_export/feature-table.tsv", skip=1, row.names=1, check.names=FALSE)
otu <- otu[!grepl("^Unassigned", rownames(otu)), , drop=FALSE]
otu <- round(as.matrix(otu))                       # ensure integer counts
meta <- read.delim("metadata_pooled.tsv", row.names=1, check.names=FALSE)
meta <- meta[rownames(meta) != "#q2:types", , drop=FALSE]
meta <- meta[colnames(otu), , drop=FALSE]
meta$diagnosis <- factor(meta$diagnosis, levels=c("Control","ASD"))  # Control = reference
meta$study     <- factor(meta$study)
stopifnot(all(colnames(otu) == rownames(meta)))
cat("samples:", ncol(otu), " genera:", nrow(otu), "\n\n")

## ================= ANCOM-BC2 =================
cat("===== running ANCOM-BC2 (study + diagnosis) =====\n")
tse <- TreeSummarizedExperiment(assays = list(counts = otu), colData = meta)
set.seed(1)
abc <- ancombc2(data = tse, assay_name = "counts",
                fix_formula = "study + diagnosis",
                group = "diagnosis", struc_zero = TRUE,
                prv_cut = 0.10, p_adj_method = "BH", verbose = FALSE)
ar <- abc$res
lfc_col <- grep("^lfc_diagnosis", names(ar), value=TRUE)
q_col   <- grep("^q_diagnosis",   names(ar), value=TRUE)
diff_col<- grep("^diff_diagnosis",names(ar), value=TRUE)
ancom <- data.frame(taxon=ar$taxon, lfc=ar[[lfc_col]], q=ar[[q_col]], sig=ar[[diff_col]])
ancom_hits <- ancom[ancom$sig==TRUE & !is.na(ancom$sig), ]
cat("ANCOM-BC2 significant genera:", nrow(ancom_hits), "\n")
write.csv(ancom, "stage5_ancombc2_all.csv", row.names=FALSE)

## ================= MaAsLin2 =================
cat("\n===== running MaAsLin2 (study + diagnosis) =====\n")
fit <- Maaslin2(input_data = as.data.frame(t(otu)),
                input_metadata = meta,
                output = "maaslin2_out",
                fixed_effects = c("study","diagnosis"),
                reference = c("study,PRJNA453621;diagnosis,Control"),
                normalization = "TSS", transform = "LOG",
                min_prevalence = 0.10, plot_heatmap = FALSE, plot_scatter = FALSE)
mr <- fit$results
mr_dx <- mr[mr$metadata=="diagnosis", ]
maas_hits <- mr_dx[mr_dx$qval < 0.05, ]
cat("MaAsLin2 significant genera (diagnosis):", nrow(maas_hits), "\n")

## ================= CONCORDANCE (both methods, same direction) =================
cat("\n===== CONCORDANCE: significant in BOTH, same direction =====\n")
a <- ancom_hits[, c("taxon","lfc")];      names(a) <- c("taxon","ancom_lfc")
m <- maas_hits[, c("feature","coef","qval")]; names(m) <- c("taxon","maas_coef","maas_q")
# MaAsLin2 may sanitize taxon names (dots for special chars); match on cleaned names
a$key <- make.names(a$taxon); m$key <- make.names(m$taxon)
both <- merge(a, m, by="key")
both$concordant <- sign(both$ancom_lfc) == sign(both$maas_coef)
core <- both[both$concordant, ]
core <- core[order(-abs(core$ancom_lfc)), ]
cat("Concordant core (both sig, same direction):", nrow(core), "genera\n")
if (nrow(core)>0) {
  for (i in seq_len(nrow(core))) {
    d <- ifelse(core$ancom_lfc[i]>0, "ASD-up", "ASD-down")
    cat(sprintf("  %-45s %s  (ANCOM lfc=%+.2f, MaAsLin coef=%+.2f)\n",
                substr(core$taxon.x[i],1,45), d, core$ancom_lfc[i], core$maas_coef[i]))
  }
}
write.csv(core, "stage5_concordant_core.csv", row.names=FALSE)
cat("\nSaved: stage5_ancombc2_all.csv, maaslin2_out/, stage5_concordant_core.csv\n")
