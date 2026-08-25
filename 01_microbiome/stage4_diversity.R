library(vegan)
setwd("/rds/projects/e/elhamsak-dissertation-autism/Microbiome/pooled")

## ---- load ----
raw  <- read.delim("merged_export/feature-table.tsv", skip=1, row.names=1, check.names=FALSE)
raw  <- raw[!grepl("^Unassigned", rownames(raw)), , drop=FALSE]       # drop Unassigned
corr <- read.delim("genus_batch_corrected.tsv", row.names=1, check.names=FALSE)  # MMUPHin output
meta <- read.delim("metadata_pooled.tsv", row.names=1, check.names=FALSE)
meta <- meta[rownames(meta) != "#q2:types", , drop=FALSE]
meta <- meta[colnames(raw), , drop=FALSE]
meta$diagnosis <- factor(meta$diagnosis, levels=c("Control","ASD"))
meta$study     <- factor(meta$study)
stopifnot(all(colnames(raw) == rownames(meta)))
cat("samples:", ncol(raw), " genera(raw):", nrow(raw), " genera(corrected):", nrow(corr), "\n\n")

## ================= ALPHA (raw counts) =================
alpha <- data.frame(
  observed = specnumber(t(raw)),
  shannon  = diversity(t(raw), "shannon"),
  simpson  = diversity(t(raw), "simpson")
)
alpha <- cbind(alpha, meta[rownames(alpha), c("diagnosis","study")])

cat("===== ALPHA DIVERSITY (ASD vs Control) =====\n")
for (m in c("observed","shannon","simpson")) {
  # group medians (gives DIRECTION) + linear model controlling for study
  med <- tapply(alpha[[m]], alpha$diagnosis, median)
  lm_p <- summary(lm(alpha[[m]] ~ diagnosis + study, data=alpha))$coefficients["diagnosisASD","Pr(>|t|)"]
  wil_p <- wilcox.test(alpha[[m]] ~ diagnosis, data=alpha)$p.value
  cat(sprintf("  %-9s  Control median=%.3f  ASD median=%.3f  |  lm(+study) p=%.4g  Wilcoxon p=%.4g\n",
              m, med["Control"], med["ASD"], lm_p, wil_p))
}
write.csv(alpha, "stage4_alpha_values.csv")

## ================= BETA =================
run_beta <- function(dmat, label) {
  set.seed(1)
  # PERMANOVA: study first (covariate), then diagnosis
  pm <- adonis2(dmat ~ study + diagnosis, data=meta, permutations=999, by="terms")
  r2 <- pm["diagnosis","R2"]; p <- pm["diagnosis","Pr(>F)"]
  # PERMDISP on diagnosis
  bd <- betadisper(dmat, meta$diagnosis)
  pd <- permutest(bd, permutations=999)$tab["Groups","Pr(>F)"]
  cat(sprintf("  %-14s  PERMANOVA diagnosis R2=%.4f p=%.4g  |  PERMDISP p=%.4g %s\n",
              label, r2, p, pd,
              ifelse(pd>0.05,"(dispersion OK -> location shift)","(WARNING: dispersion differs)")))
}

cat("\n===== BETA DIVERSITY (PERMANOVA diagnosis after study; PERMDISP) =====\n")
# Bray-Curtis on CORRECTED table (per pipeline: corrected used ONLY here)
meta_c <- meta[colnames(corr), , drop=FALSE]
set.seed(1)
d_bray <- vegdist(t(corr), method="bray")
pm <- adonis2(d_bray ~ study + diagnosis, data=meta_c, permutations=999, by="terms")
bd <- betadisper(d_bray, meta_c$diagnosis); pd <- permutest(bd, permutations=999)$tab["Groups","Pr(>F)"]
cat(sprintf("  %-14s  PERMANOVA diagnosis R2=%.4f p=%.4g  |  PERMDISP p=%.4g %s\n",
            "Bray(corrected)", pm["diagnosis","R2"], pm["diagnosis","Pr(>F)"], pd,
            ifelse(pd>0.05,"(dispersion OK)","(WARNING dispersion)")))

# Aitchison (CLR) on RAW, prevalence-filtered
keep <- rowSums(raw > 0) >= 0.10 * ncol(raw)
rawf <- raw[keep, ]
logx <- log(rawf + 1)
clr  <- sweep(logx, 2, colMeans(logx), "-")
d_ait <- dist(t(clr))            # Euclidean on CLR = Aitchison
run_beta(d_ait, "Aitchison(raw)")

cat("\nSaved alpha values -> stage4_alpha_values.csv\n")
