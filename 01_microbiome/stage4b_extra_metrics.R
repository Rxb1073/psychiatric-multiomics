library(vegan)
setwd("/rds/projects/e/elhamsak-dissertation-autism/Microbiome/pooled/unifrac")

meta <- read.delim("../metadata_pooled.tsv", row.names=1, check.names=FALSE)
meta <- meta[rownames(meta) != "#q2:types", , drop=FALSE]
meta$diagnosis <- factor(meta$diagnosis, levels=c("Control","ASD"))
meta$study     <- factor(meta$study)

## ---- ALPHA: Faith's PD + Pielou evenness (same treatment as Stage 4) ----
alpha_test <- function(file, label){
  a <- read.delim(file, row.names=1, check.names=FALSE)
  v <- setNames(a[[1]], rownames(a))
  m <- meta[names(v), , drop=FALSE]
  d <- data.frame(val=as.numeric(v), diagnosis=m$diagnosis, study=m$study)
  med  <- tapply(d$val, d$diagnosis, median)
  lm_p <- summary(lm(val ~ diagnosis + study, data=d))$coefficients["diagnosisASD","Pr(>|t|)"]
  wil_p<- wilcox.test(val ~ diagnosis, data=d)$p.value
  cat(sprintf("  %-16s n=%d  Control median=%.3f  ASD median=%.3f  |  lm(+study) p=%.4g  Wilcoxon p=%.4g\n",
              label, nrow(d), med["Control"], med["ASD"], lm_p, wil_p))
}
cat("===== ALPHA (rarefied 4030) =====\n")
alpha_test("faith_export/alpha-diversity.tsv", "Faith's PD")
alpha_test("even_export/alpha-diversity.tsv",  "Pielou evenness")

## ---- BETA: Jaccard (same treatment as Stage 4) ----
cat("\n===== BETA: Jaccard =====\n")
d <- as.dist(as.matrix(read.delim("jaccard_export/distance-matrix.tsv", row.names=1, check.names=FALSE)))
m <- meta[labels(d), , drop=FALSE]
set.seed(1)
pm <- adonis2(d ~ study + diagnosis, data=m, permutations=999, by="terms")
bd <- betadisper(d, m$diagnosis)
pd <- permutest(bd, permutations=999)$tab["Groups","Pr(>F)"]
cat(sprintf("  study     R2=%.4f p=%.4g\n", pm["study","R2"], pm["study","Pr(>F)"]))
cat(sprintf("  diagnosis R2=%.4f p=%.4g\n", pm["diagnosis","R2"], pm["diagnosis","Pr(>F)"]))
cat(sprintf("  PERMDISP  p=%.4g %s\n", pd,
    ifelse(pd>0.05,"(dispersion OK)","(dispersion differs -> caution)")))
print(round(tapply(bd$distances, m$diagnosis, mean), 4))
