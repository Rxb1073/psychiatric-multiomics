library(vegan)
setwd("/rds/projects/e/elhamsak-dissertation-autism/Microbiome/pooled/unifrac")

meta <- read.delim("../metadata_pooled.tsv", row.names=1, check.names=FALSE)
meta <- meta[rownames(meta) != "#q2:types", , drop=FALSE]
meta$diagnosis <- factor(meta$diagnosis, levels=c("Control","ASD"))
meta$study <- factor(meta$study)

run <- function(dfile, label){
  d <- as.dist(as.matrix(read.delim(dfile, row.names=1, check.names=FALSE)))
  m <- meta[labels(d), , drop=FALSE]
  set.seed(1)
  pm <- adonis2(d ~ study + diagnosis, data=m, permutations=999, by="terms")
  bd <- betadisper(d, m$diagnosis)
  pd <- permutest(bd, permutations=999)$tab["Groups","Pr(>F)"]
  cat(sprintf("\n%-20s n=%d\n", label, length(labels(d))))
  cat(sprintf("  PERMANOVA study     R2=%.4f p=%.4g\n", pm["study","R2"], pm["study","Pr(>F)"]))
  cat(sprintf("  PERMANOVA diagnosis R2=%.4f p=%.4g\n", pm["diagnosis","R2"], pm["diagnosis","Pr(>F)"]))
  cat(sprintf("  PERMDISP p=%.4g %s\n", pd,
      ifelse(pd>0.05,"(dispersion OK -> location shift)","(dispersion differs -> caution)")))
  print(round(tapply(bd$distances, m$diagnosis, mean), 4))
}

cat("===== COMPLEMENTARY SEPP UniFrac (phylogenetic beta, rarefied 4030) =====\n")
run("unweighted_export/distance-matrix.tsv", "Unweighted UniFrac")
run("weighted_export/distance-matrix.tsv",   "Weighted UniFrac")
