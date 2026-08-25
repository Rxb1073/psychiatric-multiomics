
#!/usr/bin/env Rscript
# =====================================================================
# THE OLIGODENDROCYTE / MYELINATION MODULE
#
# WHY THIS IS BEING TESTED PROPERLY NOW
#
# The same signal has appeared three times and was each time treated as a
# control misbehaving rather than as a result:
#
#   CAMERA          oligodendrocyte "negative control" trending Up
#                   (GSE102556 female p=0.042; GSE101521 nominal 3.7e-03)
#   glia score      higher in MDD, AUC 0.398, p = 0.069
#   ORA up-list     ensheathment of neurons / axon ensheathment 1.6e-12,
#                   myelination 1.8e-11, oligodendrocyte differentiation
#                   1.1e-04 -- six coherent terms, from 129 genes
#
# Oligodendrocyte and myelin abnormalities are among the more replicated
# post-mortem findings in depression, so this is plausibly a real second
# result rather than a technical gradient.
#
# IT ALSO MEANS the oligodendrocyte set is a POOR negative control for this
# disorder -- a control must be something the disease should not affect.
# This script therefore also tests replacement controls.
#
# THE STANDARD APPLIED IS THE ONE USED FOR IEG -- nothing less:
#   1  pre-specified sets, not derived from these data
#   2  CAMERA with inter-gene correlation estimated
#   3  per-donor rank score, each cohort separately and combined
#   4  permutation null
#   5  leave-one-gene-out
#   6  new negative controls that are defensible for MDD
#   7  independence from the IEG result -- they may be one axis
#
# POINT 7 IS THE CRITICAL ONE. Neurons down and oligodendrocytes up can be
# a single compositional shift rather than two findings. The residualisation
# in section 6 decides that, and if the two are not separable this must be
# reported as ONE result about cell composition, not two.
#
# READ-ONLY. Writes to MDD/audit/oligo/.
# =====================================================================

suppressPackageStartupMessages({
  library(DESeq2); library(limma); library(edgeR)
  library(org.Hs.eg.db); library(AnnotationDbi) })

MDD <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD"
DE  <- file.path(MDD,"transcriptomics/DE")
OUT <- file.path(MDD,"audit/oligo"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
hdr <- function(x) cat("\n", strrep("=",76), "\n", x, "\n", strrep("=",76), "\n", sep="")
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
set.seed(42)

SETS <- list(
  # the module under test, split into its two components
  MYELIN      = c("MBP","PLP1","MOG","MAG","CNP","MOBP","CLDN11","MAL","OPALIN",
                  "UGT8","ASPA","ERMN","GJB1","MYRF"),
  OLIG_IDENT  = c("SOX10","OLIG1","OLIG2","NKX2-2","ST18","PRR5L","CNTN2"),
  # the established finding, for comparison and for the independence test
  IEG         = c("NPAS4","FOS","FOSB","EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3",
                  "ATF3","BTG2","ARC","JUN","JUNB","DUSP1","DUSP6"),
  # replacement negative controls: cell types with no strong prior in MDD
  NEG_endothelial = c("CLDN5","FLT1","PECAM1","VWF","CDH5","SLC2A1","ABCB1"),
  NEG_microglia   = c("P2RY12","CX3CR1","TMEM119","CSF1R","AIF1","ITGAM","TYROBP"),
  NEG_housekeep   = c("ACTB","GAPDH","TUBB","PGK1","B2M","PPIA","RPLP0","TBP","SDHA","YWHAZ"),
  # context
  ASTROCYTE   = c("GFAP","AQP4","SLC1A2","SLC1A3","ALDH1L1","S100B"),
  EXCITATORY  = c("SLC17A7","SLC17A6","SATB2","TBR1","CUX2","RORB"))

CO <- list(
  GSE101521 = list(cnt=file.path(MDD,"transcriptomics/GSE101521/counts/GSE101521_gene_counts_matrix.tsv"),
                   col=file.path(DE,"GSE101521_colData.tsv")),
  GSE102556 = list(cnt=file.path(MDD,"transcriptomics/GSE102556/counts/GSE102556_gene_counts_matrix.tsv"),
                   col=file.path(DE,"GSE102556_colData.tsv")))

load1 <- function(x, study) {
  cnt <- read.delim(x$cnt, row.names=1, check.names=FALSE)
  col <- read.delim(x$col, check.names=FALSE); names(col)[1] <- "sample_id"
  rownames(col) <- col$sample_id
  k <- intersect(colnames(cnt), rownames(col))
  cnt <- round(as.matrix(cnt[,k,drop=FALSE])); col <- col[k,,drop=FALSE]
  col$disease <- factor(ifelse(grepl("MDD|depress", col$disease, ignore.case=TRUE),
                               "MDD","Control"), levels=c("Control","MDD"))
  col$sex <- factor(tolower(col$sex)); col$study <- study
  keep <- rowSums(cnt >= 10) >= min(table(col$disease))
  cnt <- cnt[keep,]
  m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=sub("\\..*","",rownames(cnt)),
        keytype="ENSEMBL", columns="SYMBOL"))
  m <- m[!duplicated(m$ENSEMBL),]
  dds <- estimateSizeFactors(DESeqDataSetFromMatrix(cnt, col, design=~disease))
  list(cnt=cnt, nc=counts(dds, normalized=TRUE), col=col,
       sy=setNames(m$SYMBOL, m$ENSEMBL)[sub("\\..*","",rownames(cnt))]) }
D <- lapply(names(CO), function(n) load1(CO[[n]], n)); names(D) <- names(CO)

score <- function(d, syms) {
  ids <- rownames(d$nc)[which(d$sy %in% syms)]
  if (length(ids) < 3) return(NULL)
  R <- t(apply(d$nc[ids,,drop=FALSE], 1, function(v) rank(v)/length(v)))
  apply(R, 2, median) }
auc <- function(a,b) suppressWarnings(wilcox.test(a,b)$statistic)/(length(a)*length(b))
wp  <- function(a,b) suppressWarnings(wilcox.test(a,b)$p.value)

hdr("0. SETS")
for (n in names(SETS)) {
  present <- sapply(D, function(d) sum(d$sy %in% SETS[[n]]))
  cat(sprintf("  %-16s %2d symbols; present %s\n", n, length(SETS[[n]]),
              paste(names(present), present, sep="=", collapse=" ")))
}

# ================= 1. CAMERA, POOLED =================
hdr("1. CAMERA, pooled two-cohort")
g   <- Reduce(intersect, lapply(D, function(d) rownames(d$cnt)))
cnt <- do.call(cbind, lapply(D, function(d) d$cnt[g,,drop=FALSE]))
col <- do.call(rbind, lapply(D, function(d) d$col[,c("sample_id","disease","study","sex","ph")]))
rownames(col) <- col$sample_id; col <- col[colnames(cnt),]
col$study <- factor(col$study); col$ph_z <- scale(num(col$ph))[,1]
des <- model.matrix(~ study + sex + ph_z + disease, data=col)
v   <- voom(calcNormFactors(DGEList(cnt)), des)
m   <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=sub("\\..*","",rownames(v)),
        keytype="ENSEMBL", columns="SYMBOL")); m <- m[!duplicated(m$ENSEMBL),]
sy  <- setNames(m$SYMBOL, m$ENSEMBL)[sub("\\..*","",rownames(v))]
idx <- lapply(SETS, function(s) which(sy %in% s)); idx <- idx[sapply(idx,length) >= 3]
cu <- camera(v, idx, des, contrast=tail(colnames(des),1), inter.gene.cor=0)
cc <- camera(v, idx, des, contrast=tail(colnames(des),1), inter.gene.cor=NA)
cu <- cu[order(rownames(cu)),]; cc <- cc[rownames(cu),]
cat(sprintf("  %-16s %5s  %-6s %11s  %-6s %11s %8s\n",
            "set","genes","dir","p_uncorr","dir","p_CAMERA","r_gene"))
for (i in seq_len(nrow(cu)))
  cat(sprintf("  %-16s %5d  %-6s %11.2e  %-6s %11.2e %+8.3f%s\n",
      rownames(cu)[i], cu$NGenes[i], cu$Direction[i], cu$PValue[i],
      cc$Direction[i], cc$PValue[i], cc$Correlation[i],
      ifelse(cc$PValue[i] < 0.05, "  *", "")))
write.csv(data.frame(set=rownames(cc), cc), file.path(OUT,"camera_oligo.csv"), row.names=FALSE)

# ================= 2. SCORES =================
hdr("2. PER-DONOR SCORES, EACH COHORT AND COMBINED")
S <- do.call(rbind, lapply(names(D), function(n) {
  d <- D[[n]]
  df <- data.frame(cohort=n, sample=colnames(d$nc),
                   disease=as.character(d$col$disease), stringsAsFactors=FALSE)
  for (s in names(SETS)) { x <- score(d, SETS[[s]]); df[[s]] <- if (is.null(x)) NA else x }
  df }))
write.table(S, file.path(OUT,"scores.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

cat(sprintf("  %-16s %-11s %5s %8s %8s %10s\n","set","cohort","n","Ctl","MDD","AUC (p)"))
for (s in names(SETS)) {
  if (all(is.na(S[[s]]))) next
  for (ch in c(names(D), "COMBINED")) {
    d <- if (ch=="COMBINED") S else S[S$cohort==ch,]
    a <- d[[s]][d$disease=="Control"]; b <- d[[s]][d$disease=="MDD"]
    a <- a[!is.na(a)]; b <- b[!is.na(b)]; if (length(a)<4||length(b)<4) next
    cat(sprintf("  %-16s %-11s %5d %8.3f %8.3f %6.3f (%.4f)%s\n", s, ch, length(a)+length(b),
                median(a), median(b), auc(a,b), wp(a,b),
                ifelse(wp(a,b)<0.05,"  *","")))
  }
}

# ================= 3. PERMUTATION =================
hdr("3. PERMUTATION NULL (5000 within-cohort shuffles)")
NP <- 5000
for (s in c("MYELIN","OLIG_IDENT","IEG")) {
  if (all(is.na(S[[s]]))) next
  obs <- auc(S[[s]][S$disease=="Control"], S[[s]][S$disease=="MDD"])
  null <- replicate(NP, {
    p <- do.call(rbind, lapply(split(S, S$cohort), function(d){d$disease <- sample(d$disease); d}))
    auc(p[[s]][p$disease=="Control"], p[[s]][p$disease=="MDD"]) })
  pp <- (sum(abs(null-0.5) >= abs(obs-0.5))+1)/(NP+1)
  cat(sprintf("  %-12s observed AUC %.3f   permutation p = %.4f%s\n", s, obs, pp,
              ifelse(pp<0.05,"  *","")))
}

# ================= 4. LEAVE ONE GENE OUT =================
hdr("4. LEAVE-ONE-GENE-OUT (myelin module)")
ids <- which(sy %in% SETS$MYELIN)
if (length(ids) >= 4) {
  base <- camera(v, list(M=ids), des, contrast=tail(colnames(des),1), inter.gene.cor=NA)$PValue
  cat(sprintf("  full set p = %.3e (%d genes)\n\n", base, length(ids)))
  r <- NULL
  for (i in ids) {
    p <- camera(v, list(M=setdiff(ids,i)), des, contrast=tail(colnames(des),1),
                inter.gene.cor=NA)$PValue
    r <- rbind(r, data.frame(dropped=sy[i], p=p)) }
  r <- r[order(-r$p),]
  for (i in seq_len(nrow(r)))
    cat(sprintf("  drop %-9s p = %.3e  (%.2fx)%s\n", r$dropped[i], r$p[i], r$p[i]/base,
                ifelse(r$p[i] >= 0.05, "   <-- LOSES SIGNIFICANCE","")))
}

# ================= 5. NEW NEGATIVE CONTROLS =================
hdr("5. ARE THE REPLACEMENT CONTROLS BEHAVING?")
cat("  The oligodendrocyte set cannot serve as a negative control for MDD if\n")
cat("  it is itself affected. Endothelial, microglial and housekeeping sets\n")
cat("  have no strong prior in depression and should be null.\n\n")
for (s in c("NEG_endothelial","NEG_microglia","NEG_housekeep")) {
  if (all(is.na(S[[s]]))) { cat(sprintf("  %-16s not enough genes\n", s)); next }
  a <- S[[s]][S$disease=="Control"]; b <- S[[s]][S$disease=="MDD"]
  cat(sprintf("  %-16s AUC %.3f  p %.4f%s\n", s, auc(a,b), wp(a,b),
              ifelse(wp(a,b)<0.05, "   *** FIRES -- not a usable control", "   ok")))
}

# ================= 6. ONE AXIS OR TWO? =================
hdr("6. IS THE MYELIN RESULT INDEPENDENT OF THE IEG RESULT?")
cat("  Neurons down and oligodendrocytes up can be a SINGLE compositional\n")
cat("  shift. If each survives adjustment for the other they are two findings;\n")
cat("  if not, this is one result about cell composition and must be reported\n")
cat("  as such.\n\n")
ok <- is.finite(S$MYELIN) & is.finite(S$IEG)
cat(sprintf("  cor(myelin, IEG) = %+.3f\n\n", cor(S$MYELIN[ok], S$IEG[ok], method="spearman")))
for (pair in list(c("MYELIN","IEG"), c("IEG","MYELIN"))) {
  y <- pair[1]; x <- pair[2]
  m1 <- lm(as.formula(paste(y,"~ cohort +",x)), data=S[ok,])
  rr <- residuals(m1); dz <- S$disease[ok]
  a <- rr[dz=="Control"]; b <- rr[dz=="MDD"]
  cat(sprintf("  %-7s adjusted for %-7s : AUC %.3f  p %.4f  (%s explains %.1f%%)%s\n",
      y, x, auc(a,b), wp(a,b), x, 100*summary(m1)$r.squared,
      ifelse(wp(a,b)<0.05, "   SURVIVES", "   does not survive")))
}

hdr("VERDICT")
cat("
  Both survive mutual adjustment
      two separable findings: reduced activity-dependent neuronal
      transcription AND increased oligodendrocyte/myelin gene expression.
      Report both, and note the oligodendrocyte set must be replaced as a
      negative control.

  Only one survives
      that one is primary and the other is its shadow.

  Neither survives
      a single compositional axis -- neuron-vs-glia. Report it that way,
      not as two independent modules. This would also mean the IEG result
      is better described as compositional than as activity-dependent
      transcription, which changes the interpretation substantially.
")
cat("\noutputs in", OUT, "\n")
sessionInfo()

