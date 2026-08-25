
#!/usr/bin/env Rscript
# =====================================================================
# Two remaining questions about the IEG score.
#
# WHERE THINGS STAND
#   per-sample IEG score:  GSE101521 p=0.0126 AUC 0.690
#                          GSE102556 p=0.0366 AUC 0.677
#                          combined  p=0.0013 AUC 0.681
#   all 10 quantile comparisons positive (binomial p=0.002)
#   15/16 genes control-median-higher in BOTH cohorts (p=2.7e-07 jointly)
#
# TWO PROBLEMS
#   1  Trimming the top 10% of controls left only 4/16 (GSE101521) and 0/16
#      (GSE102556) genes nominally significant. But that trim was run PER
#      GENE, where the tests were already underpowered, not on the score.
#      The like-for-like test was never done.
#
#   2  The negative-control set (calbindin + oligodendrocyte) came out
#      HIGHER in MDD at p=0.0505, AUC 0.390 -- the mirror image of the IEG
#      result. That is the neuron-down / glia-up axis, and a control moving
#      the opposite way at p=0.05 is not clean.
#
# THIS SCRIPT
#   A  trim the top 10%, 20% and 30% of controls and re-test THE SCORE
#   B  residualise the IEG score on the glia score, per cohort, and re-test.
#      If an IEG-specific effect survives, it is not simply the global axis.
#   C  the two scores' correlation, and the glia score's own behaviour
#   D  a permutation null for the score, to check the Wilcoxon p is honest
#
# READ-ONLY. Writes to MDD/audit/primary2/.
# =====================================================================

suppressPackageStartupMessages({ library(DESeq2) })

MDD <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD"
DE  <- file.path(MDD, "transcriptomics/DE")
OUT <- file.path(MDD, "audit/primary2")
hdr <- function(x) cat("\n", strrep("=",76), "\n", x, "\n", strrep("=",76), "\n", sep="")
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
set.seed(42)

IEG  <- c("NPAS4","FOS","FOSB","EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3",
          "ATF3","BTG2","ARC","JUN","JUNB","DUSP1","DUSP6")
GLIA <- c("MBP","PLP1","MOG","MAG","CNP","SOX10","OLIG1","OLIG2","MOBP","CLDN11")
CALB <- c("CALB1","CALB2")

CO <- list(
  GSE101521 = list(cnt=file.path(MDD,"transcriptomics/GSE101521/counts/GSE101521_gene_counts_matrix.tsv"),
                   col=file.path(DE,"GSE101521_colData.tsv")),
  GSE102556 = list(cnt=file.path(MDD,"transcriptomics/GSE102556/counts/GSE102556_gene_counts_matrix.tsv"),
                   col=file.path(DE,"GSE102556_colData.tsv")))

norm1 <- function(x) {
  cnt <- read.delim(x$cnt, row.names=1, check.names=FALSE)
  col <- read.delim(x$col, check.names=FALSE); names(col)[1] <- "sample_id"
  rownames(col) <- col$sample_id
  k <- intersect(colnames(cnt), rownames(col))
  cnt <- round(as.matrix(cnt[,k,drop=FALSE])); col <- col[k,,drop=FALSE]
  col$disease <- factor(ifelse(grepl("MDD|depress", col$disease, ignore.case=TRUE),
                               "MDD","Control"), levels=c("Control","MDD"))
  keep <- rowSums(cnt >= 10) >= min(table(col$disease))
  dds  <- estimateSizeFactors(DESeqDataSetFromMatrix(cnt[keep,], col, design=~disease))
  list(nc=counts(dds, normalized=TRUE), col=col) }
D <- lapply(names(CO), function(n) norm1(CO[[n]])); names(D) <- names(CO)

sym2id <- function(nc, syms) {
  m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
        keys=syms, keytype="SYMBOL", columns="ENSEMBL"))
  m <- m[!is.na(m$ENSEMBL), ]; e <- sub("\\..*","", rownames(nc))
  out <- sapply(syms, function(s) { h <- rownames(nc)[e %in% m$ENSEMBL[m$SYMBOL==s]]
                                    if (length(h)) h[1] else NA_character_ })
  out[!is.na(out)] }

# per-sample median percentile rank across a gene set
score <- function(nc, syms) {
  ids <- sym2id(nc, syms); if (length(ids) < 2) return(NULL)
  R <- t(apply(nc[ids, , drop=FALSE], 1, function(v) rank(v)/length(v)))
  apply(R, 2, median) }

auc <- function(a, b) suppressWarnings(wilcox.test(a,b)$statistic)/(length(a)*length(b))
wp  <- function(a, b) suppressWarnings(wilcox.test(a,b)$p.value)

hdr("0. SCORES")
S <- list()
for (n in names(D)) {
  nc <- D[[n]]$nc; dz <- D[[n]]$col$disease
  S[[n]] <- data.frame(sample=colnames(nc), disease=as.character(dz),
                       ieg=score(nc, IEG), glia=score(nc, GLIA),
                       calb=score(nc, CALB), stringsAsFactors=FALSE)
  cat(sprintf("  %-10s n=%d  (IEG %d genes, glia %d genes)\n", n, ncol(nc),
              length(sym2id(nc,IEG)), length(sym2id(nc,GLIA))))
}

# ================== A. TRIM THE SCORE ==================
hdr("A. TRIM THE TOP CONTROLS -- ON THE SCORE, NOT PER GENE")
cat("  The per-gene trim was underpowered. This is the like-for-like test.\n")
cat("  Removing the highest-IEG controls should abolish the effect if it is\n")
cat("  a tail phenomenon; it should survive, attenuated, if the whole\n")
cat("  distribution is shifted.\n")

for (n in names(D)) {
  d <- S[[n]]; a0 <- d$ieg[d$disease=="Control"]; b <- d$ieg[d$disease=="MDD"]
  cat(sprintf("\n  --- %s ---\n", n))
  cat(sprintf("  %-22s %5s %10s %8s %10s\n","controls kept","n","median C","AUC","p"))
  for (tr in c(0, 0.10, 0.20, 0.30)) {
    a <- if (tr == 0) a0 else a0[a0 <= quantile(a0, 1-tr)]
    cat(sprintf("  %-22s %5d %10.3f %8.3f %10.4f%s\n",
                sprintf("top %.0f%% removed", 100*tr), length(a), median(a),
                auc(a,b), wp(a,b), ifelse(wp(a,b)<0.05,"  *","")))
  }
}
cat("\n  --- both cohorts combined ---\n")
allS <- do.call(rbind, S)
a0 <- allS$ieg[allS$disease=="Control"]; b <- allS$ieg[allS$disease=="MDD"]
cat(sprintf("  %-22s %5s %10s %8s %10s\n","controls kept","n","median C","AUC","p"))
for (tr in c(0, 0.10, 0.20, 0.30)) {
  a <- if (tr == 0) a0 else a0[a0 <= quantile(a0, 1-tr)]
  cat(sprintf("  %-22s %5d %10.3f %8.3f %10.4f%s\n",
              sprintf("top %.0f%% removed", 100*tr), length(a), median(a),
              auc(a,b), wp(a,b), ifelse(wp(a,b)<0.05,"  *","")))
}

# ================== B. RESIDUALISE ON GLIA ==================
hdr("B. IS THE IEG EFFECT SPECIFIC, OR THE GLOBAL NEURON/GLIA AXIS?")
cat("  The negative-control set was HIGHER in MDD (p=0.0505, AUC 0.390).\n")
cat("  Regress the IEG score on the glia score and re-test the residual.\n")
cat("  A surviving effect is IEG-specific; a vanishing one was the axis.\n")

for (n in c(names(D), "COMBINED")) {
  d <- if (n=="COMBINED") allS else S[[n]]
  cat(sprintf("\n  --- %s (n=%d) ---\n", n, nrow(d)))
  ok <- is.finite(d$ieg) & is.finite(d$glia)
  r  <- cor(d$ieg[ok], d$glia[ok], method="spearman")
  cat(sprintf("    cor(IEG score, glia score) = %+.3f\n", r))

  for (v in c("ieg","glia","calb")) {
    if (!v %in% names(d) || all(is.na(d[[v]]))) next
    a <- d[[v]][d$disease=="Control"]; b2 <- d[[v]][d$disease=="MDD"]
    cat(sprintf("    %-16s Control %.3f  MDD %.3f   AUC %.3f   p %.4f%s\n",
                paste0(v," score"), median(a), median(b2), auc(a,b2), wp(a,b2),
                ifelse(wp(a,b2)<0.05,"  *","")))
  }
  fit <- lm(ieg ~ glia, data=d[ok,])
  res <- residuals(fit)
  dz  <- d$disease[ok]
  a <- res[dz=="Control"]; b2 <- res[dz=="MDD"]
  cat(sprintf("    %-16s Control %+.3f MDD %+.3f  AUC %.3f   p %.4f%s\n",
              "IEG | glia", median(a), median(b2), auc(a,b2), wp(a,b2),
              ifelse(wp(a,b2)<0.05,"  *  SURVIVES","")))
  cat(sprintf("    glia explains %.1f%% of the variance in the IEG score\n",
              100*summary(fit)$r.squared))
}

# ================== C. PERMUTATION NULL ==================
hdr("C. PERMUTATION NULL FOR THE SCORE")
cat("  Shuffle the diagnosis labels 5000 times, within cohort, and recompute\n")
cat("  the combined AUC. Confirms the Wilcoxon p without distributional\n")
cat("  assumptions and accounts for the score being cohort-internal.\n")

obs <- auc(allS$ieg[allS$disease=="Control"], allS$ieg[allS$disease=="MDD"])
NP <- 5000
null <- replicate(NP, {
  perm <- do.call(rbind, lapply(S, function(d) { d$disease <- sample(d$disease); d }))
  auc(perm$ieg[perm$disease=="Control"], perm$ieg[perm$disease=="MDD"]) })
pperm <- (sum(abs(null-0.5) >= abs(obs-0.5)) + 1)/(NP+1)
cat(sprintf("\n  observed AUC          %.3f\n", obs))
cat(sprintf("  null AUC              median %.3f, 2.5-97.5%% [%.3f, %.3f]\n",
            median(null), quantile(null,0.025), quantile(null,0.975)))
cat(sprintf("  permutation p         %.4f  (Wilcoxon gave 0.0013)\n", pperm))

# same for the residualised score
ok <- is.finite(allS$ieg) & is.finite(allS$glia)
allS$resid <- NA; allS$resid[ok] <- residuals(lm(ieg ~ glia, data=allS[ok,]))
obs_r <- auc(allS$resid[allS$disease=="Control" & ok], allS$resid[allS$disease=="MDD" & ok])
null_r <- replicate(NP, { g <- sample(allS$disease[ok])
  auc(allS$resid[ok][g=="Control"], allS$resid[ok][g=="MDD"]) })
pr <- (sum(abs(null_r-0.5) >= abs(obs_r-0.5)) + 1)/(NP+1)
cat(sprintf("\n  residualised score: observed AUC %.3f, permutation p %.4f\n", obs_r, pr))

write.table(allS, file.path(OUT,"IEG_scores_with_glia.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

hdr("VERDICT")
cat(sprintf("
  raw IEG score        AUC %.3f
  residualised on glia AUC %.3f, permutation p %.4f

  Residualised effect still significant
      the IEG reduction is specific, not a by-product of the neuron/glia
      axis. Report it, and report the glia score alongside as the honest
      caveat -- both move, but IEG moves independently of glia.

  Residualised effect gone
      the IEG score and the glia score are two readings of one compositional
      axis. The finding must then be described as a neuronal-vs-glial shift,
      not as immediate-early gene regulation.

  Effect surviving 20-30%% control trimming
      settles the tail question: it is a distributional shift.
", obs, obs_r, pr))
cat("\noutputs in", OUT, "\n")
sessionInfo()

