
#!/usr/bin/env Rscript
# =====================================================================
# Is the IEG difference a real distributional shift, or a few extreme controls?
#
# WHY THIS IS NEEDED. Removing the 4 highest-Cook's donors took the IEGs from
# 6 significant to 0 (NPAS4 -2.389 -> -0.382). But 3 of those 4 were CONTROLS
# with the highest FOS in the study. For a gene that is genuinely lower in MDD,
# the high controls ARE the contrast -- so Cook's-based removal deletes the
# samples carrying the signal and will always look devastating. It is not a
# neutral robustness test for this kind of gene.
#
# The medians tell a different story from the log2FC:
#     GSE101521   Control 565.7  MDD 360.7   ratio 1.57x
#     GSE102556   Control 240.3  MDD 206.8   ratio 1.16x
# but the pooled log2FC was -1.474 (2.8x), far larger than either median shift.
# The difference lives in the UPPER TAIL of the controls, and DESeq2's
# mean-based negative-binomial estimate is dominated by it.
#
# IEGs are activity-dependent and respond to circumstances in the final hours
# of life, so post-mortem IEG data is right-skewed by construction. The mean
# is the wrong summary. This script uses rank-based statistics instead.
#
#   1  Wilcoxon per gene, within each cohort (no distributional assumption)
#   2  the same after trimming the top 10% of controls -- the direct test of
#      "is this just a few extreme controls?"
#   3  a per-sample IEG score (median rank across the 16 genes), so each donor
#      contributes one number and no single gene or sample can dominate
#   4  quantile comparison, to show WHERE in the distribution the shift sits
#
# READ-ONLY. Writes to MDD/audit/primary2/.
# =====================================================================

suppressPackageStartupMessages({ library(DESeq2) })

MDD <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD"
DE  <- file.path(MDD, "transcriptomics/DE")
OUT <- file.path(MDD, "audit/primary2")
hdr <- function(x) cat("\n", strrep("=",76), "\n", x, "\n", strrep("=",76), "\n", sep="")
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

IEG <- c("NPAS4","FOS","FOSB","EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3",
         "ATF3","BTG2","ARC","JUN","JUNB","DUSP1","DUSP6")
# negative control: a set that should NOT differ, to calibrate the same tests
CTRLSET <- c("CALB1","CALB2","MBP","PLP1","MOG","MAG","CNP","SOX10")

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

# symbol -> ensembl present in each cohort
sym2id <- function(nc, syms) {
  m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
        keys=syms, keytype="SYMBOL", columns="ENSEMBL"))
  m <- m[!is.na(m$ENSEMBL), ]
  e <- sub("\\..*","", rownames(nc))
  out <- sapply(syms, function(s) {
    ids <- m$ENSEMBL[m$SYMBOL==s]
    hit <- rownames(nc)[e %in% ids]
    if (length(hit)) hit[1] else NA_character_ })
  out[!is.na(out)] }

hdr("0. SETUP")
for (n in names(D))
  cat(sprintf("  %-10s %d genes x %d donors (%s)\n", n, nrow(D[[n]]$nc), ncol(D[[n]]$nc),
      paste(names(table(D[[n]]$col$disease)), table(D[[n]]$col$disease), sep="=", collapse=" ")))
cat("  all tests below use DESeq2 size-factor normalised counts\n")

# ================= 1. WILCOXON PER GENE =================
hdr("1. WILCOXON PER GENE, WITHIN EACH COHORT")
cat("  Rank-based: asks whether the TYPICAL MDD brain is lower than the\n")
cat("  TYPICAL control, with no assumption about the tail.\n")

wtab <- list()
for (n in names(D)) {
  nc <- D[[n]]$nc; dz <- D[[n]]$col$disease
  ids <- sym2id(nc, IEG)
  cat(sprintf("\n  --- %s (%d/%d IEGs present) ---\n", n, length(ids), length(IEG)))
  cat(sprintf("  %-8s %10s %10s %7s %10s %10s\n",
              "gene","med Ctrl","med MDD","ratio","Wilcoxon p","BH"))
  rows <- NULL
  for (s in names(ids)) {
    v <- nc[ids[s], ]
    a <- v[dz=="Control"]; b <- v[dz=="MDD"]
    p <- suppressWarnings(wilcox.test(a, b)$p.value)
    rows <- rbind(rows, data.frame(cohort=n, gene=s, med_ctrl=median(a), med_mdd=median(b),
                                   ratio=median(a)/median(b), p=p, stringsAsFactors=FALSE))
  }
  rows$bh <- p.adjust(rows$p, "BH"); rows <- rows[order(rows$p), ]
  for (i in seq_len(nrow(rows)))
    cat(sprintf("  %-8s %10.1f %10.1f %7.2f %10.2e %10.2e%s\n", rows$gene[i],
                rows$med_ctrl[i], rows$med_mdd[i], rows$ratio[i], rows$p[i], rows$bh[i],
                ifelse(rows$bh[i]<0.05, "  *", "")))
  cat(sprintf("    %d/%d at BH<0.05;  %d/%d with control median higher\n",
              sum(rows$bh<0.05), nrow(rows), sum(rows$ratio>1), nrow(rows)))
  wtab[[n]] <- rows
}
write.table(do.call(rbind, wtab), file.path(OUT,"IEG_wilcoxon_per_gene.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ================= 2. TRIM THE CONTROL TAIL =================
hdr("2. THE DECISIVE TEST: TRIM THE TOP 10% OF CONTROLS")
cat("  If the effect is only a few extreme controls, removing the highest\n")
cat("  controls per gene should abolish it. If the whole distribution is\n")
cat("  shifted, it should survive with a smaller effect.\n")

for (n in names(D)) {
  nc <- D[[n]]$nc; dz <- D[[n]]$col$disease
  ids <- sym2id(nc, IEG)
  cat(sprintf("\n  --- %s ---\n", n))
  cat(sprintf("  %-8s %11s %11s %11s %11s\n","gene","p full","p trimmed","ratio full","ratio trim"))
  nsurv <- 0
  for (s in names(ids)) {
    v <- nc[ids[s], ]; a <- v[dz=="Control"]; b <- v[dz=="MDD"]
    p1 <- suppressWarnings(wilcox.test(a,b)$p.value)
    cut <- quantile(a, 0.90)
    a2  <- a[a <= cut]
    p2  <- suppressWarnings(wilcox.test(a2,b)$p.value)
    if (!is.na(p2) && p2 < 0.05) nsurv <- nsurv + 1
    cat(sprintf("  %-8s %11.2e %11.2e %11.2f %11.2f%s\n", s, p1, p2,
                median(a)/median(b), median(a2)/median(b),
                ifelse(!is.na(p2) && p2<0.05, "  *", "")))
  }
  cat(sprintf("    %d/%d still nominally significant after trimming the top 10%% of controls\n",
              nsurv, length(ids)))
}

# ================= 3. PER-SAMPLE IEG SCORE =================
hdr("3. PER-SAMPLE IEG SCORE (one number per donor)")
cat("  Each donor's median percentile rank across the 16 IEGs. No single gene\n")
cat("  or sample can dominate, and it is scale-free.\n")

scores <- NULL
for (n in names(D)) {
  nc <- D[[n]]$nc; dz <- D[[n]]$col$disease
  ids <- sym2id(nc, IEG)
  R <- t(apply(nc[ids, , drop=FALSE], 1, function(v) rank(v)/length(v)))
  sc <- apply(R, 2, median)
  a <- sc[dz=="Control"]; b <- sc[dz=="MDD"]
  p <- suppressWarnings(wilcox.test(a,b)$p.value)
  auc <- suppressWarnings(wilcox.test(a,b)$statistic) / (length(a)*length(b))
  cat(sprintf("\n  %-10s Control median %.3f (n=%d)   MDD median %.3f (n=%d)\n",
              n, median(a), length(a), median(b), length(b)))
  cat(sprintf("             Wilcoxon p = %.4f   AUC = %.3f  (0.5 = no difference)\n", p, auc))
  scores <- rbind(scores, data.frame(cohort=n, sample=names(sc), score=sc,
                                     disease=as.character(dz), stringsAsFactors=FALSE))
}
write.table(scores, file.path(OUT,"IEG_sample_scores.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

cat("\n  --- both cohorts combined (score is cohort-internal, so poolable) ---\n")
a <- scores$score[scores$disease=="Control"]; b <- scores$score[scores$disease=="MDD"]
p <- wilcox.test(a,b)$p.value
auc <- wilcox.test(a,b)$statistic/(length(a)*length(b))
cat(sprintf("  Control median %.3f (n=%d)   MDD median %.3f (n=%d)\n",
            median(a), length(a), median(b), length(b)))
cat(sprintf("  Wilcoxon p = %.4f   AUC = %.3f\n", p, auc))

# negative-control set, same score, for calibration
cat("\n  --- same score built from the NEGATIVE CONTROL genes ---\n")
ns <- NULL
for (n in names(D)) {
  nc <- D[[n]]$nc; dz <- D[[n]]$col$disease
  ids <- sym2id(nc, CTRLSET); if (length(ids) < 3) next
  R <- t(apply(nc[ids, , drop=FALSE], 1, function(v) rank(v)/length(v)))
  ns <- rbind(ns, data.frame(score=apply(R,2,median), disease=as.character(dz)))
}
if (!is.null(ns)) {
  a2 <- ns$score[ns$disease=="Control"]; b2 <- ns$score[ns$disease=="MDD"]
  cat(sprintf("  Control median %.3f   MDD median %.3f   p = %.4f   AUC = %.3f\n",
              median(a2), median(b2), wilcox.test(a2,b2)$p.value,
              wilcox.test(a2,b2)$statistic/(length(a2)*length(b2))))
  cat("  (this should be non-significant with AUC near 0.5)\n")
}

# ================= 4. WHERE IS THE SHIFT? =================
hdr("4. WHERE IN THE DISTRIBUTION DOES THE SHIFT SIT?")
cat("  If only the top tail differs, the lower quantiles will match.\n")
for (n in names(D)) {
  nc <- D[[n]]$nc; dz <- D[[n]]$col$disease
  ids <- sym2id(nc, IEG)
  cat(sprintf("\n  --- %s: quantiles of the per-sample IEG score ---\n", n))
  R <- t(apply(nc[ids, , drop=FALSE], 1, function(v) rank(v)/length(v)))
  sc <- apply(R, 2, median)
  qs <- c(0.10,0.25,0.50,0.75,0.90)
  a <- quantile(sc[dz=="Control"], qs); b <- quantile(sc[dz=="MDD"], qs)
  cat(sprintf("  %-10s %8s %8s %8s %8s %8s\n","quantile","10%","25%","50%","75%","90%"))
  cat(sprintf("  %-10s %8.3f %8.3f %8.3f %8.3f %8.3f\n","Control", a[1],a[2],a[3],a[4],a[5]))
  cat(sprintf("  %-10s %8.3f %8.3f %8.3f %8.3f %8.3f\n","MDD",     b[1],b[2],b[3],b[4],b[5]))
  cat(sprintf("  %-10s %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f\n","diff",
              a[1]-b[1],a[2]-b[2],a[3]-b[3],a[4]-b[4],a[5]-b[5]))
}

hdr("VERDICT")
cat("
  Positive differences at EVERY quantile
      the whole distribution is shifted, not just the tail. The Cook's result
      then reflects DESeq2 protecting against influential points on a skewed
      gene, not a finding carried by four donors. Report medians and the
      rank-based p alongside the log2FC, and note the mean-based estimate is
      tail-influenced.

  Differences ONLY at the 75th/90th percentile
      the effect is confined to high-IEG controls. That is still a real group
      difference, but it should be described as 'a subset of controls show high
      immediate-early expression not seen in MDD' rather than as a uniform
      reduction.

  Nothing significant once the control tail is trimmed
      the effect is a handful of extreme controls and must be reported as such.

  The negative-control score is the calibration: it should be non-significant
  with AUC near 0.5 under all three readings.
")
cat("\noutputs in", OUT, "\n")
sessionInfo()

