
#!/usr/bin/env Rscript
# =====================================================================
# What is wrong with the pipeline? RNA quality and abundance.
#
# THREE CONCERNS, ONE DIAGNOSTIC
#
#  1  RIN IS NOT IN THE POOLED MODEL.
#     The primary design is ~ study + sex + ph_z + disease. RIN was
#     dropped because it is not fully populated in both cohorts (the
#     Columbia paper notes two samples without it). But the dominant
#     GSEA signal is RIBOSOME UP / ECM-COLLAGEN DOWN, which is a textbook
#     RNA-degradation signature. Both cohorts deposit RIN. This is the
#     most likely unmodelled confound in the whole arm.
#
#  2  THE RANKING CORRELATES WITH ABUNDANCE (rho = +0.189).
#     Ten times stronger than the length correlation (-0.002). Abundant
#     transcripts tend to be UP in MDD, and ribosomal proteins are among
#     the most abundant in the transcriptome. Unexplained, and the most
#     likely driver of the 1,143-term GO BP result.
#
#  3  THE COHORTS DIFFER IN DEPTH AND ASSIGNMENT RATE.
#     GSE101521 ~12M reads at 29.7% assigned; GSE102556 ~24M at 41.0%.
#     Study-as-covariate handles mean shifts, not depth-dependent
#     dispersion.
#
# WHAT THIS RUNS
#   A  how complete is RIN, and does it differ by diagnosis?
#   B  refit the pooled model WITH RIN; compare gene lists and the IEG set
#   C  does the ribosome/ECM axis survive RIN adjustment?
#   D  abundance: is the ranking's dependence on baseMean removed by RIN?
#   E  the IEG score against RIN, depth and assignment rate directly
#   F  benchmark the hand-curated IEG set against MSigDB Hallmark and GO
#
# The finding to protect is the IEG result. If it collapses under RIN
# adjustment that is a real problem; if only the ribosome/ECM axis
# collapses, that axis was RNA quality and should not be reported.
#
# READ-ONLY. Writes to MDD/audit/qualitycheck/.
# =====================================================================

suppressPackageStartupMessages({
  library(DESeq2); library(limma); library(edgeR)
  library(org.Hs.eg.db); library(AnnotationDbi) })

MDD <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD"
DE  <- file.path(MDD, "transcriptomics/DE")
P2  <- file.path(MDD, "audit/primary2")
OUT <- file.path(MDD, "audit/qualitycheck"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
hdr <- function(x) cat("\n", strrep("=",76), "\n", x, "\n", strrep("=",76), "\n", sep="")
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
set.seed(42)

IEG  <- c("NPAS4","FOS","FOSB","EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3",
          "ATF3","BTG2","ARC","JUN","JUNB","DUSP1","DUSP6")
RIBO <- "^RP[LS][0-9]"          # ribosomal proteins  -- the UP arm of the axis
COLL <- "^COL[0-9]"             # collagens           -- the DOWN arm

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
  col$study <- study
  list(cnt=cnt, col=col) }
D <- lapply(names(CO), function(n) load1(CO[[n]], n)); names(D) <- names(CO)

# ===================== A. RIN COMPLETENESS =====================
hdr("A. IS RIN USABLE?")
for (n in names(D)) {
  col <- D[[n]]$col; r <- num(col$rin)
  cat(sprintf("\n  %-10s RIN present %d/%d", n, sum(!is.na(r)), length(r)))
  if (sum(!is.na(r))) {
    cat(sprintf("   range %.1f-%.1f  median %.1f\n", min(r,na.rm=TRUE),
                max(r,na.rm=TRUE), median(r,na.rm=TRUE)))
    ok <- !is.na(r)
    a <- r[ok & col$disease=="Control"]; b <- r[ok & col$disease=="MDD"]
    if (length(a)>2 && length(b)>2) {
      p <- wilcox.test(a,b)$p.value
      cat(sprintf("             Control median %.2f vs MDD %.2f   p = %.3f%s\n",
                  median(a), median(b), p,
                  ifelse(p<0.05, "   <-- RIN DIFFERS BY DIAGNOSIS","")))
    }
  } else cat("   -- unusable\n")
  ls_ <- colSums(D[[n]]$cnt)
  cat(sprintf("             library size median %.1fM  range %.1f-%.1fM\n",
              median(ls_)/1e6, min(ls_)/1e6, max(ls_)/1e6))
  a <- ls_[col$disease=="Control"]; b <- ls_[col$disease=="MDD"]
  cat(sprintf("             depth by diagnosis: Control %.1fM vs MDD %.1fM  p = %.3f\n",
              median(a)/1e6, median(b)/1e6, wilcox.test(a,b)$p.value))
}

# assemble pooled
g   <- Reduce(intersect, lapply(D, function(d) rownames(d$cnt)))
cnt <- do.call(cbind, lapply(D, function(d) d$cnt[g,,drop=FALSE]))
col <- do.call(rbind, lapply(D, function(d)
        d$col[, c("sample_id","disease","study","sex","ph","rin")]))
rownames(col) <- col$sample_id; col <- col[colnames(cnt), ]
col$study <- factor(col$study); col$sex <- factor(tolower(col$sex))
col$ph_z  <- scale(num(col$ph))[,1]
col$rin_n <- num(col$rin)
col$depth_z <- scale(log10(colSums(cnt)))[,1]

usable_rin <- sum(!is.na(col$rin_n))
cat(sprintf("\n  POOLED: RIN present in %d/%d donors\n", usable_rin, nrow(col)))
if (usable_rin < nrow(col)) {
  cat(sprintf("  %d donors lack RIN. Two options tested below:\n", nrow(col)-usable_rin))
  cat("    (i)  drop them  (ii) mean-impute within cohort\n") }

# ===================== B. REFIT WITH RIN =====================
hdr("B. POOLED MODEL, WITH AND WITHOUT RIN")
keep <- rowSums(cnt >= 10) >= min(table(col$disease))
cat(sprintf("  %d genes pass the filter\n", sum(keep)))

sym <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
        keys=sub("\\..*","",rownames(cnt)[keep]), keytype="ENSEMBL", columns="SYMBOL"))
sym <- sym[!duplicated(sym$ENSEMBL),]
SY  <- setNames(sym$SYMBOL, sym$ENSEMBL)

fitmod <- function(cn, cl, form, lab) {
  d <- DESeq(DESeqDataSetFromMatrix(cn, cl, design=form), quiet=TRUE)
  r <- as.data.frame(results(d, contrast=c("disease","MDD","Control")))
  r$gene <- rownames(r); r$symbol <- SY[sub("\\..*","",r$gene)]
  n05 <- sum(r$padj<0.05, na.rm=TRUE)
  ne  <- sum(r$padj<0.05 & abs(r$log2FoldChange)>0.3, na.rm=TRUE)
  cat(sprintf("  %-34s n=%3d  FDR<.05 %3d  +|LFC|>0.3 %3d  %s\n",
              lab, ncol(cn), n05, ne, deparse(form)))
  r }

# baseline (current primary)
r_base <- fitmod(cnt[keep,], col, ~ study + sex + ph_z + disease, "current primary (no RIN)")

# RIN, complete cases
ok <- !is.na(col$rin_n)
col2 <- col[ok,]; col2$rin_z <- scale(col2$rin_n)[,1]
k2 <- rowSums(cnt[,ok] >= 10) >= min(table(col2$disease))
r_rin <- fitmod(cnt[k2, ok], col2, ~ study + sex + ph_z + rin_z + disease,
                "+ RIN (complete cases)")

# RIN, mean-imputed within cohort
col3 <- col
for (s in levels(col3$study)) {
  i <- col3$study==s
  m <- mean(col3$rin_n[i], na.rm=TRUE)
  col3$rin_n[i & is.na(col3$rin_n)] <- m }
col3$rin_z <- scale(col3$rin_n)[,1]
r_imp <- fitmod(cnt[keep,], col3, ~ study + sex + ph_z + rin_z + disease,
                "+ RIN (cohort mean-imputed)")

# + depth
r_dep <- fitmod(cnt[keep,], col3, ~ study + sex + ph_z + rin_z + depth_z + disease,
                "+ RIN + sequencing depth")

# ===================== C. THE AXIS AND THE IEGs =====================
hdr("C. WHAT HAPPENS TO EACH GENE FAMILY?")
famstat <- function(r, pat, lab) {
  k <- grepl(pat, r$symbol)
  if (sum(k, na.rm=TRUE) < 5) return(c(NA,NA))
  c(median(r$log2FoldChange[k], na.rm=TRUE), sum(k, na.rm=TRUE)) }
iegstat <- function(r) {
  k <- r$symbol %in% IEG
  c(median(r$log2FoldChange[k], na.rm=TRUE), sum(r$padj[k]<0.05, na.rm=TRUE)) }

MODS <- list("no RIN"=r_base, "+RIN complete"=r_rin,
             "+RIN imputed"=r_imp, "+RIN +depth"=r_dep)
cat(sprintf("  %-16s %10s %10s %12s %10s\n",
            "model","ribo LFC","coll LFC","IEG med LFC","IEG sig"))
for (n in names(MODS)) {
  r <- MODS[[n]]
  rb <- famstat(r, RIBO, "ribo"); cl_ <- famstat(r, COLL, "coll"); ie <- iegstat(r)
  cat(sprintf("  %-16s %+10.3f %+10.3f %+12.3f %10d\n", n, rb[1], cl_[1], ie[1], ie[2]))
}
cat("\n  ribosome UP + collagen DOWN is the RNA-degradation signature.\n")
cat("  If those two collapse toward zero once RIN is in the model, the axis\n")
cat("  was RNA quality. If the IEG column holds, the finding is safe.\n")

cat("\n  per-gene IEG comparison, no-RIN vs +RIN(imputed):\n")
cat(sprintf("  %-8s %10s %10s %11s %11s\n","gene","LFC noRIN","LFC +RIN","padj noRIN","padj +RIN"))
for (s in IEG) {
  i1 <- which(r_base$symbol==s)[1]; i2 <- which(r_imp$symbol==s)[1]
  if (is.na(i1)||is.na(i2)) next
  cat(sprintf("  %-8s %+10.3f %+10.3f %11s %11s\n", s,
      r_base$log2FoldChange[i1], r_imp$log2FoldChange[i2],
      ifelse(is.na(r_base$padj[i1]),"NA",sprintf("%.2e",r_base$padj[i1])),
      ifelse(is.na(r_imp$padj[i2]),"NA",sprintf("%.2e",r_imp$padj[i2]))))
}

# ===================== D. ABUNDANCE DEPENDENCE =====================
hdr("D. DOES THE ABUNDANCE DEPENDENCE SURVIVE RIN ADJUSTMENT?")
cat("  rho(ranking statistic, log10 baseMean). Was +0.189 without RIN.\n\n")
for (n in names(MODS)) {
  r <- MODS[[n]]
  rk <- sign(r$log2FoldChange) * -log10(pmax(r$pvalue, .Machine$double.xmin))
  ok2 <- is.finite(rk) & is.finite(r$baseMean)
  rho <- cor(rk[ok2], log10(r$baseMean[ok2]+1), method="spearman")
  cat(sprintf("  %-16s rho = %+.3f%s\n", n, rho,
              ifelse(abs(rho)>0.15, "   <-- still present", "")))
}
cat("\n  A drop toward zero means abundance dependence was RNA quality.\n")
cat("  If it persists, something else is driving it (normalisation,\n")
cat("  library complexity, or composition).\n")

# ===================== E. IEG SCORE vs TECHNICAL VARIABLES =====================
hdr("E. IS THE IEG SCORE ITSELF EXPLAINED BY TECHNICAL VARIABLES?")
dds <- estimateSizeFactors(DESeqDataSetFromMatrix(cnt[keep,], col, design=~disease))
nc  <- counts(dds, normalized=TRUE)
ids <- rownames(nc)[SY[sub("\\..*","",rownames(nc))] %in% IEG]
R <- t(apply(nc[ids,,drop=FALSE], 1, function(v) rank(v)/length(v)))
col$ieg <- apply(R, 2, median)
cat(sprintf("  IEG score built from %d genes\n\n", length(ids)))
for (v in c("rin_n","depth_z","ph_z")) {
  ok2 <- is.finite(col[[v]]) & is.finite(col$ieg)
  if (sum(ok2) < 10) next
  cat(sprintf("  rho(IEG score, %-8s) = %+.3f   (n=%d)\n", v,
              cor(col$ieg[ok2], col[[v]][ok2], method="spearman"), sum(ok2)))
}
ok2 <- is.finite(col$rin_n)
if (sum(ok2) > 20) {
  m <- lm(ieg ~ rin_n + depth_z + ph_z, data=col[ok2,])
  col$ieg_resid <- NA; col$ieg_resid[ok2] <- residuals(m)
  a <- col$ieg_resid[ok2 & col$disease=="Control"]
  b <- col$ieg_resid[ok2 & col$disease=="MDD"]
  auc <- wilcox.test(a,b)$statistic/(length(a)*length(b))
  cat(sprintf("\n  IEG score residualised on RIN + depth + pH:\n"))
  cat(sprintf("    Control %+.3f  MDD %+.3f   AUC %.3f   p %.4f%s\n",
              median(a), median(b), auc, wilcox.test(a,b)$p.value,
              ifelse(wilcox.test(a,b)$p.value<0.05, "   SURVIVES","")))
  cat(sprintf("    technical variables explain %.1f%% of the IEG score\n",
              100*summary(m)$r.squared))
}

# ===================== F. BENCHMARK THE GENE SET =====================
hdr("F. IS THE HAND-CURATED IEG SET DEFENSIBLE?")
cat("  The 16-gene set was curated by us. Compare against GO terms that\n")
cat("  define the same biology independently.\n\n")
GOSETS <- c("GO:0080135","GO:0034097","GO:0006357")   # regulation of cellular response to stress etc.
for (go in c("GO:0032870","GO:0071375","GO:0035690")) {
  gg <- try(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=go,
             keytype="GO", columns="SYMBOL")), silent=TRUE)
  if (inherits(gg,"try-error") || !nrow(gg)) next
  s <- unique(gg$SYMBOL); ov <- intersect(s, IEG)
  cat(sprintf("  %s: %d genes, %d overlap with our IEG set (%s)\n",
              go, length(s), length(ov), paste(head(ov,6), collapse=", ")))
}
d0 <- estimateSizeFactors(DESeqDataSetFromMatrix(cnt[keep,], col, design=~disease))
lc <- log2(counts(d0, normalized=TRUE)+1)
iegm <- lc[ids,,drop=FALSE]
cc <- cor(t(iegm))
cat(sprintf("\n  mean pairwise correlation WITHIN the IEG set: %+.3f\n",
            mean(cc[upper.tri(cc)])))
rnd <- replicate(200, { j <- sample(nrow(lc), length(ids))
  m <- cor(t(lc[j,,drop=FALSE])); mean(m[upper.tri(m)]) })
cat(sprintf("  mean for 200 random gene sets of the same size: %+.3f (95%% up to %+.3f)\n",
            mean(rnd), quantile(rnd, 0.95)))
cat("  A much higher within-set correlation confirms the set is a genuine\n")
cat("  co-regulated module, not an arbitrary collection.\n")

hdr("VERDICT")
cat("
  Read the table in section C first.

  ribosome and collagen collapse, IEG holds
      RNA quality was driving the GSEA axis. Drop the ribosome/ECM terms
      from the integration and report the IEG result unchanged.

  IEG also collapses
      the finding was partly RNA quality. That is serious and must be
      reported.

  nothing changes
      RIN was not a confound here; the ribosome/ECM axis needs another
      explanation, and section D says whether abundance is it.
")
save.image(file.path(OUT,"qualitycheck.RData"))
cat("\noutputs in", OUT, "\n")
sessionInfo()

