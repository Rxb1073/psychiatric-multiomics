
#!/usr/bin/env Rscript
# =====================================================================
# Why is FOS padj = NA in the pooled two-cohort model?
#
# FOS had log2FC -1.474 in the pooled analysis but padj = NA. DESeq2
# returns NA for two reasons:
#
#   (a) Cook's distance outlier filtering -- one sample has an extreme
#       count that dominates the fit, so DESeq2 refuses to report a p-value
#   (b) independent filtering -- low mean expression, which is implausible
#       here since FOS is abundant
#
# This matters because FOS is the single most down-regulated gene in the
# GSEA ranking and one of the headline IEGs. If a single donor is driving
# it, that should be known before it goes in a chapter.
#
# Prints the per-donor normalised counts so the outlier is visible, then
# re-runs FOS three ways: default, with Cook's filtering off, and with the
# outlier donor removed.
#
# READ-ONLY. Writes to MDD/audit/primary2/.
# =====================================================================

suppressPackageStartupMessages({ library(DESeq2) })

MDD <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD"
DE  <- file.path(MDD, "transcriptomics/DE")
OUT <- file.path(MDD, "audit/primary2")
hdr <- function(x) cat("\n", strrep("=",76), "\n", x, "\n", strrep("=",76), "\n", sep="")
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

IEG <- c("FOS","FOSB","NPAS4","EGR1","EGR2","NR4A1","NR4A3","ATF3","JUN","JUNB",
         "ARC","BTG2","DUSP1","DUSP6","EGR3","NR4A2")

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

g   <- Reduce(intersect, lapply(D, function(d) rownames(d$cnt)))
cnt <- do.call(cbind, lapply(D, function(d) d$cnt[g,,drop=FALSE]))
col <- do.call(rbind, lapply(D, function(d)
        d$col[, c("sample_id","disease","study","sex","ph")]))
rownames(col) <- col$sample_id; col <- col[colnames(cnt), ]
col$study <- factor(col$study); col$sex <- factor(tolower(col$sex))
col$ph_z  <- scale(num(col$ph))[,1]
keep <- rowSums(cnt >= 10) >= min(table(col$disease))
form <- ~ study + sex + ph_z + disease

hdr("1. FIT THE POOLED MODEL")
cat(sprintf("  %d genes x %d donors\n", sum(keep), ncol(cnt)))
dds <- DESeq(DESeqDataSetFromMatrix(cnt[keep,], col, design=form), quiet=TRUE)

# map symbols
map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
        keys=sub("\\..*","",rownames(dds)), keytype="ENSEMBL", columns="SYMBOL"))
map <- map[!duplicated(map$ENSEMBL),]
sym <- setNames(map$SYMBOL, map$ENSEMBL)
ens <- sub("\\..*","", rownames(dds))
gid <- function(s) rownames(dds)[which(ens %in% names(sym)[which(sym==s)])][1]

# ================== 2. WHY NA? ==================
hdr("2. WHICH IEGs HAVE padj = NA, AND WHY?")
res  <- results(dds)
cooks <- assays(dds)[["cooks"]]
p <- ncol(attr(dds,"modelMatrix")); n <- ncol(dds)
thr <- qf(0.99, p, n - p)
cat(sprintf("  Cook's cut-off for this design: F(0.99, %d, %d) = %.2f\n\n", p, n-p, thr))
cat(sprintf("  %-8s %10s %9s %10s %10s %10s  %s\n",
            "symbol","baseMean","log2FC","pvalue","padj","maxCooks","verdict"))
for (s in IEG) {
  id <- gid(s); if (is.na(id)) next
  r  <- res[id, ]
  mc <- max(cooks[id, ], na.rm=TRUE)
  v <- if (!is.na(r$padj)) "ok" else
       if (mc > thr) "COOK'S OUTLIER" else
       if (is.na(r$pvalue)) "all-zero / filtered" else "independent filtering"
  cat(sprintf("  %-8s %10.1f %+9.3f %10.2e %10s %10.2f  %s\n",
              s, r$baseMean, r$log2FoldChange, r$pvalue,
              ifelse(is.na(r$padj), "NA", sprintf("%.2e", r$padj)), mc, v))
}

# ================== 3. FOS PER DONOR ==================
hdr("3. FOS NORMALISED COUNTS, PER DONOR")
fid <- gid("FOS")
nc  <- counts(dds, normalized=TRUE)[fid, ]
ck  <- cooks[fid, ]
d <- data.frame(sample=names(nc), study=col$study, disease=col$disease,
                norm_count=as.numeric(nc), cooks=as.numeric(ck),
                stringsAsFactors=FALSE)
d <- d[order(-d$cooks), ]
cat(sprintf("  Cook's cut-off %.2f; samples above it force padj = NA\n\n", thr))
cat(sprintf("  %-12s %-10s %-8s %12s %9s %s\n",
            "sample","study","group","norm count","Cook's",""))
for (i in seq_len(min(12, nrow(d))))
  cat(sprintf("  %-12s %-10s %-8s %12.1f %9.2f %s\n", d$sample[i], d$study[i],
              as.character(d$disease[i]), d$norm_count[i], d$cooks[i],
              ifelse(d$cooks[i] > thr, "  <-- ABOVE CUT-OFF", "")))

cat("\n  group summaries (normalised counts):\n")
for (st in levels(col$study)) for (dg in levels(col$disease)) {
  k <- d$study==st & d$disease==dg
  if (!any(k)) next
  v <- d$norm_count[k]
  cat(sprintf("    %-10s %-8s n=%2d  median %8.1f  range %8.1f - %8.1f\n",
              st, dg, sum(k), median(v), min(v), max(v)))
}

out <- d[d$cooks > thr, ]
cat(sprintf("\n  samples above the Cook's cut-off: %d\n", nrow(out)))
if (nrow(out)) for (i in seq_len(nrow(out)))
  cat(sprintf("    %s (%s, %s) normalised count %.1f\n",
              out$sample[i], out$study[i], out$disease[i], out$norm_count[i]))

# ================== 4. THREE REFITS ==================
hdr("4. FOS UNDER THREE TREATMENTS")

r1 <- results(dds)[fid, ]
cat(sprintf("  %-34s log2FC %+7.3f  p %9.2e  padj %10s\n", "default (Cook's filtering on)",
            r1$log2FoldChange, r1$pvalue, ifelse(is.na(r1$padj),"NA",sprintf("%.2e", r1$padj))))

r2 <- results(dds, cooksCutoff=FALSE)[fid, ]
cat(sprintf("  %-34s log2FC %+7.3f  p %9.2e  padj %10s\n", "Cook's filtering OFF",
            r2$log2FoldChange, r2$pvalue, ifelse(is.na(r2$padj),"NA",sprintf("%.2e", r2$padj))))

if (nrow(out)) {
  k2   <- !(colnames(cnt) %in% out$sample)
  col2 <- col[k2, ]; c2 <- cnt[, k2]
  kp2  <- rowSums(c2 >= 10) >= min(table(col2$disease))
  dd2  <- DESeq(DESeqDataSetFromMatrix(c2[kp2,], col2, design=form), quiet=TRUE)
  r3   <- results(dd2)[fid, ]
  cat(sprintf("  %-34s log2FC %+7.3f  p %9.2e  padj %10s   (n=%d)\n",
              sprintf("outlier donor(s) removed (%d)", nrow(out)),
              r3$log2FoldChange, r3$pvalue,
              ifelse(is.na(r3$padj),"NA",sprintf("%.2e", r3$padj)), ncol(c2)))

  cat("\n  effect on the other IEGs of removing the outlier donor(s):\n")
  ra <- results(dd2)
  cat(sprintf("    %-8s %10s %10s %10s %10s\n","symbol","lfc full","lfc drop","padj full","padj drop"))
  for (s in IEG) {
    id <- gid(s); if (is.na(id) || !(id %in% rownames(ra))) next
    a <- res[id,]; b <- ra[id,]
    cat(sprintf("    %-8s %+10.3f %+10.3f %10s %10s\n", s,
                a$log2FoldChange, b$log2FoldChange,
                ifelse(is.na(a$padj),"NA",sprintf("%.2e",a$padj)),
                ifelse(is.na(b$padj),"NA",sprintf("%.2e",b$padj))))
  }
  cat(sprintf("\n  IEGs at FDR<0.05: full %d, outlier removed %d\n",
              sum(sapply(IEG, function(s){id<-gid(s); !is.na(id) && !is.na(res[id,"padj"]) && res[id,"padj"]<0.05})),
              sum(sapply(IEG, function(s){id<-gid(s); !is.na(id) && id %in% rownames(ra) && !is.na(ra[id,"padj"]) && ra[id,"padj"]<0.05}))))
}

hdr("HOW TO READ THIS")
cat("
  Cook's distance measures how much one sample moves the fitted coefficients.
  DESeq2 sets padj to NA when any sample exceeds the cut-off, as a deliberate
  refusal to report a p-value that one observation could be driving.

  padj NA but the effect unchanged with the outlier removed
      the gene is real; the NA is DESeq2 being cautious about a single extreme
      count. Report it with the outlier-removed p-value and say so.

  effect collapses when the outlier is removed
      the gene was driven by that donor and should not be reported.

  IEGs are extremely variable between individuals -- the SCZ chapter records a
  ~20x within-control range for NPAS4, and the MDD robustness run found 153x.
  A Cook's flag on an IEG is therefore expected rather than alarming, but it
  has to be checked rather than assumed.
")
cat("\n(no files written; this is a diagnostic)\n")
sessionInfo()

