
#!/usr/bin/env Rscript
# =====================================================================
# STEP 7 -- EXPLORATORY ANALYSIS, two-cohort primary
#
# The SCZ chapter presents a table of the form
#     set | variance on PC1 | leading PC1 association | diagnosis on PC1
# and that table is what makes the claim "diagnosis is not a source of
# global variance" quotable. MDD has the headline figure (PC1 = 66.7%,
# diagnosis absent) from sanity.py but not the supporting statistics.
#
# This produces, for the two-cohort primary and for each cohort alone:
#   1  variance explained by PC1-PC5
#   2  every PC tested against every covariate
#       -- Kruskal-Wallis for categorical, Pearson correlation for continuous,
#          the SCZ arm's convention
#   3  sample-distance matrices and hierarchical clustering
#   4  outlier screening at 2 and 3 SD of mean between-sample distance
#   5  covariate balance table by diagnosis
#
# Note the pooled PCA will be dominated by cohort, since the two cohorts
# differ in library prep and depth. That is expected and is itself the
# result -- the per-cohort PCAs are the informative ones for diagnosis.
#
# READ-ONLY. Writes to MDD/audit/exploratory/.
# =====================================================================

suppressPackageStartupMessages({ library(DESeq2) })

MDD <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD"
DE  <- file.path(MDD, "transcriptomics/DE")
OUT <- file.path(MDD, "audit/exploratory"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
hdr <- function(x) cat("\n", strrep("=",76), "\n", x, "\n", strrep("=",76), "\n", sep="")
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
set.seed(42)

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
  col$sex   <- factor(tolower(col$sex))
  col$study <- study
  col$depth <- colSums(cnt)
  list(cnt=cnt, col=col) }
D <- lapply(names(CO), function(n) load1(CO[[n]], n)); names(D) <- names(CO)

# pooled
g   <- Reduce(intersect, lapply(D, function(d) rownames(d$cnt)))
pc  <- do.call(cbind, lapply(D, function(d) d$cnt[g,,drop=FALSE]))
pcol<- do.call(rbind, lapply(D, function(d)
        d$col[, c("sample_id","disease","study","sex","age","pmi","ph","rin","depth")]))
rownames(pcol) <- pcol$sample_id; pcol <- pcol[colnames(pc),]
pcol$study <- factor(pcol$study)

SETS <- list(GSE101521=list(cnt=D$GSE101521$cnt, col=D$GSE101521$col),
             GSE102556=list(cnt=D$GSE102556$cnt, col=D$GSE102556$col),
             POOLED   =list(cnt=pc, col=pcol))

# ---------------------------------------------------------------- core
analyse <- function(cnt, col, label) {
  hdr(paste0(label, "  (n = ", ncol(cnt), ")"))
  keep <- rowSums(cnt >= 10) >= min(table(col$disease))
  dds  <- DESeqDataSetFromMatrix(cnt[keep,], col, design=~1)
  vsd  <- assay(vst(dds, blind=TRUE))
  cat(sprintf("  %d genes after filtering; variance-stabilised\n", sum(keep)))

  # ---- 1. PCA ----
  pca <- prcomp(t(vsd), scale.=FALSE)
  ve  <- (pca$sdev^2)/sum(pca$sdev^2)
  np  <- min(5, ncol(pca$x))
  cat(sprintf("\n  variance explained: %s\n",
      paste(sprintf("PC%d %.1f%%", 1:np, 100*ve[1:np]), collapse="  ")))

  # ---- 2. PC vs covariate association ----
  cands <- c("disease","study","sex","age","pmi","ph","rin","depth")
  cands <- cands[cands %in% names(col)]
  cat("\n  PC-covariate association (Kruskal-Wallis categorical, Pearson continuous)\n")
  cat(sprintf("  %-10s", "covariate"))
  for (i in 1:np) cat(sprintf(" %10s", paste0("PC",i)))
  cat("\n")
  res <- NULL
  for (v in cands) {
    x <- col[[v]]
    if (v %in% c("disease","study","sex")) {
      x <- factor(x); if (nlevels(x) < 2) next
      ps <- sapply(1:np, function(i) tryCatch(kruskal.test(pca$x[,i] ~ x)$p.value, error=function(e) NA))
    } else {
      x <- num(x); if (sum(!is.na(x)) < 10 || sd(x, na.rm=TRUE)==0) next
      ps <- sapply(1:np, function(i) tryCatch(cor.test(pca$x[,i], x)$p.value, error=function(e) NA))
    }
    res <- rbind(res, data.frame(set=label, covariate=v, PC=paste0("PC",1:np), p=ps))
    cat(sprintf("  %-10s", v))
    for (i in 1:np)
      cat(sprintf(" %10s", ifelse(is.na(ps[i]), "-",
        ifelse(ps[i] < 0.001, sprintf("%.1e",ps[i]), sprintf("%.3f",ps[i])))))
    cat(ifelse(any(ps < 0.05, na.rm=TRUE), "   <--\n", "\n"))
  }

  lead <- res[res$PC=="PC1" & !is.na(res$p), ]
  lead <- lead[order(lead$p), ]
  dis  <- res$p[res$covariate=="disease" & res$PC=="PC1"]
  cat(sprintf("\n  SUMMARY  variance on PC1 %.1f%%  |  leading PC1 association: %s (p = %.2e)  |  diagnosis on PC1 p = %.3f\n",
              100*ve[1], lead$covariate[1], lead$p[1], dis))
  cat(ifelse(dis > 0.05,
    "  Diagnosis is NOT associated with PC1 -- any effect is modest and not a\n  global source of variance.\n",
    "  *** diagnosis IS associated with PC1\n"))

  # ---- 3. sample distances + clustering ----
  dm <- dist(t(vsd))
  m  <- as.matrix(dm)
  mean_d <- mean(m[upper.tri(m)]); sd_d <- sd(m[upper.tri(m)])
  cat(sprintf("\n  sample-distance matrix: mean %.1f, sd %.1f\n", mean_d, sd_d))
  hc <- hclust(dm, method="ward.D2")
  k2 <- cutree(hc, k=2)
  cat("  hierarchical clustering (ward.D2), 2 clusters vs:\n")
  for (v in intersect(c("disease","study","sex"), names(col))) {
    tb <- table(k2, col[[v]])
    p  <- tryCatch(fisher.test(tb)$p.value, error=function(e) chisq.test(tb)$p.value)
    cat(sprintf("    %-8s p = %.3f%s\n", v, p,
        ifelse(p<0.05, "   <-- clusters track this", "")))
  }

  # ---- 4. outliers ----
  md <- rowMeans(m)
  o2 <- names(md)[md > mean_d + 2*sd_d]
  o3 <- names(md)[md > mean_d + 3*sd_d]
  cat(sprintf("\n  outliers: %d beyond 2 SD, %d beyond 3 SD\n", length(o2), length(o3)))
  if (length(o3)) for (s in o3)
    cat(sprintf("    3SD: %s (%s, %s)\n", s, col$disease[match(s,rownames(col))],
                if ("study" %in% names(col)) col$study[match(s,rownames(col))] else ""))
  else if (length(o2)) for (s in head(o2,6))
    cat(sprintf("    2SD: %s (%s)\n", s, col$disease[match(s,rownames(col))]))
  cat("  No samples excluded on this basis -- at these donor numbers exclusion\n")
  cat("  costs more power than it removes noise (the SCZ arm's decision).\n")

  # PCA coordinates for figures
  pcout <- data.frame(sample=rownames(pca$x), pca$x[,1:np,drop=FALSE],
                      col[rownames(pca$x), intersect(c("disease","study","sex"), names(col)),
                          drop=FALSE], mean_dist=md[rownames(pca$x)])
  write.table(pcout, file.path(OUT, sprintf("pca_%s.tsv", label)),
              sep="\t", quote=FALSE, row.names=FALSE)
  res
}

ALL <- do.call(rbind, lapply(names(SETS), function(n)
        analyse(SETS[[n]]$cnt, SETS[[n]]$col, n)))
write.table(ALL, file.path(OUT,"pc_covariate_associations.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ================= 5. COVARIATE BALANCE =================
hdr("5. COVARIATE BALANCE BY DIAGNOSIS")
cat("  Reported per cohort. Pseudoreplication distorts balance testing in both\n")
cat("  directions, so this is run only on donor-correct data.\n")
for (n in names(D)) {
  col <- D[[n]]$col
  cat(sprintf("\n  --- %s ---\n", n))
  cat(sprintf("  %-8s %14s %14s %10s\n","covariate","Control","MDD","p"))
  for (v in c("age","pmi","ph","rin","depth")) {
    if (!(v %in% names(col))) next
    x <- num(col[[v]]); if (sum(!is.na(x)) < 8) next
    a <- x[col$disease=="Control"]; b <- x[col$disease=="MDD"]
    a <- a[!is.na(a)]; b <- b[!is.na(b)]
    if (length(a)<4 || length(b)<4) next
    p <- suppressWarnings(wilcox.test(a,b)$p.value)
    fmt <- if (v=="depth") function(z) sprintf("%.1fM", z/1e6) else function(z) sprintf("%.2f", z)
    cat(sprintf("  %-8s %14s %14s %10.3f%s\n", v,
        sprintf("%s (%s)", fmt(median(a)), length(a)),
        sprintf("%s (%s)", fmt(median(b)), length(b)), p,
        ifelse(p<0.05, "   <-- DIFFERS","")))
  }
  tb <- table(col$sex, col$disease)
  cat(sprintf("  %-8s %14s %14s %10.3f\n", "sex",
      paste(names(tb[,1]), tb[,1], collapse=" ", sep="="),
      paste(names(tb[,2]), tb[,2], collapse=" ", sep="="),
      fisher.test(tb)$p.value))
}

hdr("HOW TO READ THIS")
cat("
  The pooled PCA will be dominated by study, because the cohorts differ in
  library preparation and depth. That is expected and is itself the reason
  study is a covariate in every model -- it is not a problem to fix.

  The informative panels are the per-cohort PCAs. If diagnosis is absent from
  PC1-PC5 within each cohort, the effect is modest and distributed rather than
  a global axis, which is what motivates set-level testing. It is the same
  pattern the SCZ arm reported.

  Clustering that tracks study rather than diagnosis says the same thing.
")
cat("\noutputs in", OUT, "\n")
sessionInfo()

