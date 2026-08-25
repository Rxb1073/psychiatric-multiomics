
#!/usr/bin/env Rscript
# =====================================================================
# MDD PRIMARY ANALYSIS on the two technically-matched cohorts
#
# The three-cohort meta-analysis returned 0 genes under random effects,
# and GSE80655 was the reason:
#
#     cohort       region        library          p<.05 excess   r vs others
#     GSE80655     DLPFC         poly-A (Tn-seq)      0.34x      +0.06 / -0.03
#     GSE101521    BA9           rRNA-depleted        1.56x      -
#     GSE102556    dlPFC BA8/9   rRNA-depleted        1.00x      +0.23 with 101521
#
# An excess BELOW chance means GSE80655 is noisier than its n suggests, and
# it is essentially uncorrelated with both other cohorts. It is also the
# cohort that carried the pseudoreplication and comes from a four-diagnosis
# deposit.
#
# BA9 is part of dlPFC, so GSE101521 and GSE102556 share BOTH region and
# library chemistry -- which is presumably why they are the two that agree.
#
# This mirrors the supervisor's advice on the ASD arm: drop the problematic
# cohort, and estimate effects within each valid cohort before combining.
#
# Runs BOTH combination strategies, because they answer slightly different
# questions and should agree if the effect is real:
#
#   META-ANALYSIS  summary statistics, models between-cohort variance,
#                  conservative, makes no assumption of a shared dispersion
#   POOLING        individual-level data with study as a covariate, more
#                  powerful when the effect genuinely is shared
#
# GSE80655 is retained only as a sensitivity check, never as primary.
#
# READ-ONLY. Writes to MDD/audit/primary2/.
# =====================================================================

suppressPackageStartupMessages({ library(DESeq2); library(limma); library(edgeR) })

MDD <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD"
DE  <- file.path(MDD, "transcriptomics/DE")
OUT <- file.path(MDD, "audit/primary2"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
hdr <- function(x) cat("\n", strrep("=",76), "\n", x, "\n", strrep("=",76), "\n", sep="")
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
FDR <- 0.05; LFC <- 0.3

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
  list(cnt=cnt, col=col)
}
D <- lapply(names(CO), function(n) load1(CO[[n]], n)); names(D) <- names(CO)

hdr("0. COHORTS")
for (n in names(D)) {
  cat(sprintf("  %-10s %d genes x %d donors  (%s)\n", n, nrow(D[[n]]$cnt), ncol(D[[n]]$cnt),
      paste(names(table(D[[n]]$col$disease)), table(D[[n]]$col$disease), sep="=", collapse=" ")))
}
tot <- sum(sapply(D, function(d) ncol(d$cnt)))
alld <- unlist(lapply(D, function(d) as.character(d$col$disease)))
cat(sprintf("  TOTAL %d donors (%s)\n", tot,
    paste(names(table(alld)), table(alld), sep="=", collapse=" ")))

# which covariates are populated in BOTH
usable <- function(col, v) {
  if (!(v %in% names(col))) return(FALSE)
  z <- num(col[[v]]); sum(!is.na(z)) == nrow(col) && sd(z) > 0 }
COMMON <- Filter(function(v) all(sapply(D, function(d) usable(d$col, v))),
                 c("age","pmi","ph","rin"))
cat("\n  numeric covariates populated in both:",
    if (length(COMMON)) paste(COMMON, collapse=", ") else "(none)", "\n")
cat("  sex is available in both and is always included.\n")

# ================== 1. PER-COHORT ==================
hdr("1. PER-COHORT DE (each with its own covariates)")
fits <- list()
for (n in names(D)) {
  col <- D[[n]]$col; cnt <- D[[n]]$cnt
  col$sex <- factor(tolower(col$sex)); use <- "sex"
  for (v in c("age","pmi","ph","rin")) if (usable(col, v)) {
    col[[paste0(v,"_z")]] <- scale(num(col[[v]]))[,1]; use <- c(use, paste0(v,"_z")) }
  keep <- rowSums(cnt >= 10) >= min(table(col$disease))
  f <- as.formula(paste("~", paste(c(use,"disease"), collapse=" + ")))
  r <- as.data.frame(results(DESeq(DESeqDataSetFromMatrix(cnt[keep,], col, design=f), quiet=TRUE),
                             contrast=c("disease","MDD","Control")))
  r$gene <- rownames(r); fits[[n]] <- r
  cat(sprintf("  %-10s design %-46s FDR<.05 %3d   p<.05 %5d (%.2fx)\n", n, deparse(f),
      sum(r$padj<FDR, na.rm=TRUE), sum(r$pvalue<0.05, na.rm=TRUE),
      sum(r$pvalue<0.05, na.rm=TRUE)/(0.05*sum(!is.na(r$pvalue)))))
  write.table(r[order(r$pvalue),], file.path(OUT, sprintf("de_%s.tsv", n)),
              sep="\t", quote=FALSE, row.names=FALSE)
}

# ================== 2. META-ANALYSIS ==================
hdr("2. META-ANALYSIS of the two cohorts")
common <- Reduce(intersect, lapply(fits, function(f)
             f$gene[!is.na(f$log2FoldChange) & !is.na(f$lfcSE) & f$lfcSE>0]))
B  <- sapply(fits, function(f) f$log2FoldChange[match(common,f$gene)])
SE <- sapply(fits, function(f) f$lfcSE[match(common,f$gene)])
W  <- 1/SE^2
b_fe <- rowSums(B*W)/rowSums(W); s_fe <- sqrt(1/rowSums(W))
q_fe <- p.adjust(2*pnorm(-abs(b_fe/s_fe)), "BH")
Q <- rowSums(W*(B-b_fe)^2); df <- ncol(B)-1
I2 <- pmax(0,(Q-df)/Q)*100
tau2 <- pmax(0,(Q-df)/(rowSums(W)-rowSums(W^2)/rowSums(W)))
Wr <- 1/(SE^2+tau2); b_re <- rowSums(B*Wr)/rowSums(Wr); s_re <- sqrt(1/rowSums(Wr))
q_re <- p.adjust(2*pnorm(-abs(b_re/s_re)), "BH")

meta <- data.frame(gene=common, beta_fe=b_fe, q_fe, beta_re=b_re, q_re, I2,
                   lfc_101521=B[,1], lfc_102556=B[,2], stringsAsFactors=FALSE)
meta <- meta[order(meta$q_fe), ]
write.table(meta, file.path(OUT,"meta_2cohort.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("  genes in both: %d\n", length(common)))
cat(sprintf("  FIXED  effect: %d at FDR<%.2f\n", sum(meta$q_fe<FDR), FDR))
cat(sprintf("  RANDOM effect: %d at FDR<%.2f\n", sum(meta$q_re<FDR), FDR))
cat(sprintf("  median I2 = %.1f%%   I2>50%%: %.1f%% of genes\n",
            median(meta$I2), 100*mean(meta$I2>50)))
ok <- is.finite(B[,1]) & is.finite(B[,2])
cat(sprintf("  genome-wide log2FC correlation between cohorts: r = %+.3f\n", cor(B[ok,1],B[ok,2])))
cat(sprintf("  direction concordance: %.1f%% (chance 50%%)\n",
            100*mean(sign(B[ok,1])==sign(B[ok,2]))))

# ================== 3. POOLED ==================
hdr("3. POOLED DE on individual-level data")
g <- Reduce(intersect, lapply(D, function(d) rownames(d$cnt)))
cnt <- do.call(cbind, lapply(D, function(d) d$cnt[g,,drop=FALSE]))
col <- do.call(rbind, lapply(D, function(d) {
  c0 <- d$col[, intersect(c("sample_id","disease","study","sex","age","pmi","ph","rin"),
                          names(d$col)), drop=FALSE]
  for (v in c("age","pmi","ph","rin")) if (!(v %in% names(c0))) c0[[v]] <- NA
  c0[, c("sample_id","disease","study","sex","age","pmi","ph","rin")] }))
rownames(col) <- col$sample_id; col <- col[colnames(cnt), ]
col$study <- factor(col$study); col$sex <- factor(tolower(col$sex))
use <- "sex"
for (v in COMMON) { col[[paste0(v,"_z")]] <- scale(num(col[[v]]))[,1]; use <- c(use, paste0(v,"_z")) }
keep <- rowSums(cnt >= 10) >= min(table(col$disease))
f <- as.formula(paste("~ study +", paste(c(use,"disease"), collapse=" + ")))
cat(sprintf("  %d genes x %d donors;  design %s\n", sum(keep), ncol(cnt), deparse(f)))
pooled <- as.data.frame(results(DESeq(DESeqDataSetFromMatrix(cnt[keep,], col, design=f),
                        quiet=TRUE), contrast=c("disease","MDD","Control")))
pooled$gene <- rownames(pooled)
np <- sum(pooled$padj<FDR, na.rm=TRUE)
ne <- sum(pooled$padj<FDR & abs(pooled$log2FoldChange)>LFC, na.rm=TRUE)
nom <- sum(pooled$pvalue<0.05, na.rm=TRUE); tst <- sum(!is.na(pooled$pvalue))
cat(sprintf("  FDR<.05 %3d  |  +|LFC|>0.3 %3d  |  p<.05 %5d  |  excess %.2fx\n",
            np, ne, nom, nom/(0.05*tst)))
write.table(pooled[order(pooled$pvalue),], file.path(OUT,"pooled_2cohort.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# ================== 4. DO THEY AGREE? ==================
hdr("4. META vs POOLED")
mp <- merge(meta[,c("gene","beta_fe","q_fe","q_re")],
            pooled[,c("gene","log2FoldChange","padj")], by="gene")
ok <- is.finite(mp$beta_fe) & is.finite(mp$log2FoldChange)
cat(sprintf("  effect-size correlation: r = %+.3f over %d genes\n",
            cor(mp$beta_fe[ok], mp$log2FoldChange[ok]), sum(ok)))
sm <- mp$gene[mp$q_fe<FDR]; sp <- mp$gene[!is.na(mp$padj) & mp$padj<FDR]
cat(sprintf("  meta (fixed) %d | pooled %d | overlap %d\n",
            length(sm), length(sp), length(intersect(sm,sp))))
if (length(sp) && length(sm))
  cat(sprintf("  %.0f%% of pooled hits are also meta hits\n",
              100*length(intersect(sm,sp))/length(sp)))
cat("\n  Agreement between the two supports the result. Pooled-only hits are\n")
cat("  candidates for being driven by one cohort -- check their I2 in the meta table.\n")

# ================== 5. SET-LEVEL ==================
hdr("5. SET-LEVEL TESTING on the pooled two-cohort data")
SETS_SYM <- list(
  IEG = c("NPAS4","FOS","FOSB","EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3",
          "ATF3","BTG2","ARC","JUN","JUNB","DUSP1","DUSP6"),
  GABAergic = c("SST","PVALB","GAD1","GAD2","VIP","NPY","CCK","LHX6","RELN","CNR1"),
  NEGCTRL_calbindin = c("CALB1","CALB2"),
  NEGCTRL_oligodendrocyte = c("MBP","PLP1","MOG","MAG","CNP","SOX10","OLIG1","OLIG2","MOBP","CLDN11"),
  Astrocyte = c("GFAP","AQP4","SLC1A2","SLC1A3","ALDH1L1","S100B"))
SETS <- lapply(SETS_SYM, function(s) {
  m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
        keys=s, keytype="SYMBOL", columns="ENSEMBL")); unique(na.omit(m$ENSEMBL)) })
des <- model.matrix(f, data=col)
d <- DGEList(cnt[keep,]); d <- calcNormFactors(d); v <- voom(d, des)
idx <- lapply(SETS, function(x) which(sub("\\..*","",rownames(v)) %in% x))
idx <- idx[sapply(idx,length) >= 2]
cu <- camera(v, idx, des, contrast=tail(colnames(des),1), inter.gene.cor=0)
cc <- camera(v, idx, des, contrast=tail(colnames(des),1), inter.gene.cor=NA)
cu <- cu[order(rownames(cu)),]; cc <- cc[rownames(cu),]
cat(sprintf("  %-24s %5s  %-6s %10s  %-6s %10s %8s\n",
            "set","genes","dir","p_uncorr","dir","p_CAMERA","r_gene"))
for (i in seq_len(nrow(cu)))
  cat(sprintf("  %-24s %5d  %-6s %10.2e  %-6s %10.2e %+8.3f%s\n",
      rownames(cu)[i], cu$NGenes[i], cu$Direction[i], cu$PValue[i],
      cc$Direction[i], cc$PValue[i], cc$Correlation[i],
      ifelse(grepl("NEGCTRL", rownames(cu)[i]) & cc$PValue[i] < 0.05, "  *FIRED*","")))
write.csv(data.frame(set=rownames(cc), cc), file.path(OUT,"camera_pooled2.csv"), row.names=FALSE)

# ================== 6. IEG TABLE ==================
hdr("6. IEG MODULE: per cohort, meta, pooled")
m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
      keys=SETS_SYM$IEG, keytype="SYMBOL", columns="ENSEMBL"))
sym <- setNames(m$SYMBOL, m$ENSEMBL)
ie <- meta[sub("\\..*","",meta$gene) %in% names(sym), ]
ie$symbol <- sym[sub("\\..*","",ie$gene)]
ie$pooled_lfc  <- pooled$log2FoldChange[match(ie$gene, pooled$gene)]
ie$pooled_padj <- pooled$padj[match(ie$gene, pooled$gene)]
ie <- ie[order(ie$q_fe), ]
cat(sprintf("  %-8s %8s %8s %9s %9s %7s %9s %10s\n",
            "symbol","lfc_101","lfc_102","beta_FE","q_FE","I2","pool_lfc","pool_padj"))
for (i in seq_len(nrow(ie)))
  cat(sprintf("  %-8s %+8.3f %+8.3f %+9.3f %9.2e %6.1f%% %+9.3f %10.2e\n",
      ie$symbol[i], ie$lfc_101521[i], ie$lfc_102556[i], ie$beta_fe[i],
      ie$q_fe[i], ie$I2[i], ie$pooled_lfc[i], ie$pooled_padj[i]))
cat(sprintf("\n  negative in both cohorts: %d / %d\n",
            sum(ie$lfc_101521<0 & ie$lfc_102556<0, na.rm=TRUE), nrow(ie)))
cat(sprintf("  FDR<0.05: meta-fixed %d, meta-random %d, pooled %d\n",
            sum(ie$q_fe<FDR), sum(ie$q_re<FDR), sum(ie$pooled_padj<FDR, na.rm=TRUE)))
write.table(ie, file.path(OUT,"IEG_2cohort.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

hdr("HOW TO READ THIS")
cat("
  Why two cohorts rather than three
      GSE80655 correlates with neither other cohort (r = +0.06, -0.03) and has
      a nominal p<0.05 excess BELOW chance (0.34x), meaning it is noisier than
      its sample size implies. It also differs in library chemistry: poly-A
      selection against rRNA depletion in the other two. Excluding it is a
      documented design decision, not a fishing expedition -- and it mirrors
      the advice given on the ASD arm.

  Meta versus pooled
      They should agree. Meta-analysis is conservative and models between-
      cohort variance; pooling is more powerful but assumes a shared effect
      and dispersion. Genes significant when pooled but NOT in the meta-
      analysis, with high I2, are being carried by one cohort.

  Report whichever is more conservative as primary, and say that both were run.
")
cat("\noutputs in", OUT, "\n")
sessionInfo()
