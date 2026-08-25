
#!/usr/bin/env Rscript
# =====================================================================
# MDD cross-cohort META-ANALYSIS
#
# Supervisor's advice, given on the ASD arm and equally applicable here:
#
#   "Instead of pooling all samples into one DESeq2 model, you analyze each
#    cohort separately, then combine the results statistically... inverse-
#    variance fixed- or random-effects meta-analysis. This is generally
#    considered a more robust approach when studies differ in protocol,
#    brain region, or sequencing platform."
#
# The MDD arm has only ever POOLED (one DESeq2 model with study as a
# covariate, plus a ComBat-seq variant). That is precisely the design she
# advised against, and the pooled result failing leave-one-dataset-out
# (retaining 0% / 7% / 33%) is the symptom she was describing.
#
# The three cohorts differ in exactly the ways that motivate meta-analysis:
#   GSE80655   DLPFC,      poly-A,           non-stranded  (-s 0)
#   GSE101521  BA9,        rRNA-depleted,    reverse       (-s 2)
#   GSE102556  dlPFC BA8/9,rRNA-depleted,    forward       (-s 1)
#
# What this produces, mirroring the SCZ arm:
#   per-cohort DE with cohort-appropriate covariates
#   inverse-variance FIXED-effect meta-analysis
#   Cochran's Q and I^2 per gene
#   DerSimonian-Laird RANDOM-effects sensitivity analysis
#   leave-one-cohort-out
#   direction concordance against the 25% three-way chance expectation
#
# GSE80655 uses the DONOR-CORRECTED matrix (47 donors, not 60 runs).
#
# READ-ONLY on inputs. Writes to MDD/audit/meta/.
# =====================================================================

suppressPackageStartupMessages({ library(DESeq2) })

MDD <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD"
DE  <- file.path(MDD, "transcriptomics/DE")
OUT <- file.path(MDD, "audit/meta"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
hdr <- function(x) cat("\n", strrep("=",76), "\n", x, "\n", strrep("=",76), "\n", sep="")
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
FDR <- 0.05; LFC <- 0.3

COHORTS <- list(
  GSE80655 = list(
    cnt = file.path(MDD,"audit/GSE80655_MDD_DLPFC_counts_donorlevel.tsv"),
    col = file.path(MDD,"audit/GSE80655_MDD_DLPFC_colData_donorlevel.tsv"),
    note = "DLPFC, 47 donors, DONOR-CORRECTED"),
  GSE101521 = list(
    cnt = file.path(MDD,"transcriptomics/GSE101521/counts/GSE101521_gene_counts_matrix.tsv"),
    col = file.path(DE,"GSE101521_colData.tsv"),
    note = "BA9, 59 donors, unmedicated by design"),
  GSE102556 = list(
    cnt = file.path(MDD,"transcriptomics/GSE102556/counts/GSE102556_gene_counts_matrix.tsv"),
    col = file.path(DE,"GSE102556_colData.tsv"),
    note = "dlPFC BA8/9, 48 donors")
)

# ===================== 1. PER-COHORT DE =====================
hdr("1. PER-COHORT DIFFERENTIAL EXPRESSION")
cat("Covariates are chosen per cohort from what is populated -- the SCZ arm's\n")
cat("approach. Forcing a common model would drop cohorts that lack a field.\n")

fits <- list()
for (nm in names(COHORTS)) {
  cf <- COHORTS[[nm]]
  if (!file.exists(cf$cnt) || !file.exists(cf$col)) {
    cat("\n  ", nm, ": MISSING input, skipped\n"); next }
  cnt <- read.delim(cf$cnt, row.names=1, check.names=FALSE)
  col <- read.delim(cf$col, check.names=FALSE); names(col)[1] <- "sample_id"
  rownames(col) <- col$sample_id
  k <- intersect(colnames(cnt), rownames(col))
  cnt <- round(as.matrix(cnt[,k,drop=FALSE])); col <- col[k,,drop=FALSE]
  col$disease <- factor(ifelse(grepl("MDD|depress", col$disease, ignore.case=TRUE),
                               "MDD","Control"), levels=c("Control","MDD"))

  use <- character(0)
  if ("sex" %in% names(col) && length(unique(col$sex)) > 1) {
    col$sex <- factor(tolower(col$sex)); use <- c(use,"sex") }
  for (v in c("age","pmi","ph","rin")) if (v %in% names(col)) {
    z <- num(col[[v]])
    if (sum(!is.na(z))==nrow(col) && sd(z)>0) {
      col[[paste0(v,"_z")]] <- scale(z)[,1]; use <- c(use, paste0(v,"_z")) } }

  keep <- rowSums(cnt >= 10) >= min(table(col$disease))
  form <- as.formula(paste("~", paste(c(use,"disease"), collapse=" + ")))
  dds  <- DESeq(DESeqDataSetFromMatrix(cnt[keep,], col, design=form), quiet=TRUE)
  r    <- as.data.frame(results(dds, contrast=c("disease","MDD","Control")))
  r$gene <- rownames(r)
  fits[[nm]] <- r

  n05 <- sum(r$padj < FDR, na.rm=TRUE)
  ne  <- sum(r$padj < FDR & abs(r$log2FoldChange) > LFC, na.rm=TRUE)
  nom <- sum(r$pvalue < 0.05, na.rm=TRUE); tst <- sum(!is.na(r$pvalue))
  cat(sprintf("\n  %-10s n=%2d (%s)\n", nm, nrow(col),
      paste(names(table(col$disease)), table(col$disease), sep="=", collapse=" ")))
  cat(sprintf("    %s\n    design %s\n", cf$note, deparse(form)))
  cat(sprintf("    tested %5d | FDR<.05 %3d | +|LFC|>0.3 %3d | p<.05 %5d | excess %.2fx\n",
              tst, n05, ne, nom, nom/(0.05*tst)))
  write.table(r[order(r$pvalue),], file.path(OUT, sprintf("de_%s.tsv", nm)),
              sep="\t", quote=FALSE, row.names=FALSE)
}
if (length(fits) < 2) stop("need at least two cohorts")

# ===================== 2. META-ANALYSIS =====================
hdr("2. INVERSE-VARIANCE META-ANALYSIS")
common <- Reduce(intersect, lapply(fits, function(f)
             f$gene[!is.na(f$log2FoldChange) & !is.na(f$lfcSE) & f$lfcSE > 0]))
cat(sprintf("  genes testable in all %d cohorts: %d\n", length(fits), length(common)))

B  <- sapply(fits, function(f) f$log2FoldChange[match(common, f$gene)])
SE <- sapply(fits, function(f) f$lfcSE[match(common, f$gene)])
W  <- 1/SE^2

beta_fe <- rowSums(B*W)/rowSums(W)
se_fe   <- sqrt(1/rowSums(W))
z_fe    <- beta_fe/se_fe
p_fe    <- 2*pnorm(-abs(z_fe))
q_fe    <- p.adjust(p_fe, "BH")

# heterogeneity
Q  <- rowSums(W*(B - beta_fe)^2)
df <- ncol(B) - 1
p_Q <- pchisq(Q, df, lower.tail=FALSE)
I2  <- pmax(0, (Q - df)/Q)*100

# DerSimonian-Laird random effects
tau2 <- pmax(0, (Q - df)/(rowSums(W) - rowSums(W^2)/rowSums(W)))
Wr   <- 1/(SE^2 + tau2)
beta_re <- rowSums(B*Wr)/rowSums(Wr)
se_re   <- sqrt(1/rowSums(Wr))
p_re    <- 2*pnorm(-abs(beta_re/se_re))
q_re    <- p.adjust(p_re, "BH")

meta <- data.frame(gene=common, beta_fe, se_fe, p_fe, q_fe,
                   beta_re, se_re, p_re, q_re, Q, p_Q, I2,
                   n_cohorts=ncol(B), stringsAsFactors=FALSE)
for (i in seq_along(fits))
  meta[[paste0("lfc_", names(fits)[i])]] <- B[,i]
meta <- meta[order(meta$p_fe), ]
write.table(meta, file.path(OUT,"meta_analysis_full.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

sig_fe <- meta$gene[meta$q_fe < FDR]
sig_re <- meta$gene[meta$q_re < FDR]
cat(sprintf("\n  FIXED effect : %d genes at FDR<%.2f\n", length(sig_fe), FDR))
cat(sprintf("  RANDOM effect: %d genes at FDR<%.2f\n", length(sig_re), FDR))
cat(sprintf("  surviving both: %d\n", length(intersect(sig_fe, sig_re))))

cat(sprintf("\n  heterogeneity across all %d genes:\n", nrow(meta)))
cat(sprintf("    median I2 = %.1f%%   genes with I2>50%%: %d (%.1f%%)\n",
            median(meta$I2, na.rm=TRUE), sum(meta$I2>50, na.rm=TRUE),
            100*mean(meta$I2>50, na.rm=TRUE)))
cat(sprintf("    Cochran Q significant (p<0.05): %d (%.1f%%)\n",
            sum(meta$p_Q<0.05, na.rm=TRUE), 100*mean(meta$p_Q<0.05, na.rm=TRUE)))
if (length(sig_fe)) {
  s <- meta[meta$gene %in% sig_fe, ]
  cat(sprintf("\n  among the %d fixed-effect hits: median I2 = %.1f%%, %d with I2>50%%\n",
              nrow(s), median(s$I2), sum(s$I2>50)))
  cat("  (in SCZ, 5 of 7 fixed-effect genes lost significance under random\n")
  cat("   effects and those lost were precisely the high-I2 genes)\n")
}

# ===================== 3. LEAVE ONE COHORT OUT =====================
hdr("3. LEAVE-ONE-COHORT-OUT")
for (drop in names(fits)) {
  idx <- setdiff(seq_along(fits), match(drop, names(fits)))
  b <- B[,idx,drop=FALSE]; w <- W[,idx,drop=FALSE]
  bf <- rowSums(b*w)/rowSums(w); sf <- sqrt(1/rowSums(w))
  q  <- p.adjust(2*pnorm(-abs(bf/sf)), "BH")
  kept <- sum(meta$gene[q < FDR] %in% sig_fe)
  cat(sprintf("  drop %-10s -> %4d genes; retains %3d of the %d full-set hits (%.0f%%)\n",
              drop, sum(q<FDR), kept, length(sig_fe),
              ifelse(length(sig_fe), 100*kept/length(sig_fe), 0)))
}

# ===================== 4. DIRECTION CONCORDANCE =====================
hdr("4. DIRECTION CONCORDANCE")
sgn <- sign(B)
allsame <- rowSums(sgn == sgn[,1]) == ncol(sgn)
chance <- 100 * 2 / 2^ncol(sgn)
cat(sprintf("  genes agreeing in direction across all %d cohorts: %d / %d (%.1f%%)\n",
            ncol(sgn), sum(allsame, na.rm=TRUE), nrow(sgn),
            100*mean(allsame, na.rm=TRUE)))
cat(sprintf("  chance expectation for %d cohorts: %.1f%%\n", ncol(sgn), chance))

cat("\n  pairwise genome-wide log2FC correlation:\n")
for (i in 1:(ncol(B)-1)) for (j in (i+1):ncol(B)) {
  ok <- is.finite(B[,i]) & is.finite(B[,j])
  cat(sprintf("    %-10s vs %-10s r = %+.3f\n",
              colnames(B)[i], colnames(B)[j], cor(B[ok,i], B[ok,j])))
}

# ===================== 5. IEG MODULE =====================
hdr("5. THE IEG MODULE IN THE META-ANALYSIS")
ieg <- c("NPAS4","FOS","FOSB","EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3",
         "ATF3","BTG2","ARC","JUN","JUNB","DUSP1","DUSP6")
sym <- NULL
if (requireNamespace("org.Hs.eg.db", quietly=TRUE)) {
  m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
        keys=ieg, keytype="SYMBOL", columns="ENSEMBL"))
  sym <- setNames(m$SYMBOL, m$ENSEMBL)
}
if (!is.null(sym)) {
  im <- meta[sub("\\..*","",meta$gene) %in% names(sym), ]
  im$symbol <- sym[sub("\\..*","",im$gene)]
  im <- im[order(im$p_fe), ]
  cat(sprintf("  %-8s %9s %9s %10s %10s %7s\n",
              "symbol","beta_FE","beta_RE","q_FE","q_RE","I2"))
  for (i in seq_len(nrow(im)))
    cat(sprintf("  %-8s %+9.3f %+9.3f %10.2e %10.2e %6.1f%%\n", im$symbol[i],
                im$beta_fe[i], im$beta_re[i], im$q_fe[i], im$q_re[i], im$I2[i]))
  cat(sprintf("\n  IEGs at FDR<0.05: fixed %d/%d, random %d/%d; median I2 %.1f%%\n",
              sum(im$q_fe<FDR), nrow(im), sum(im$q_re<FDR), nrow(im), median(im$I2)))
  write.table(im, file.path(OUT,"meta_IEG.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
}

# ===================== 6. TOP GENES =====================
hdr("6. TOP META-ANALYSIS GENES")
top <- head(meta, 30)
if (!is.null(sym)) {
  m2 <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
        keys=sub("\\..*","",top$gene), keytype="ENSEMBL", columns="SYMBOL"))
  m2 <- m2[!duplicated(m2$ENSEMBL),]
  top$symbol <- m2$SYMBOL[match(sub("\\..*","",top$gene), m2$ENSEMBL)]
} else top$symbol <- top$gene
cat(sprintf("  %-12s %9s %10s %10s %7s %6s\n","symbol","beta_FE","q_FE","q_RE","I2","dir"))
for (i in seq_len(nrow(top)))
  cat(sprintf("  %-12s %+9.3f %10.2e %10.2e %6.1f%% %6s\n",
              ifelse(is.na(top$symbol[i]), top$gene[i], top$symbol[i]),
              top$beta_fe[i], top$q_fe[i], top$q_re[i], top$I2[i],
              ifelse(top$beta_fe[i]<0,"down","up")))

hdr("HOW TO READ THIS")
cat("
  Why this is preferable to pooling
      Pooling assumes one common effect and one common dispersion across
      cohorts that differ in region, library chemistry and platform. Meta-
      analysis estimates each cohort separately and weights by precision, so
      between-study variance is modelled rather than absorbed into a covariate.

  Fixed vs random effects
      Fixed effect assumes one true effect; random effects allow it to vary
      between cohorts. Genes significant under fixed but NOT random effects
      are being carried by one cohort -- which is what leave-one-dataset-out
      showed qualitatively for the pooled analysis.

  I^2
      the proportion of variance due to between-cohort differences rather than
      sampling error. High I^2 among the significant genes means the result is
      heterogeneous and the random-effects estimate should be quoted.

  Direction concordance
      compare against the printed chance expectation (25% for three cohorts).
      Well above chance means the cohorts agree in direction even where they
      disagree on significance.
")
cat("\noutputs in", OUT, "\n")
sessionInfo()

