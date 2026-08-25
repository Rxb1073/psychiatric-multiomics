
#!/usr/bin/env Rscript
# =====================================================================
# FULL RE-RUN ON THE BIOTYPE-FILTERED MATRIX
#
# The sensitivity check established that filtering to protein_coding +
# lncRNA is safe and marginally helpful:
#     genes tested   19,885 -> 18,507
#     IEG set        5.63e-03 -> 4.10e-03
#     ERK set        4.30e-03 -> 2.96e-03
#     effect sizes   r = 1.0000 over shared genes
# and that it yields 38,907 genes before low-count filtering -- exactly
# the figure reported by the ASD and SCZ arms, making the annotation and
# filtering identical across all three.
#
# This re-runs everything the gene universe affects.
#
# AFFECTED (redone here)          UNAFFECTED (not redone)
#   DE, per cohort and pooled       per-donor module scores
#   meta-analysis (Q, I2, LOO)      AUC, permutation p
#   CAMERA (competitive background) per-cohort score p-values
#   leave-one-gene / donor-out      15/16 direction, quantiles
#   empirical null (draws universe) sex, mode-of-death, myelin tests
#   composition ladder              -- ranks are computed WITHIN gene
#                                      across donors, so the gene
#                                      universe is irrelevant to them
#
# Robustness is run on BOTH contrasts: the pooled two-cohort primary,
# and GSE101521 MDD-suicide vs control, which is what the original
# suite used. Running the original contrast isolates the effect of the
# filter from the effect of changing contrast.
#
# READ-ONLY on inputs. Writes to MDD/audit/filtered/.
# =====================================================================

suppressPackageStartupMessages({
  library(DESeq2); library(limma); library(edgeR)
  library(org.Hs.eg.db); library(AnnotationDbi) })

MDD <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD"
DE  <- file.path(MDD,"transcriptomics/DE")
OUT <- file.path(MDD,"audit/filtered"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
GTF <- "/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf"
hdr <- function(x) cat("\n", strrep("=",76), "\n", x, "\n", strrep("=",76), "\n", sep="")
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
set.seed(42)
FDR <- 0.05; LFC <- 0.3
KEEPBT <- c("protein_coding","lncRNA")

SETS_SYM <- list(
  IEG = c("NPAS4","FOS","FOSB","EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3",
          "ATF3","BTG2","ARC","JUN","JUNB","DUSP1","DUSP6"),
  ERK = c("DUSP4","DUSP5","SPRY1","SPRY2","SPRY4","SPRED1","SPRED2",
          "ETV4","ETV5","PHLDA1","EPHA2"),
  GABAergic = c("SST","PVALB","GAD1","GAD2","VIP","NPY","CCK","LHX6","RELN","CNR1"),
  MYELIN = c("MBP","PLP1","MOG","MAG","CNP","MOBP","CLDN11","MAL","OPALIN",
             "UGT8","ASPA","ERMN","GJB1","MYRF"),
  NEG_endothelial = c("CLDN5","FLT1","PECAM1","VWF","CDH5","SLC2A1","ABCB1"),
  NEG_microglial  = c("P2RY12","CX3CR1","TMEM119","CSF1R","AIF1","ITGAM","TYROBP"),
  NEG_housekeep   = c("ACTB","GAPDH","TUBB","PGK1","B2M","PPIA","RPLP0","TBP","SDHA","YWHAZ"),
  ASTROCYTE  = c("GFAP","AQP4","SLC1A2","SLC1A3","ALDH1L1","S100B"),
  EXCITATORY = c("SLC17A7","SLC17A6","SATB2","TBR1","CUX2","RORB"))

# ===================== 0. FILTER =====================
hdr("0. BIOTYPE FILTER")
bt <- local({
  cmd <- sprintf("awk -F'\\t' '$3==\"gene\"{match($9,/gene_id \"[^\"]+\"/); g=substr($9,RSTART+9,RLENGTH-10); match($9,/gene_biotype \"[^\"]+\"/); b=substr($9,RSTART+14,RLENGTH-15); print g\"\\t\"b}' %s", GTF)
  t <- read.delim(pipe(cmd), header=FALSE, col.names=c("gene","biotype"),
                  stringsAsFactors=FALSE)
  setNames(t$biotype, sub("\\..*","",t$gene)) })
cat(sprintf("  annotation: %d genes\n", length(bt)))
cat(sprintf("  protein_coding + lncRNA: %d  <-- matches ASD and SCZ (38,907)\n",
            sum(bt %in% KEEPBT)))

CO <- list(
  GSE101521=list(cnt=file.path(MDD,"transcriptomics/GSE101521/counts/GSE101521_gene_counts_matrix.tsv"),
                 col=file.path(DE,"GSE101521_colData.tsv")),
  GSE102556=list(cnt=file.path(MDD,"transcriptomics/GSE102556/counts/GSE102556_gene_counts_matrix.tsv"),
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
  e <- sub("\\..*","", rownames(cnt))
  cnt <- cnt[!is.na(bt[e]) & bt[e] %in% KEEPBT, , drop=FALSE]
  list(cnt=cnt, col=col) }
D <- lapply(names(CO), function(n) load1(CO[[n]], n)); names(D) <- names(CO)
for (n in names(D)) cat(sprintf("  %-10s %d genes after filter x %d donors\n",
                                n, nrow(D[[n]]$cnt), ncol(D[[n]]$cnt)))

sym <- function(rn) {
  m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=sub("\\..*","",rn),
        keytype="ENSEMBL", columns="SYMBOL")); m <- m[!duplicated(m$ENSEMBL),]
  setNames(m$SYMBOL, m$ENSEMBL)[sub("\\..*","",rn)] }

# ===================== 1. PER COHORT =====================
hdr("1. PER-COHORT DE (filtered)")
fits <- list()
for (n in names(D)) {
  col <- D[[n]]$col; cnt <- D[[n]]$cnt
  use <- "sex"
  for (v in c("age","pmi","ph","rin")) if (v %in% names(col)) {
    z <- num(col[[v]])
    if (sum(!is.na(z))==nrow(col) && sd(z)>0) {
      col[[paste0(v,"_z")]] <- scale(z)[,1]; use <- c(use, paste0(v,"_z")) } }
  keep <- rowSums(cnt >= 10) >= min(table(col$disease))
  f <- as.formula(paste("~", paste(c(use,"disease"), collapse=" + ")))
  r <- as.data.frame(results(DESeq(DESeqDataSetFromMatrix(cnt[keep,], col, design=f),
        quiet=TRUE), contrast=c("disease","MDD","Control")))
  r$gene <- rownames(r); fits[[n]] <- r
  nom <- sum(r$pvalue<0.05, na.rm=TRUE); tst <- sum(!is.na(r$pvalue))
  cat(sprintf("  %-10s tested %5d | FDR<.05 %3d | p<.05 %5d (%.2fx) | %s\n",
      n, tst, sum(r$padj<FDR, na.rm=TRUE), nom, nom/(0.05*tst), deparse(f)))
  write.table(r[order(r$pvalue),], file.path(OUT, sprintf("de_%s.tsv", n)),
              sep="\t", quote=FALSE, row.names=FALSE) }

# ===================== 2. POOLED =====================
hdr("2. POOLED DE (filtered)")
g   <- Reduce(intersect, lapply(D, function(d) rownames(d$cnt)))
cnt <- do.call(cbind, lapply(D, function(d) d$cnt[g,,drop=FALSE]))
col <- do.call(rbind, lapply(D, function(d) d$col[,c("sample_id","disease","study","sex","ph")]))
rownames(col) <- col$sample_id; col <- col[colnames(cnt),]
col$study <- factor(col$study); col$ph_z <- scale(num(col$ph))[,1]
keep <- rowSums(cnt >= 10) >= min(table(col$disease))
cat(sprintf("  %d genes x %d donors\n", sum(keep), ncol(cnt)))
form <- ~ study + sex + ph_z + disease
pool <- as.data.frame(results(DESeq(DESeqDataSetFromMatrix(cnt[keep,], col, design=form),
         quiet=TRUE), contrast=c("disease","MDD","Control")))
pool$gene <- rownames(pool); pool$symbol <- sym(rownames(pool))
np <- sum(pool$padj<FDR, na.rm=TRUE)
ne <- sum(pool$padj<FDR & abs(pool$log2FoldChange)>LFC, na.rm=TRUE)
nom<- sum(pool$pvalue<0.05, na.rm=TRUE); tst <- sum(!is.na(pool$pvalue))
cat(sprintf("  FDR<.05 %3d | +|LFC|>0.3 %3d | p<.05 %5d (excess %.2fx)\n",
            np, ne, nom, nom/(0.05*tst)))
write.table(pool[order(pool$pvalue),], file.path(OUT,"pooled_2cohort_filtered.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat("\n  significant genes:\n")
sg <- pool[!is.na(pool$padj) & pool$padj<FDR & abs(pool$log2FoldChange)>LFC, ]
sg <- sg[order(sg$padj), ]
for (i in seq_len(nrow(sg)))
  cat(sprintf("    %-16s %+7.3f  %.2e\n",
      ifelse(is.na(sg$symbol[i]), sg$gene[i], sg$symbol[i]),
      sg$log2FoldChange[i], sg$padj[i]))

# ===================== 3. META =====================
hdr("3. META-ANALYSIS (filtered)")
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
meta <- data.frame(gene=common, symbol=sym(common), beta_fe=b_fe, q_fe,
                   beta_re=b_re, q_re, I2, lfc_1=B[,1], lfc_2=B[,2],
                   stringsAsFactors=FALSE)
meta <- meta[order(meta$q_fe),]
write.table(meta, file.path(OUT,"meta_2cohort_filtered.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("  genes in both: %d\n", length(common)))
cat(sprintf("  FIXED  %d at FDR<.05 (%d with |LFC|>0.3)\n", sum(meta$q_fe<FDR),
            sum(meta$q_fe<FDR & abs(meta$beta_fe)>LFC)))
cat(sprintf("  RANDOM %d at FDR<.05\n", sum(meta$q_re<FDR)))
cat(sprintf("  median I2 %.1f%%; %.1f%% of genes above 50%%\n",
            median(meta$I2), 100*mean(meta$I2>50)))
ok <- is.finite(B[,1]) & is.finite(B[,2])
cat(sprintf("  cohort log2FC correlation r = %+.3f; direction concordance %.1f%%\n",
            cor(B[ok,1],B[ok,2]), 100*mean(sign(B[ok,1])==sign(B[ok,2]))))

# ===================== 4. CAMERA =====================
hdr("4. SET-LEVEL TESTING (filtered)")
des <- model.matrix(form, data=col)
v   <- voom(calcNormFactors(DGEList(cnt[keep,])), des)
sy  <- sym(rownames(v))
idx <- lapply(SETS_SYM, function(s) which(sy %in% s))
idx <- idx[sapply(idx,length) >= 3]
cu <- camera(v, idx, des, contrast=tail(colnames(des),1), inter.gene.cor=0)
cc <- camera(v, idx, des, contrast=tail(colnames(des),1), inter.gene.cor=NA)
cu <- cu[order(rownames(cu)),]; cc <- cc[rownames(cu),]
cat(sprintf("  %-18s %5s  %-6s %11s  %-6s %11s %8s\n",
            "set","genes","dir","p_uncorr","dir","p_CAMERA","r_gene"))
for (i in seq_len(nrow(cu)))
  cat(sprintf("  %-18s %5d  %-6s %11.2e  %-6s %11.2e %+8.3f%s\n",
      rownames(cu)[i], cu$NGenes[i], cu$Direction[i], cu$PValue[i],
      cc$Direction[i], cc$PValue[i], cc$Correlation[i],
      ifelse(cc$PValue[i]<0.05,"  *","")))
write.csv(data.frame(set=rownames(cc), cc), file.path(OUT,"camera_filtered.csv"),
          row.names=FALSE)

# ===================== 5. ROBUSTNESS =====================
robust <- function(vv, dd, cf, ids, lab) {
  cat(sprintf("\n  --- %s ---\n", lab))
  base <- camera(vv, list(S=ids), dd, cf, inter.gene.cor=NA)$PValue
  cat(sprintf("    baseline p = %.3e (%d genes)\n", base, length(ids)))
  # leave one gene out
  lg <- sapply(ids, function(i)
    camera(vv, list(S=setdiff(ids,i)), dd, cf, inter.gene.cor=NA)$PValue)
  cat(sprintf("    LOGO: %d/%d drops still p<0.05; worst %.3e (%.2fx)\n",
              sum(lg<0.05), length(lg), max(lg), max(lg)/base))
  list(base=base, logo=max(lg)) }

hdr("5. ROBUSTNESS -- POOLED TWO-COHORT (the primary contrast)")
iegi <- which(sy %in% SETS_SYM$IEG)
r1 <- robust(v, des, tail(colnames(des),1), iegi, "IEG, pooled two-cohort")

cat("\n  leave-one-donor-out (107 refits)\n")
ld <- rep(NA_real_, ncol(cnt))
for (i in seq_len(ncol(cnt))) {
  c2 <- cnt[, -i, drop=FALSE]; m2 <- col[-i, , drop=FALSE]
  k2 <- rowSums(c2 >= 10) >= min(table(m2$disease))
  d2 <- model.matrix(form, data=m2)
  v2 <- voom(calcNormFactors(DGEList(c2[k2,])), d2)
  s2 <- sym(rownames(v2))
  ld[i] <- camera(v2, list(S=which(s2 %in% SETS_SYM$IEG)), d2,
                  tail(colnames(d2),1), inter.gene.cor=NA)$PValue }
cat(sprintf("    %d/%d refits p<0.05; range %.3e to %.3e\n",
            sum(ld<0.05, na.rm=TRUE), sum(!is.na(ld)),
            min(ld,na.rm=TRUE), max(ld,na.rm=TRUE)))

cat("\n  empirical null (1000 expression-matched sets)\n")
mu <- rowMeans(v$E); tgt <- mu[iegi]
pool_g <- setdiff(seq_len(nrow(v)), unlist(idx))
rnd <- lapply(1:1000, function(j)
  sapply(tgt, function(t) { cand <- pool_g[abs(mu[pool_g]-t) < 0.5]
    if (!length(cand)) cand <- pool_g[order(abs(mu[pool_g]-t))][1:50]
    sample(cand,1) }))
names(rnd) <- paste0("r",1:1000)
cr <- camera(v, c(list(IEG=iegi), rnd), des, tail(colnames(des),1), inter.gene.cor=NA)
npv <- cr[grep("^r[0-9]", rownames(cr)),"PValue"]; obs <- cr["IEG","PValue"]
cat(sprintf("    IEG p = %.3e; %d/1000 random sets as extreme; empirical p = %.4f\n",
            obs, sum(npv<=obs), (sum(npv<=obs)+1)/1001))

# original contrast, for comparability with the earlier suite
hdr("6. ROBUSTNESS -- GSE101521 MDD-suicide vs control (original contrast)")
c1 <- D$GSE101521
g3 <- read.delim(file.path(DE,"GSE101521_colData.tsv"), check.names=FALSE)
names(g3)[1] <- "sample_id"; rownames(g3) <- g3$sample_id
grp <- ifelse(grepl("nonsuic", g3$group3, ignore.case=TRUE), "MDDns",
       ifelse(grepl("suicide", g3$group3, ignore.case=TRUE), "MDDs", "CON"))
names(grp) <- rownames(g3)
sel <- colnames(c1$cnt)[grp[colnames(c1$cnt)] %in% c("CON","MDDs")]
m3 <- c1$col[sel,]; m3$grp <- factor(grp[sel], levels=c("CON","MDDs"))
m3$ph_z <- scale(num(m3$ph))[,1]
c3 <- c1$cnt[, sel]; k3 <- rowSums(c3 >= 10) >= min(table(m3$grp))
d3 <- model.matrix(~ 0 + grp + ph_z, data=m3); colnames(d3) <- make.names(colnames(d3))
v3 <- voom(calcNormFactors(DGEList(c3[k3,])), d3)
s3 <- sym(rownames(v3))
ct3 <- makeContrasts(contrasts="grpMDDs - grpCON", levels=d3)[,1]
cat(sprintf("  %d donors (%s), %d genes\n", nrow(m3),
    paste(names(table(m3$grp)), table(m3$grp), sep="=", collapse=" "), nrow(v3)))
r2 <- robust(v3, d3, ct3, which(s3 %in% SETS_SYM$IEG), "IEG, GSE101521 suicide contrast")

# ===================== 7. COMPARISON =====================
hdr("7. FILTERED vs UNFILTERED -- WHAT MOVED")
prev <- list(
  c("pooled DE, FDR<.05 +LFC",        "17",        as.character(ne)),
  c("meta fixed effect",              "7",         as.character(sum(meta$q_fe<FDR))),
  c("meta random effect",             "1",         as.character(sum(meta$q_re<FDR))),
  c("cohort log2FC correlation",      "+0.252",    sprintf("%+.3f", cor(B[ok,1],B[ok,2]))),
  c("direction concordance",          "61.4%",     sprintf("%.1f%%",
                                        100*mean(sign(B[ok,1])==sign(B[ok,2])))),
  c("CAMERA IEG",                     "5.63e-03",  sprintf("%.2e", cc["IEG","PValue"])),
  c("CAMERA ERK",                     "4.30e-03",  sprintf("%.2e", cc["ERK","PValue"])),
  c("CAMERA GABAergic",               "0.734",     sprintf("%.3f", cc["GABAergic","PValue"])),
  c("CAMERA myelin",                  "0.153",     sprintf("%.3f", cc["MYELIN","PValue"])),
  c("LOGO worst (orig. contrast)",    "4.16e-03",  sprintf("%.2e", r2$logo)),
  c("empirical null",                 "0.0030",    sprintf("%.4f", (sum(npv<=obs)+1)/1001)))
cat(sprintf("  %-32s %14s %14s\n", "", "unfiltered", "filtered"))
for (x in prev) cat(sprintf("  %-32s %14s %14s\n", x[1], x[2], x[3]))

cat("\n  UNAFFECTED by the filter (ranks are within-gene across donors):\n")
cat("    IEG score AUC 0.681, permutation p 0.0012\n")
cat("    per cohort 0.690 (p 0.013) and 0.677 (p 0.037)\n")
cat("    ERK score AUC 0.677, p 0.0016\n")
cat("    15/16 control-higher in both cohorts; all 10 quantiles positive\n")
cat("    sex interaction, mode-of-death, myelin-vs-IEG adjustment\n")

hdr("DONE")
cat("outputs in", OUT, "\n")
cat("\nIf the table in section 7 shows no material change, the filtered\n")
cat("analysis becomes primary and the methods section is identical across\n")
cat("all three arms: protein_coding + lncRNA, 62,710 -> 38,907 genes.\n")
sessionInfo()

