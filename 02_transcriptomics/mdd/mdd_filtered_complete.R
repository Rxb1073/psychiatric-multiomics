
#!/usr/bin/env Rscript
# =====================================================================
# COMPLETE THE FILTERED RE-RUN
#
# mdd_filtered_rerun.R printed leave-one-gene-out, leave-one-donor-out
# and the empirical null but did not SAVE them, so the figure scripts
# cannot be regenerated from it. This produces every table the figure
# scripts read, under the SAME FILENAMES, in a mirror directory:
#
#     audit/filtered/mirror/primary2/pooled_2cohort.tsv
#                            /meta_2cohort.tsv
#                            /camera_pooled2.csv
#                            /IEG_2cohort.tsv
#     audit/filtered/mirror/robust/leave_one_gene_out.tsv
#                          /leave_one_donor_out.tsv
#                          /empirical_null.tsv
#     audit/filtered/mirror/suicide/suicide_decomposition.tsv
#     audit/filtered/mirror/meta/meta_analysis_full.tsv   (3-cohort, for S1)
#     audit/filtered/mirror/sexspecific/scores_by_sex.tsv
#
# The figure scripts can then be pointed at the mirror with one sed each.
#
# NOTE ON THE MODULE SCORES. I previously called them "unaffected" by the
# filter. Strictly they are not: DESeq2 size factors are computed across
# all genes, so filtering shifts them slightly and normalised counts move
# with them. Median-of-ratios is robust and 18,507 vs 19,885 genes gives
# near-identical factors, so any movement is negligible -- but the correct
# word is "structurally insensitive", not "unaffected". This script
# recomputes them on the filtered matrix and reports the difference, so
# the claim is measured rather than asserted.
#
# READ-ONLY on inputs. Writes only into audit/filtered/mirror/.
# =====================================================================

suppressPackageStartupMessages({
  library(DESeq2); library(limma); library(edgeR)
  library(org.Hs.eg.db); library(AnnotationDbi) })

MDD <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD"
DE  <- file.path(MDD,"transcriptomics/DE")
A   <- file.path(MDD,"audit")
FIL <- file.path(A,"filtered")
M   <- file.path(FIL,"mirror")
for (d in c("primary2","robust","suicide","meta","sexspecific"))
  dir.create(file.path(M,d), recursive=TRUE, showWarnings=FALSE)
GTF <- "/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf"
hdr <- function(x) cat("\n", strrep("=",76), "\n", x, "\n", strrep("=",76), "\n", sep="")
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
set.seed(42)
FDR <- 0.05; LFC <- 0.3; KEEPBT <- c("protein_coding","lncRNA")

SETS <- list(
  IEG = c("NPAS4","FOS","FOSB","EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3",
          "ATF3","BTG2","ARC","JUN","JUNB","DUSP1","DUSP6"),
  ERK = c("DUSP4","DUSP5","SPRY1","SPRY2","SPRY4","SPRED1","SPRED2",
          "ETV4","ETV5","PHLDA1","EPHA2"),
  GABAergic = c("SST","PVALB","GAD1","GAD2","VIP","NPY","CCK","LHX6","RELN","CNR1"),
  MYELIN = c("MBP","PLP1","MOG","MAG","CNP","MOBP","CLDN11","MAL","OPALIN",
             "UGT8","ASPA","ERMN","GJB1","MYRF"),
  NEG_calbindin = c("CALB1","CALB2"),
  NEG_oligodendrocyte = c("MBP","PLP1","MOG","MAG","CNP","SOX10","OLIG1",
                          "OLIG2","MOBP","CLDN11"),
  NEG_endothelial = c("CLDN5","FLT1","PECAM1","VWF","CDH5","SLC2A1","ABCB1"),
  NEG_microglial  = c("P2RY12","CX3CR1","TMEM119","CSF1R","AIF1","ITGAM","TYROBP"),
  NEG_housekeep   = c("ACTB","GAPDH","TUBB","PGK1","B2M","PPIA","RPLP0","TBP","SDHA","YWHAZ"),
  Astrocyte = c("GFAP","AQP4","SLC1A2","SLC1A3","ALDH1L1","S100B"),
  Excitatory = c("SLC17A7","SLC17A6","SATB2","TBR1","CUX2","RORB"))

hdr("0. BIOTYPE FILTER")
bt <- local({
  cmd <- sprintf("awk -F'\\t' '$3==\"gene\"{match($9,/gene_id \"[^\"]+\"/); g=substr($9,RSTART+9,RLENGTH-10); match($9,/gene_biotype \"[^\"]+\"/); b=substr($9,RSTART+14,RLENGTH-15); print g\"\\t\"b}' %s", GTF)
  t <- read.delim(pipe(cmd), header=FALSE, col.names=c("gene","biotype"),
                  stringsAsFactors=FALSE)
  setNames(t$biotype, sub("\\..*","",t$gene)) })
cat(sprintf("  retained biotypes: %d genes\n", sum(bt %in% KEEPBT)))

CO <- list(
  GSE101521=list(cnt=file.path(MDD,"transcriptomics/GSE101521/counts/GSE101521_gene_counts_matrix.tsv"),
                 col=file.path(DE,"GSE101521_colData.tsv")),
  GSE102556=list(cnt=file.path(MDD,"transcriptomics/GSE102556/counts/GSE102556_gene_counts_matrix.tsv"),
                 col=file.path(DE,"GSE102556_colData.tsv")),
  GSE80655 =list(cnt=file.path(A,"GSE80655_MDD_DLPFC_counts_donorlevel.tsv"),
                 col=file.path(A,"GSE80655_MDD_DLPFC_colData_donorlevel.tsv")))
load1 <- function(x, study, filt=TRUE) {
  if (!file.exists(x$cnt)) return(NULL)
  cnt <- read.delim(x$cnt, row.names=1, check.names=FALSE)
  col <- read.delim(x$col, check.names=FALSE); names(col)[1] <- "sample_id"
  rownames(col) <- col$sample_id
  k <- intersect(colnames(cnt), rownames(col))
  cnt <- round(as.matrix(cnt[,k,drop=FALSE])); col <- col[k,,drop=FALSE]
  col$disease <- factor(ifelse(grepl("MDD|depress", col$disease, ignore.case=TRUE),
                               "MDD","Control"), levels=c("Control","MDD"))
  if ("sex" %in% names(col)) col$sex <- factor(tolower(col$sex))
  col$study <- study
  if (filt) { e <- sub("\\..*","", rownames(cnt))
              cnt <- cnt[!is.na(bt[e]) & bt[e] %in% KEEPBT, , drop=FALSE] }
  list(cnt=cnt, col=col) }
D <- lapply(names(CO), function(n) load1(CO[[n]], n)); names(D) <- names(CO)
D <- D[!sapply(D, is.null)]
P <- D[c("GSE101521","GSE102556")]

symof <- function(rn) {
  m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=sub("\\..*","",rn),
        keytype="ENSEMBL", columns="SYMBOL")); m <- m[!duplicated(m$ENSEMBL),]
  setNames(m$SYMBOL, m$ENSEMBL)[sub("\\..*","",rn)] }

# ---------------- pooled two-cohort ----------------
hdr("1. POOLED DE + CAMERA")
g   <- Reduce(intersect, lapply(P, function(d) rownames(d$cnt)))
cnt <- do.call(cbind, lapply(P, function(d) d$cnt[g,,drop=FALSE]))
col <- do.call(rbind, lapply(P, function(d) d$col[,c("sample_id","disease","study","sex","ph")]))
rownames(col) <- col$sample_id; col <- col[colnames(cnt),]
col$study <- factor(col$study); col$ph_z <- scale(num(col$ph))[,1]
keep <- rowSums(cnt >= 10) >= min(table(col$disease))
form <- ~ study + sex + ph_z + disease
pool <- as.data.frame(results(DESeq(DESeqDataSetFromMatrix(cnt[keep,], col, design=form),
         quiet=TRUE), contrast=c("disease","MDD","Control")))
pool$gene <- rownames(pool)
write.table(pool[order(pool$pvalue),], file.path(M,"primary2","pooled_2cohort.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("  pooled: %d tested, %d at FDR<.05 +LFC\n", nrow(pool),
    sum(pool$padj<FDR & abs(pool$log2FoldChange)>LFC, na.rm=TRUE)))

des <- model.matrix(form, data=col)
v   <- voom(calcNormFactors(DGEList(cnt[keep,])), des)
sy  <- symof(rownames(v))
idx <- lapply(SETS, function(s) which(sy %in% s)); idx <- idx[sapply(idx,length)>=2]
cc  <- camera(v, idx, des, contrast=tail(colnames(des),1), inter.gene.cor=NA)
write.csv(data.frame(set=rownames(cc), cc), file.path(M,"primary2","camera_pooled2.csv"),
          row.names=FALSE)
cat(sprintf("  CAMERA IEG p = %.3e\n", cc["IEG","PValue"]))

# ---------------- per-cohort + meta ----------------
hdr("2. PER-COHORT DE AND META-ANALYSIS")
fit1 <- function(d) {
  cl <- d$col; ct <- d$cnt; use <- character(0)
  if ("sex" %in% names(cl) && nlevels(cl$sex)>1) use <- "sex"
  for (vv in c("age","pmi","ph","rin")) if (vv %in% names(cl)) {
    z <- num(cl[[vv]])
    if (sum(!is.na(z))==nrow(cl) && sd(z)>0) {
      cl[[paste0(vv,"_z")]] <- scale(z)[,1]; use <- c(use, paste0(vv,"_z")) } }
  k <- rowSums(ct >= 10) >= min(table(cl$disease))
  f <- as.formula(paste("~", paste(c(use,"disease"), collapse=" + ")))
  r <- as.data.frame(results(DESeq(DESeqDataSetFromMatrix(ct[k,], cl, design=f), quiet=TRUE),
        contrast=c("disease","MDD","Control")))
  r$gene <- rownames(r); r }

domet <- function(fits, out, lab) {
  common <- Reduce(intersect, lapply(fits, function(f)
              f$gene[!is.na(f$log2FoldChange) & !is.na(f$lfcSE) & f$lfcSE>0]))
  B <- sapply(fits, function(f) f$log2FoldChange[match(common,f$gene)])
  S <- sapply(fits, function(f) f$lfcSE[match(common,f$gene)]); W <- 1/S^2
  bf <- rowSums(B*W)/rowSums(W); sf <- sqrt(1/rowSums(W))
  qf <- p.adjust(2*pnorm(-abs(bf/sf)),"BH")
  Q <- rowSums(W*(B-bf)^2); df <- ncol(B)-1; I2 <- pmax(0,(Q-df)/Q)*100
  t2 <- pmax(0,(Q-df)/(rowSums(W)-rowSums(W^2)/rowSums(W)))
  Wr <- 1/(S^2+t2); br <- rowSums(B*Wr)/rowSums(Wr); sr <- sqrt(1/rowSums(Wr))
  qr <- p.adjust(2*pnorm(-abs(br/sr)),"BH")
  m <- data.frame(gene=common, beta_fe=bf, se_fe=sf, p_fe=2*pnorm(-abs(bf/sf)), q_fe=qf,
                  beta_re=br, se_re=sr, p_re=2*pnorm(-abs(br/sr)), q_re=qr,
                  Q=Q, p_Q=pchisq(Q,df,lower.tail=FALSE), I2=I2,
                  n_cohorts=ncol(B), stringsAsFactors=FALSE)
  for (i in seq_along(fits)) m[[paste0("lfc_", names(fits)[i])]] <- B[,i]
  # figure scripts expect lfc_101521 / lfc_102556 on the two-cohort table
  if (ncol(B)==2) { m$lfc_101521 <- B[,1]; m$lfc_102556 <- B[,2] }
  m <- m[order(m$p_fe),]
  write.table(m, out, sep="\t", quote=FALSE, row.names=FALSE)
  cat(sprintf("  %-24s %d genes | fixed %d | random %d\n", lab, nrow(m),
              sum(m$q_fe<FDR), sum(m$q_re<FDR)))
  m }

f2 <- lapply(P, fit1)
m2 <- domet(f2, file.path(M,"primary2","meta_2cohort.tsv"), "two-cohort")
if ("GSE80655" %in% names(D)) {
  f3 <- lapply(D, fit1)
  m3 <- domet(f3, file.path(M,"meta","meta_analysis_full.tsv"), "three-cohort (for S1)")
}

# IEG per-gene table for the two-cohort meta
ie <- m2[symof(m2$gene) %in% SETS$IEG, ]
ie$symbol <- symof(ie$gene)
ie$pooled_lfc  <- pool$log2FoldChange[match(ie$gene, pool$gene)]
ie$pooled_padj <- pool$padj[match(ie$gene, pool$gene)]
write.table(ie, file.path(M,"primary2","IEG_2cohort.tsv"), sep="\t",
            quote=FALSE, row.names=FALSE)

# ---------------- module scores ----------------
hdr("3. MODULE SCORES -- HOW MUCH DOES THE FILTER MOVE THEM?")
score1 <- function(d, syms, filt) {
  ct <- d$cnt
  if (!filt) {  # reload unfiltered for the comparison
    x <- CO[[unique(as.character(d$col$study))]]
    ct <- round(as.matrix(read.delim(x$cnt, row.names=1, check.names=FALSE)))
    ct <- ct[, rownames(d$col), drop=FALSE] }
  k <- rowSums(ct >= 10) >= min(table(d$col$disease))
  dd <- estimateSizeFactors(DESeqDataSetFromMatrix(ct[k,], d$col, design=~disease))
  nc <- counts(dd, normalized=TRUE); s <- symof(rownames(nc))
  ids <- which(s %in% syms); if (length(ids)<3) return(NULL)
  R <- t(apply(nc[ids,,drop=FALSE],1,function(x) rank(x)/length(x)))
  apply(R,2,median) }
auc <- function(a,b) suppressWarnings(wilcox.test(a,b)$statistic)/(length(a)*length(b))
wp  <- function(a,b) suppressWarnings(wilcox.test(a,b)$p.value)

S <- NULL
for (n in names(P)) {
  d <- P[[n]]
  sf <- score1(d, SETS$IEG, TRUE); su <- score1(d, SETS$IEG, FALSE)
  ef <- score1(d, SETS$ERK, TRUE)
  cat(sprintf("  %-10s cor(filtered, unfiltered IEG score) = %+.5f\n", n,
              cor(sf, su[names(sf)], method="spearman")))
  a <- sf[d$col$disease=="Control"]; b <- sf[d$col$disease=="MDD"]
  au <- su[d$col$disease=="Control"]; bu <- su[d$col$disease=="MDD"]
  cat(sprintf("             AUC filtered %.3f (p %.4f) | unfiltered %.3f (p %.4f)\n",
              auc(a,b), wp(a,b), auc(au,bu), wp(au,bu)))
  S <- rbind(S, data.frame(cohort=n, sample=names(sf),
        disease=as.character(d$col$disease),
        sex=if ("sex" %in% names(d$col)) as.character(d$col$sex) else NA,
        ieg=sf, erk=if (is.null(ef)) NA else ef, stringsAsFactors=FALSE)) }
write.table(S, file.path(M,"sexspecific","scores_by_sex.tsv"), sep="\t",
            quote=FALSE, row.names=FALSE)
a <- S$ieg[S$disease=="Control"]; b <- S$ieg[S$disease=="MDD"]
cat(sprintf("\n  COMBINED filtered: AUC %.3f, p %.4f  (unfiltered was 0.681, 0.0013)\n",
            auc(a,b), wp(a,b)))

# ---------------- robustness, saved ----------------
hdr("4. ROBUSTNESS TABLES (saved this time)")
iegi <- which(sy %in% SETS$IEG)
base <- camera(v, list(S=iegi), des, tail(colnames(des),1), inter.gene.cor=NA)$PValue
cat(sprintf("  baseline p = %.3e\n", base))

lg <- data.frame(dropped=sy[iegi], n=length(iegi)-1,
                 p=sapply(iegi, function(i)
                   camera(v, list(S=setdiff(iegi,i)), des, tail(colnames(des),1),
                          inter.gene.cor=NA)$PValue), stringsAsFactors=FALSE)
lg$dir <- "Down"; lg <- lg[order(-lg$p),]
write.table(lg, file.path(M,"robust","leave_one_gene_out.tsv"), sep="\t",
            quote=FALSE, row.names=FALSE)
cat(sprintf("  LOGO: %d/%d significant, worst %.3e\n", sum(lg$p<0.05), nrow(lg), max(lg$p)))

ld <- data.frame(dropped=colnames(cnt), group=as.character(col$disease),
                 p=NA_real_, dir="Down", stringsAsFactors=FALSE)
for (i in seq_len(ncol(cnt))) {
  c2 <- cnt[,-i,drop=FALSE]; m2c <- col[-i,,drop=FALSE]
  k2 <- rowSums(c2>=10) >= min(table(m2c$disease))
  d2 <- model.matrix(form, data=m2c)
  v2 <- voom(calcNormFactors(DGEList(c2[k2,])), d2); s2 <- symof(rownames(v2))
  ld$p[i] <- camera(v2, list(S=which(s2 %in% SETS$IEG)), d2, tail(colnames(d2),1),
                    inter.gene.cor=NA)$PValue }
ld <- ld[order(-ld$p),]
write.table(ld, file.path(M,"robust","leave_one_donor_out.tsv"), sep="\t",
            quote=FALSE, row.names=FALSE)
cat(sprintf("  LODO: %d/%d significant, range %.3e to %.3e\n",
            sum(ld$p<0.05), nrow(ld), min(ld$p), max(ld$p)))

mu <- rowMeans(v$E); tgt <- mu[iegi]
poolg <- setdiff(seq_len(nrow(v)), unlist(idx))
rnd <- lapply(1:1000, function(j) sapply(tgt, function(t) {
  cand <- poolg[abs(mu[poolg]-t)<0.5]
  if (!length(cand)) cand <- poolg[order(abs(mu[poolg]-t))][1:50]
  sample(cand,1) }))
names(rnd) <- paste0("rand",1:1000)
cr <- camera(v, c(list(IEG=iegi), rnd), des, tail(colnames(des),1), inter.gene.cor=NA)
write.table(data.frame(set=rownames(cr), p=cr$PValue),
            file.path(M,"robust","empirical_null.tsv"), sep="\t",
            quote=FALSE, row.names=FALSE)
np <- cr[grep("^rand",rownames(cr)),"PValue"]; ob <- cr["IEG","PValue"]
cat(sprintf("  empirical null: %d/1000 as extreme, p = %.4f\n",
            sum(np<=ob), (sum(np<=ob)+1)/1001))

# ---------------- mode of death ----------------
hdr("5. MODE OF DEATH (filtered)")
rt <- file.path(MDD,"transcriptomics/GSE102556/metadata/PRJNA398031_SraRunTable.csv")
if (file.exists(rt)) {
  r <- read.csv(rt, check.names=FALSE)
  cod <- r[[grep("cause_of_death", names(r), ignore.case=TRUE)[1]]][
           match(rownames(P$GSE102556$col), r$Run)]
  grp <- ifelse(P$GSE102556$col$disease=="MDD","MDDsui",
         ifelse(grepl("suicide",cod,ignore.case=TRUE),"CONsui","CONnat"))
  cl <- P$GSE102556$col; cl$grp <- factor(grp)
  ct <- P$GSE102556$cnt[, rownames(cl)]
  k <- rowSums(ct>=10) >= 4
  cl$ph_z <- scale(num(cl$ph))[,1]
  d4 <- model.matrix(~ 0 + grp + ph_z, data=cl); colnames(d4) <- make.names(colnames(d4))
  v4 <- voom(calcNormFactors(DGEList(ct[k,])), d4); s4 <- symof(rownames(v4))
  i4 <- lapply(SETS, function(x) which(s4 %in% x)); i4 <- i4[sapply(i4,length)>=2]
  out <- NULL
  ctr <- list(T1_modeofdeath_in_controls="grpCONsui - grpCONnat",
              T2_MDD_vs_suicide_controls="grpMDDsui - grpCONsui",
              MDD_vs_natural_controls   ="grpMDDsui - grpCONnat")
  for (nm in names(ctr)) {
    cs <- try(makeContrasts(contrasts=ctr[[nm]], levels=d4), silent=TRUE)
    if (inherits(cs,"try-error")) next
    cc4 <- camera(v4, i4, d4, contrast=cs[,1], inter.gene.cor=NA)
    fit <- eBayes(contrasts.fit(lmFit(v4,d4), cs)); tt <- topTable(fit, number=Inf, sort.by="none")
    gg <- intersect(names(s4)[s4 %in% SETS$IEG], rownames(tt))
    out <- rbind(out, data.frame(test="GSE102556", contrast=nm, n=nrow(cl),
      set=rownames(cc4), dir_corr=cc4$Direction, p_corr=cc4$PValue,
      ieg_median_logFC=median(tt[gg,"logFC"]),
      ieg_negative=sum(tt[gg,"logFC"]<0), stringsAsFactors=FALSE))
    cat(sprintf("  %-28s IEG %s p=%.3e  median logFC %+.3f\n", nm,
        cc4["IEG","Direction"], cc4["IEG","PValue"], median(tt[gg,"logFC"]))) }
  write.table(out, file.path(M,"suicide","suicide_decomposition.tsv"), sep="\t",
              quote=FALSE, row.names=FALSE) }

hdr("DONE")
cat("mirror written to", M, "\n")
cat("\nTo regenerate the figures against it:\n")
cat("  sed -i 's|audit\"|audit/filtered/mirror\"|' <figure script>\n")
cat("or point A <- file.path(MDD,\"audit/filtered/mirror\") at the top.\n")
sessionInfo()

