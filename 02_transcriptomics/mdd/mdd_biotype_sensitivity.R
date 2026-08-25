
#!/usr/bin/env Rscript
# =====================================================================
# BIOTYPE FILTER -- SENSITIVITY CHECK
#
# The ASD and SCZ arms restrict to protein_coding + lncRNA before
# testing. The MDD primary analysis does not.
#
# WHY THE FILTER EXISTS IN ASD
#   An artefactual enrichment of small nucleolar RNAs, arising from the
#   contrast between poly(A)-selected and rRNA-depleted libraries.
#   snoRNAs are short, structured and non-polyadenylated, so the two
#   protocols capture them very differently.
#
# WHY IT SHOULD NOT BITE HERE
#   Both primary cohorts (GSE101521, GSE102556) are rRNA-depleted, so
#   that contrast does not exist. Already verified: 77 snoRNAs tested,
#   none significant.
#
# WHY RUN THIS ANYWAY
#   1  So the claim "it would not have mattered" is evidence, not
#      assertion.
#   2  Two of the 18 significant genes are processed pseudogenes and
#      will drop. Meanwhile the smaller multiple-testing burden
#      (19,885 -> ~18,500 rows) may admit marginal genes. The net is
#      not predictable and should be measured.
#
# WHAT IT DOES NOT DO
#   It does not re-run the arm. Changing the gene universe changes the
#   testing burden, the CAMERA competitive background, the GSEA universe
#   and the size factors, so a full re-run means redoing DE, meta,
#   CAMERA, scores, all five robustness tests, enrichment and every
#   figure. This runs the primary contrast both ways and reports the
#   difference, which is what the chapter needs.
#
# READ-ONLY. Writes to MDD/audit/biotype/.
# =====================================================================

suppressPackageStartupMessages({
  library(DESeq2); library(limma); library(edgeR)
  library(org.Hs.eg.db); library(AnnotationDbi) })

MDD <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD"
DE  <- file.path(MDD,"transcriptomics/DE")
OUT <- file.path(MDD,"audit/biotype"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
GTF <- "/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf"
hdr <- function(x) cat("\n", strrep("=",76), "\n", x, "\n", strrep("=",76), "\n", sep="")
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
set.seed(42)
FDR <- 0.05; LFC <- 0.3
KEEP_BIOTYPES <- c("protein_coding","lncRNA")

IEG <- c("NPAS4","FOS","FOSB","EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3",
         "ATF3","BTG2","ARC","JUN","JUNB","DUSP1","DUSP6")
ERK <- c("DUSP4","DUSP5","SPRY1","SPRY2","SPRY4","SPRED1","SPRED2",
         "ETV4","ETV5","PHLDA1","EPHA2")
NEG <- list(endothelial=c("CLDN5","FLT1","PECAM1","VWF","CDH5","SLC2A1","ABCB1"),
            microglial =c("P2RY12","CX3CR1","TMEM119","CSF1R","AIF1","ITGAM","TYROBP"),
            housekeep  =c("ACTB","GAPDH","TUBB","PGK1","B2M","PPIA","RPLP0","TBP","SDHA","YWHAZ"))

# ================= 0. BIOTYPES FROM THE GTF =================
hdr("0. GENE BIOTYPES FROM THE SHARED GTF")
cat("  parsing gene records (awk)\n")
bt <- local({
  cmd <- sprintf(
   "awk -F'\\t' '$3==\"gene\"{match($9,/gene_id \"[^\"]+\"/); g=substr($9,RSTART+9,RLENGTH-10); match($9,/gene_biotype \"[^\"]+\"/); b=substr($9,RSTART+14,RLENGTH-15); print g\"\\t\"b}' %s",
   GTF)
  t <- read.delim(pipe(cmd), header=FALSE, col.names=c("gene","biotype"),
                  stringsAsFactors=FALSE)
  setNames(t$biotype, sub("\\..*","",t$gene)) })
cat(sprintf("  biotypes for %d genes\n", length(bt)))
tb <- sort(table(bt), decreasing=TRUE)
cat("\n  most frequent biotypes in the annotation:\n")
for (i in seq_len(min(8,length(tb))))
  cat(sprintf("    %-28s %6d\n", names(tb)[i], tb[i]))

# ================= 1. COHORTS =================
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
  list(cnt=cnt, col=col) }
D <- lapply(names(CO), function(n) load1(CO[[n]], n)); names(D) <- names(CO)

g   <- Reduce(intersect, lapply(D, function(d) rownames(d$cnt)))
cnt <- do.call(cbind, lapply(D, function(d) d$cnt[g,,drop=FALSE]))
col <- do.call(rbind, lapply(D, function(d) d$col[,c("sample_id","disease","study","sex","ph")]))
rownames(col) <- col$sample_id; col <- col[colnames(cnt),]
col$study <- factor(col$study); col$ph_z <- scale(num(col$ph))[,1]

ens <- sub("\\..*","", rownames(cnt))
BT  <- bt[ens]
hdr("1. WHAT THE FILTER REMOVES")
cat(sprintf("  matrix: %d genes x %d donors\n", nrow(cnt), ncol(cnt)))
cat(sprintf("  biotype annotated: %d (%.1f%%)\n",
            sum(!is.na(BT)), 100*mean(!is.na(BT))))
keep_bt <- !is.na(BT) & BT %in% KEEP_BIOTYPES
cat(sprintf("  protein_coding + lncRNA: %d genes retained, %d removed\n",
            sum(keep_bt), sum(!keep_bt)))
rem <- sort(table(BT[!keep_bt]), decreasing=TRUE)
cat("\n  removed, by biotype:\n")
for (i in seq_len(min(8,length(rem))))
  cat(sprintf("    %-28s %6d\n", names(rem)[i], rem[i]))

# ================= 2. DE BOTH WAYS =================
hdr("2. DIFFERENTIAL EXPRESSION, UNFILTERED vs FILTERED")
sym <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=unique(ens),
        keytype="ENSEMBL", columns="SYMBOL"))
sym <- sym[!duplicated(sym$ENSEMBL),]
SY  <- setNames(sym$SYMBOL, sym$ENSEMBL)

run_de <- function(mask, lab) {
  cn <- cnt[mask,,drop=FALSE]
  keep <- rowSums(cn >= 10) >= min(table(col$disease))
  cn <- cn[keep,,drop=FALSE]
  dds <- DESeq(DESeqDataSetFromMatrix(cn, col, design=~ study + sex + ph_z + disease),
               quiet=TRUE)
  r <- as.data.frame(results(dds, contrast=c("disease","MDD","Control")))
  r$gene <- rownames(r); e <- sub("\\..*","", r$gene)
  r$symbol <- SY[e]; r$biotype <- bt[e]
  sig <- !is.na(r$padj) & r$padj < FDR
  sigl<- sig & abs(r$log2FoldChange) > LFC
  cat(sprintf("  %-12s tested %5d | FDR<.05 %3d | +|LFC|>0.3 %3d\n",
              lab, nrow(r), sum(sig), sum(sigl)))
  write.table(r[order(r$pvalue),], file.path(OUT, sprintf("de_%s.tsv", lab)),
              sep="\t", quote=FALSE, row.names=FALSE)
  r }

U <- run_de(rep(TRUE, nrow(cnt)), "unfiltered")
FI <- run_de(keep_bt,              "filtered")

sigset <- function(r) r$gene[!is.na(r$padj) & r$padj<FDR & abs(r$log2FoldChange)>LFC]
su <- sigset(U); sf <- sigset(FI)
cat(sprintf("\n  unfiltered %d | filtered %d | shared %d\n",
            length(su), length(sf), length(intersect(su,sf))))
lost <- setdiff(su,sf); gained <- setdiff(sf,su)
if (length(lost)) {
  cat("\n  LOST under the filter:\n")
  for (x in lost) {
    i <- match(x, U$gene)
    cat(sprintf("    %-16s %-22s log2FC %+6.3f  padj %.2e\n",
        ifelse(is.na(U$symbol[i]), x, U$symbol[i]), U$biotype[i],
        U$log2FoldChange[i], U$padj[i])) } }
if (length(gained)) {
  cat("\n  GAINED under the filter (smaller testing burden):\n")
  for (x in gained) {
    i <- match(x, FI$gene)
    cat(sprintf("    %-16s %-22s log2FC %+6.3f  padj %.2e (unfiltered padj %.3f)\n",
        ifelse(is.na(FI$symbol[i]), x, FI$symbol[i]), FI$biotype[i],
        FI$log2FoldChange[i], FI$padj[i],
        ifelse(is.na(U$padj[match(x,U$gene)]), NA, U$padj[match(x,U$gene)]))) } }

# effect-size agreement
mg <- intersect(U$gene, FI$gene)
a <- U$log2FoldChange[match(mg,U$gene)]; b <- FI$log2FoldChange[match(mg,FI$gene)]
ok <- is.finite(a)&is.finite(b)
cat(sprintf("\n  effect sizes over %d shared genes: r = %+.4f\n", sum(ok), cor(a[ok],b[ok])))

# ================= 3. snoRNA CHECK =================
hdr("3. IS THE snoRNA ARTEFACT PRESENT?")
sn <- grepl("snoRNA|snRNA|scaRNA", U$biotype)
cat(sprintf("  small nucleolar / nuclear RNAs tested: %d\n", sum(sn, na.rm=TRUE)))
if (sum(sn, na.rm=TRUE)) {
  s2 <- U[sn & !is.na(U$padj), ]
  cat(sprintf("  significant at FDR<0.05: %d\n", sum(s2$padj<FDR)))
  cat(sprintf("  median log2FC %+.3f; %d of %d negative\n",
              median(s2$log2FoldChange), sum(s2$log2FoldChange<0), nrow(s2)))
  nom <- sum(s2$pvalue < 0.05); exp_ <- 0.05*nrow(s2)
  cat(sprintf("  nominal p<0.05: %d against %.1f expected (%.2fx)\n", nom, exp_, nom/exp_))
  cat("\n  In ASD, snoRNAs were 42x enriched among significant genes and all in\n")
  cat("  one direction. If the ratio here is near 1 and none is significant, the\n")
  cat("  artefact is absent -- as expected, since both cohorts are rRNA-depleted.\n")
}

# ================= 4. SET-LEVEL BOTH WAYS =================
hdr("4. SET-LEVEL TESTING, UNFILTERED vs FILTERED")
camera_run <- function(mask, lab) {
  cn <- cnt[mask,,drop=FALSE]
  keep <- rowSums(cn >= 10) >= min(table(col$disease))
  cn <- cn[keep,,drop=FALSE]
  des <- model.matrix(~ study + sex + ph_z + disease, data=col)
  v <- voom(calcNormFactors(DGEList(cn)), des)
  sy <- SY[sub("\\..*","", rownames(v))]
  sets <- c(list(IEG=which(sy %in% IEG), ERK=which(sy %in% ERK)),
            lapply(NEG, function(x) which(sy %in% x)))
  sets <- sets[sapply(sets,length) >= 3]
  cc <- camera(v, sets, des, contrast=tail(colnames(des),1), inter.gene.cor=NA)
  cat(sprintf("\n  --- %s (%d genes tested) ---\n", lab, nrow(v)))
  for (i in seq_len(nrow(cc)))
    cat(sprintf("    %-14s %4d genes  %-6s p = %.3e%s\n", rownames(cc)[i],
        cc$NGenes[i], cc$Direction[i], cc$PValue[i],
        ifelse(cc$PValue[i]<0.05," *","")))
  cc }
cu <- camera_run(rep(TRUE,nrow(cnt)), "unfiltered")
cf <- camera_run(keep_bt,             "filtered")

hdr("VERDICT")
cat(sprintf("
  DE           unfiltered %d significant   filtered %d   shared %d
  effect sizes r = %.4f over shared genes
  IEG set      unfiltered p = %.3e        filtered p = %.3e
  ERK set      unfiltered p = %.3e        filtered p = %.3e

  If the set-level p-values are close and the shared gene set is most of
  both, the filter does not matter here and the chapter can say so with a
  number rather than an assertion:

    \"Restricting to protein-coding and lncRNA biotypes, as in the other
     arms, gives N significant genes and a set-level p of X; the finding
     is unchanged.\"

  If they diverge, re-run the arm filtered for cross-arm comparability.
", length(su), length(sf), length(intersect(su,sf)), cor(a[ok],b[ok]),
   cu["IEG","PValue"], cf["IEG","PValue"],
   ifelse("ERK" %in% rownames(cu), cu["ERK","PValue"], NA),
   ifelse("ERK" %in% rownames(cf), cf["ERK","PValue"], NA)))
cat("\noutputs in", OUT, "\n")
sessionInfo()

