
#!/usr/bin/env Rscript
# =====================================================================
# Two questions.
#
# 1  IS THE IEG SCORE CONFOUNDED BY SEQUENCING DEPTH?
#    The previous run found rho(IEG score, depth) = -0.531. But it built
#    the score by ranking across BOTH cohorts pooled, and the cohorts
#    differ two-fold in depth (12M vs 24M). That alone would induce the
#    correlation. The main analysis ranked WITHIN cohort. This recomputes
#    the correlation within each cohort separately, where cohort-level
#    depth differences cannot contribute.
#
# 2  IS ERK ACTIVITY UP IN FEMALE MDD, AS LABONTE REPORTED?
#    Labonte et al. 2017 (Nat Med), the source study for GSE102556, found
#    DUSP6 downregulated in female MDD prefrontal cortex and inferred
#    INCREASED ERK signalling -- confirmed as raised pERK protein in
#    females of both species. DUSP6 is a phosphatase that inactivates ERK,
#    so less DUSP6 means more active ERK.
#
#    pERK cannot be measured from RNA-seq. But ERK activity induces a set
#    of negative-feedback genes -- DUSP4/5/6, SPRY1/2/4, ETV4/5, PHLDA1 --
#    which are transcriptional readouts of ERK OUTPUT. If ERK output were
#    high, these should be UP.
#
#    Our IEG result predicts the opposite: FOS and EGR1 are canonical
#    ERK target genes and both are DOWN. So the transcriptional data and
#    the phospho-protein data may disagree, and it is worth showing
#    exactly how rather than asserting either.
#
#    CAVEAT stated up front: inferring kinase activity from transcript
#    levels is indirect. A negative-feedback gene can be low because ERK
#    output is low, OR because the feedback loop itself is disrupted.
#    This test constrains the interpretation; it does not settle it.
#
# 3  CIRCULARITY CHECK. FOS, JUN, DUSP1 and DUSP6 are annotated members
#    of KEGG MAPK signalling. The earlier "MAPK pathway down" result may
#    just be re-detecting the IEGs. Re-tested with the IEGs removed.
#
# READ-ONLY. Writes to MDD/audit/erk/.
# =====================================================================

suppressPackageStartupMessages({
  library(DESeq2); library(limma); library(edgeR)
  library(org.Hs.eg.db); library(AnnotationDbi) })

MDD <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD"
DE  <- file.path(MDD, "transcriptomics/DE")
OUT <- file.path(MDD, "audit/erk"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
hdr <- function(x) cat("\n", strrep("=",76), "\n", x, "\n", strrep("=",76), "\n", sep="")
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
set.seed(42)

IEG <- c("NPAS4","FOS","FOSB","EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3",
         "ATF3","BTG2","ARC","JUN","JUNB","DUSP1","DUSP6")
# ERK negative-feedback genes: induced BY ERK activity, so a readout of ERK OUTPUT
ERK_OUT <- c("DUSP4","DUSP5","DUSP6","SPRY1","SPRY2","SPRY4","SPRED1","SPRED2",
             "ETV4","ETV5","PHLDA1","EPHA2","CCND1","MYC")
ERK_OUT_NOIEG <- setdiff(ERK_OUT, IEG)     # drops DUSP6, the only overlap

CO <- list(
  GSE101521 = list(cnt=file.path(MDD,"transcriptomics/GSE101521/counts/GSE101521_gene_counts_matrix.tsv"),
                   col=file.path(DE,"GSE101521_colData.tsv")),
  GSE102556 = list(cnt=file.path(MDD,"transcriptomics/GSE102556/counts/GSE102556_gene_counts_matrix.tsv"),
                   col=file.path(DE,"GSE102556_colData.tsv")))

load1 <- function(x) {
  cnt <- read.delim(x$cnt, row.names=1, check.names=FALSE)
  col <- read.delim(x$col, check.names=FALSE); names(col)[1] <- "sample_id"
  rownames(col) <- col$sample_id
  k <- intersect(colnames(cnt), rownames(col))
  cnt <- round(as.matrix(cnt[,k,drop=FALSE])); col <- col[k,,drop=FALSE]
  col$disease <- factor(ifelse(grepl("MDD|depress", col$disease, ignore.case=TRUE),
                               "MDD","Control"), levels=c("Control","MDD"))
  col$sex <- factor(tolower(col$sex))
  col$depth <- colSums(cnt)
  keep <- rowSums(cnt >= 10) >= min(table(col$disease))
  dds <- estimateSizeFactors(DESeqDataSetFromMatrix(cnt[keep,], col, design=~disease))
  list(cnt=cnt[keep,], nc=counts(dds, normalized=TRUE), col=col) }
D <- lapply(CO, load1)

SYM <- function(rn) {
  m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=sub("\\..*","",rn),
        keytype="ENSEMBL", columns="SYMBOL"))
  m <- m[!duplicated(m$ENSEMBL),]
  setNames(m$SYMBOL, m$ENSEMBL)[sub("\\..*","",rn)] }

score <- function(nc, syms, sy) {
  ids <- rownames(nc)[which(sy %in% syms)]
  if (length(ids) < 3) return(NULL)
  R <- t(apply(nc[ids,,drop=FALSE], 1, function(v) rank(v)/length(v)))
  apply(R, 2, median) }
auc <- function(a,b) suppressWarnings(wilcox.test(a,b)$statistic)/(length(a)*length(b))
wp  <- function(a,b) suppressWarnings(wilcox.test(a,b)$p.value)

# ================= 1. DEPTH, WITHIN COHORT =================
hdr("1. IS THE IEG SCORE CONFOUNDED BY DEPTH? -- WITHIN-COHORT TEST")
cat("  Pooled across cohorts the correlation was -0.531, but the cohorts\n")
cat("  differ two-fold in depth, which alone would produce that. Within\n")
cat("  cohort, depth cannot carry a cohort difference.\n")
for (n in names(D)) {
  sy <- SYM(rownames(D[[n]]$nc))
  s  <- score(D[[n]]$nc, IEG, sy)
  cl <- D[[n]]$col
  cat(sprintf("\n  --- %s (n=%d, depth %.1f-%.1fM) ---\n", n, length(s),
              min(cl$depth)/1e6, max(cl$depth)/1e6))
  for (v in c("depth","rin","ph")) {
    x <- if (v=="depth") log10(cl$depth) else num(cl[[v]])
    ok <- is.finite(x) & is.finite(s)
    if (sum(ok) < 10) next
    cat(sprintf("    rho(IEG score, %-5s) = %+.3f  (n=%d)%s\n", v,
                cor(s[ok], x[ok], method="spearman"), sum(ok),
                ifelse(abs(cor(s[ok],x[ok],method="spearman"))>0.4,"   <-- strong","")))
  }
  a <- s[cl$disease=="Control"]; b <- s[cl$disease=="MDD"]
  cat(sprintf("    IEG score by diagnosis: AUC %.3f  p %.4f\n", auc(a,b), wp(a,b)))
}
cat("\n  If the within-cohort depth correlations are small, the -0.531 was an\n")
cat("  artefact of pooling and the main analysis is unaffected.\n")

# ================= 2. SEX-STRATIFIED, PER COHORT =================
hdr("2. SEX-STRATIFIED DE, PER COHORT")
cat("  Labonte reported DUSP6 DOWN in female MDD prefrontal cortex.\n")
cat("  GSE102556 is that study's data, so this should reproduce.\n")

de_sub <- function(d, mask, lab) {
  cl <- d$col[mask,,drop=FALSE]; ct <- d$cnt[,mask,drop=FALSE]
  if (length(unique(cl$disease)) < 2 || min(table(cl$disease)) < 4) {
    cat(sprintf("    %-28s too few (%s)\n", lab,
        paste(table(cl$disease), collapse="/"))); return(NULL) }
  k <- rowSums(ct >= 10) >= min(table(cl$disease))
  dd <- DESeq(DESeqDataSetFromMatrix(ct[k,], droplevels(cl), design=~disease), quiet=TRUE)
  r <- as.data.frame(results(dd, contrast=c("disease","MDD","Control")))
  r$symbol <- SYM(rownames(r)); r }

RES <- list()
for (n in names(D)) {
  cl <- D[[n]]$col
  cat(sprintf("\n  --- %s ---\n", n))
  print(table(cl$sex, cl$disease, dnn=c("sex","disease")))
  for (s in c("all", levels(cl$sex))) {
    mask <- if (s=="all") rep(TRUE, nrow(cl)) else cl$sex==s
    r <- de_sub(D[[n]], mask, s); if (is.null(r)) next
    RES[[paste(n,s)]] <- r
    d6 <- r[which(r$symbol=="DUSP6")[1], ]
    cat(sprintf("    %-8s n=%2d   DUSP6 log2FC %+6.3f  p %8.4f  padj %8s\n",
        s, sum(mask), d6$log2FoldChange, d6$pvalue,
        ifelse(is.na(d6$padj),"NA",sprintf("%.3f", d6$padj))))
  }
}

# ================= 3. ERK OUTPUT SIGNATURE =================
hdr("3. ERK OUTPUT SIGNATURE -- DOES THE TRANSCRIPTOME SAY ERK IS UP?")
cat("  These genes are INDUCED BY ERK activity (negative feedback):\n")
cat("   ", paste(ERK_OUT, collapse=", "), "\n")
cat("  ERK output HIGH  -> these UP.  ERK output LOW -> these DOWN.\n")
cat("  Labonte's mouse work implies ERK UP in female MDD.\n\n")

for (nm in names(RES)) {
  r <- RES[[nm]]
  e <- r[r$symbol %in% ERK_OUT & !is.na(r$log2FoldChange), ]
  e2 <- r[r$symbol %in% ERK_OUT_NOIEG & !is.na(r$log2FoldChange), ]
  i <- r[r$symbol %in% IEG & !is.na(r$log2FoldChange), ]
  if (!nrow(e)) next
  cat(sprintf("  %-22s ERK-out median %+6.3f (%d/%d down) | excl. IEG %+6.3f | IEG %+6.3f\n",
      nm, median(e$log2FoldChange), sum(e$log2FoldChange<0), nrow(e),
      ifelse(nrow(e2), median(e2$log2FoldChange), NA),
      ifelse(nrow(i),  median(i$log2FoldChange),  NA)))
}

cat("\n  per-gene, both cohorts, all donors:\n")
cat(sprintf("  %-9s %22s %22s\n", "gene", "GSE101521 lfc (p)", "GSE102556 lfc (p)"))
for (g in ERK_OUT) {
  out <- sprintf("  %-9s", g)
  for (n in names(D)) {
    r <- RES[[paste(n,"all")]]
    if (is.null(r)) { out <- paste0(out, sprintf(" %22s","-")); next }
    j <- which(r$symbol==g)[1]
    out <- paste0(out, if (is.na(j)) sprintf(" %22s","not tested")
                  else sprintf(" %+13.3f (%.3f)", r$log2FoldChange[j], r$pvalue[j]))
  }
  cat(out, "\n")
}

# score-based, per cohort per sex
cat("\n  ERK-output SCORE (excluding IEG overlap), by cohort and sex:\n")
for (n in names(D)) {
  sy <- SYM(rownames(D[[n]]$nc)); cl <- D[[n]]$col
  s <- score(D[[n]]$nc, ERK_OUT_NOIEG, sy)
  if (is.null(s)) { cat(sprintf("    %-10s too few genes\n", n)); next }
  for (sx in c("all", levels(cl$sex))) {
    m <- if (sx=="all") rep(TRUE,nrow(cl)) else cl$sex==sx
    a <- s[m & cl$disease=="Control"]; b <- s[m & cl$disease=="MDD"]
    if (length(a)<4 || length(b)<4) next
    cat(sprintf("    %-10s %-7s Control %.3f  MDD %.3f   AUC %.3f  p %.4f%s\n",
        n, sx, median(a), median(b), auc(a,b), wp(a,b),
        ifelse(wp(a,b)<0.05,"  *","")))
  }
}

# ================= 4. MAPK CIRCULARITY =================
hdr("4. IS 'MAPK PATHWAY DOWN' INDEPENDENT OF THE IEGs?")
cat("  FOS, JUN, DUSP1 and DUSP6 are KEGG MAPK members, so the earlier\n")
cat("  result may simply be re-detecting the IEG set. Re-run CAMERA on the\n")
cat("  pooled two-cohort data with MAPK tested twice: intact, and with the\n")
cat("  16 IEGs removed.\n")

g   <- Reduce(intersect, lapply(D, function(d) rownames(d$cnt)))
cnt <- do.call(cbind, lapply(D, function(d) d$cnt[g,,drop=FALSE]))
col <- do.call(rbind, lapply(names(D), function(n) {
  x <- D[[n]]$col[, c("sample_id","disease","sex","ph")]; x$study <- n; x }))
rownames(col) <- col$sample_id; col <- col[colnames(cnt),]
col$study <- factor(col$study); col$ph_z <- scale(num(col$ph))[,1]
des <- model.matrix(~ study + sex + ph_z + disease, data=col)
v <- voom(calcNormFactors(DGEList(cnt)), des)
sy <- SYM(rownames(v))

mapk <- try(suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
        keys="04010", keytype="PATH", columns="SYMBOL")), silent=TRUE)
MAPK <- if (inherits(mapk,"try-error")) character(0) else unique(mapk$SYMBOL)
if (!length(MAPK)) {
  cat("  KEGG PATH lookup unavailable; using a curated core MAPK/ERK list\n")
  MAPK <- c("MAP2K1","MAP2K2","MAPK1","MAPK3","RAF1","BRAF","ARAF","KRAS","NRAS",
            "HRAS","SOS1","SOS2","GRB2","SHC1","EGFR","FGFR1","NTRK2","RPS6KA1",
            "RPS6KA2","RPS6KA3","ELK1","ELK4","MKNK1","MKNK2","DUSP1","DUSP4",
            "DUSP5","DUSP6","FOS","JUN","MYC") }
cat(sprintf("  MAPK set: %d symbols; %d overlap the IEG set\n",
            length(MAPK), length(intersect(MAPK, IEG))))
sets <- list(MAPK_intact = which(sy %in% MAPK),
             MAPK_no_IEG = which(sy %in% setdiff(MAPK, IEG)),
             IEG_only    = which(sy %in% IEG))
sets <- sets[sapply(sets, length) >= 3]
cc <- camera(v, sets, des, contrast=tail(colnames(des),1), inter.gene.cor=NA)
cat(sprintf("\n  %-14s %6s %-6s %12s\n","set","genes","dir","p_CAMERA"))
for (i in seq_len(nrow(cc)))
  cat(sprintf("  %-14s %6d %-6s %12.3e\n", rownames(cc)[i], cc$NGenes[i],
              cc$Direction[i], cc$PValue[i]))
cat("\n  MAPK_no_IEG still significant -> the pathway signal is independent.\n")
cat("  MAPK_no_IEG non-significant   -> 'MAPK down' was the IEG result relabelled.\n")

hdr("HOW TO READ THIS")
cat("
  Section 1  within-cohort depth correlation near zero means the pooled
             -0.531 was an artefact of combining cohorts of different depth.

  Section 3  ERK-output genes DOWN would mean the transcriptome indicates
             REDUCED ERK output, which conflicts with the raised pERK protein
             Labonte reported in females. Both can be true: pERK is
             phosphorylation state, these genes are transcriptional feedback,
             and a disrupted feedback loop can uncouple them. Report the
             transcriptional result as what it is and cite theirs as the
             protein-level counterpart rather than claiming agreement.

  Section 4  decides whether MAPK can be described as an independent finding
             or must be presented as the IEG result in pathway language.
")
cat("\noutputs in", OUT, "\n")
sessionInfo()

