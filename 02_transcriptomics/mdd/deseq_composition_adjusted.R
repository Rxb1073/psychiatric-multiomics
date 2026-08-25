#!/usr/bin/env Rscript
# MDD composition-adjusted DE (mediation ladder) — closes the cross-arm gap.
# Ports scz_composition_adjusted.R / ASD pooled_DE_composition_adjusted.R to MDD.
# MEDIATION (astrocyte) vs CONFOUND (quality) — same interpretive frame as ASD/SCZ.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(library(DESeq2))
logf <- function(...) cat(sprintf(...), "\n"); CASE <- "MDD"
M   <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics"
DE  <- file.path(M,"DE"); DEC <- file.path(M,"Deconvolution")
OR  <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
OUT <- file.path(DE,"COMPOSITION_adjusted"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
SIGNATURE <- "MB"; if (nzchar(Sys.getenv("DECONV_SIG"))) SIGNATURE <- Sys.getenv("DECONV_SIG")
logf("=== MDD COMPOSITION-ADJUSTED DE | signature: %s ===", SIGNATURE)
bt <- read.delim(file.path(OR,"gene_id2biotype.tsv"), header=FALSE, col.names=c("gid","bt")); btv <- setNames(bt$bt,bt$gid)
len <- read.delim(file.path(OR,"gene_lengths.tsv"), header=FALSE, col.names=c("gid","len")); lenv <- setNames(as.numeric(len$len),len$gid)
gm <- read.delim(file.path(OR,"gene_id2name.tsv"), header=FALSE, col.names=c("gid","nm")); nmv <- setNames(gm$nm,gm$gid)
# ---- counts + design ----
cm <- read.delim(file.path(DE,"merged_counts.tsv"), check.names=FALSE, stringsAsFactors=FALSE)
gid <- cm[[1]]; cm[[1]] <- NULL; rownames(cm) <- gid
cnt <- as.matrix(cm); cnt[is.na(cnt)] <- 0; storage.mode(cnt) <- "integer"
mt <- read.delim(file.path(DE,"rna_colData.tsv"), check.names=FALSE, stringsAsFactors=FALSE)
names(mt)[names(mt)=="sample-id"] <- "srr"; rownames(mt) <- mt$srr
k <- intersect(colnames(cnt), rownames(mt)); cnt <- cnt[,k,drop=FALSE]; mt <- mt[k,,drop=FALSE]
g <- rownames(cnt); g <- g[!is.na(btv[g]) & btv[g] %in% c("protein_coding","lncRNA")]
cnt <- cnt[g,,drop=FALSE]
logf("pooled MDD: %d samples x %d genes", ncol(cnt), nrow(cnt))
# ---- quality proxies ----
lnc <- g[btv[g]=="lncRNA"]; pc <- g[btv[g]=="protein_coding"]
mt$q_lnc <- colSums(cnt[lnc,,drop=FALSE])/colSums(cnt)
Lp <- lenv[pc]; ok <- !is.na(Lp)&Lp>0; W <- cnt[pc,,drop=FALSE][ok,,drop=FALSE]; Lo <- Lp[ok]
mt$q_len <- apply(W,2,function(x){s<-sum(x); if(is.na(s)||s==0) NA_real_ else sum(x*Lo)/s})
# ---- proportions ----
props <- read.csv(file.path(DEC,"cell_type_proportions_long.csv"), stringsAsFactors=FALSE, check.names=FALSE)
names(props)[names(props)=="sample-id"] <- "srr"
props <- subset(props, signature==SIGNATURE)
wide <- reshape(props[,c("srr","cell_type","proportion")], idvar="srr", timevar="cell_type", direction="wide")
names(wide) <- sub("^proportion\\.","",names(wide)); rownames(wide) <- wide$srr
if (!"Astrocytes" %in% names(wide)) stop("Astrocytes absent from deconv output")
mt$astro <- wide[rownames(mt),"Astrocytes"]
ncol_ <- if ("Neurons" %in% names(wide)) "Neurons" else "Excitatory"
mt$neuron <- wide[rownames(mt),ncol_]; logf("neuronal column used: %s", ncol_)
zwithin <- function(x,grp) ave(x,grp,FUN=function(v) as.numeric(scale(v)))
mt$astro_z <- zwithin(mt$astro,mt$study); mt$neuron_z <- zwithin(mt$neuron,mt$study)
mt$q_len_z <- zwithin(mt$q_len,mt$study); mt$q_lnc_z <- zwithin(mt$q_lnc,mt$study)
mt$disease <- relevel(factor(mt$disease),"Control"); mt$study <- factor(mt$study); mt$sex <- factor(mt$sex)
mt$pmi_z <- as.numeric(scale(suppressWarnings(as.numeric(mt$pmi))))
mt$ph_z  <- as.numeric(scale(suppressWarnings(as.numeric(mt$ph))))
logf("\ncor(astro_z,q_len_z)=%+.3f | cor(astro_z,q_lnc_z)=%+.3f",
  cor(mt$astro_z,mt$q_len_z,use="complete.obs"), cor(mt$astro_z,mt$q_lnc_z,use="complete.obs"))
logf("astrocyte z by disease: Control %.3f | MDD %.3f (Wilcoxon p=%.4f)",
  median(mt$astro_z[mt$disease=="Control"],na.rm=TRUE), median(mt$astro_z[mt$disease=="MDD"],na.rm=TRUE),
  suppressWarnings(wilcox.test(astro_z~disease,data=mt)$p.value))
logf("neuron z by disease   : Control %.3f | MDD %.3f (Wilcoxon p=%.4f)",
  median(mt$neuron_z[mt$disease=="Control"],na.rm=TRUE), median(mt$neuron_z[mt$disease=="MDD"],na.rm=TRUE),
  suppressWarnings(wilcox.test(neuron_z~disease,data=mt)$p.value))
keep <- complete.cases(mt[,c("disease","study","sex","pmi_z","ph_z","astro_z","neuron_z","q_len_z","q_lnc_z")])
if (any(!keep)) logf("dropping %d samples with NA in model vars", sum(!keep))
cnt <- cnt[,rownames(mt)[keep],drop=FALSE]; mt <- mt[keep,]
sg <- min(table(mt$disease)); cnt <- cnt[rowSums(cnt>=10)>=sg,,drop=FALSE]
logf("after filtering: %d samples x %d genes\n", ncol(cnt), nrow(cnt))
cat("study x disease:\n"); print(table(mt$study, mt$disease))
fit <- function(form,lab){
  dds <- suppressMessages(DESeq(DESeqDataSetFromMatrix(cnt, mt, form), quiet=TRUE))
  r <- results(dds, contrast=c("disease",CASE,"Control"), alpha=0.05)
  s <- subset(as.data.frame(r), !is.na(padj)&padj<0.05)
  logf("%-56s : %4d genes (up=%d, down=%d)", lab, nrow(s), sum(s$log2FoldChange>0), sum(s$log2FoldChange<0))
  list(res=r, sig=rownames(s)) }
logf("\n=========== POOLED MDD MODELS (n=%d) ===========", ncol(cnt))
m0 <- fit(~ study + sex + pmi_z + ph_z + disease,                              "M0 baseline (study+sex+pmi+ph)")
m1 <- fit(~ study + sex + pmi_z + ph_z + astro_z + disease,                    "M1 + astrocyte  [MEDIATION]")
m2 <- fit(~ study + sex + pmi_z + ph_z + q_len_z + q_lnc_z + disease,          "M2 + quality proxies [CONFOUND]")
m3 <- fit(~ study + sex + pmi_z + ph_z + astro_z + q_len_z + q_lnc_z + disease,"M3 + both")
m4 <- fit(~ study + sex + pmi_z + ph_z + astro_z + neuron_z + disease,         "M4 + astrocyte + neuronal [COMPOSITION]")
base <- m0$sig
logf("\n=========== FATE OF THE %d BASELINE GENES ===========", length(base))
for (nm in c("M1 astro","M2 quality","M3 both","M4 composition")) {
  s <- switch(nm,"M1 astro"=m1$sig,"M2 quality"=m2$sig,"M3 both"=m3$sig,"M4 composition"=m4$sig)
  logf("  retained under %-16s : %4d / %d (%.0f%%)", nm, length(intersect(base,s)), length(base),
       100*length(intersect(base,s))/max(1,length(base))) }
surv <- Reduce(intersect, list(base,m1$sig,m2$sig,m3$sig,m4$sig))
logf("\n  survive EVERY model: %d", length(surv))
if (length(surv)) { r <- as.data.frame(m3$res)[surv,]
  r$gene <- ifelse(is.na(nmv[surv]), surv, nmv[surv]); r$biotype <- btv[surv]
  print(head(r[order(r$padj), c("gene","biotype","log2FoldChange","padj")], 30)) }
IEG <- c("NPAS4","FOS","FOSB","EGR1","EGR2","NR4A1","NR4A2","NR4A3","ATF3","BTG2","ARC","DUSP1","JUNB")
ieg_ids <- names(nmv)[nmv %in% IEG]
logf("\n  --- IEG fate (padj per model) ---")
logf("  %-8s %8s %8s %8s %8s %9s", "gene","M0","M1","M2","M3","LFC(M3)")
for (id in ieg_ids) {
  if (!id %in% rownames(m0$res)) next
  p <- function(m){ v <- as.data.frame(m$res)[id,"padj"]; if (is.na(v)) "  NA" else sprintf("%.4f", v) }
  l3 <- as.data.frame(m3$res)[id,"log2FoldChange"]
  logf("  %-8s %8s %8s %8s %8s %9.3f", nmv[id], p(m0), p(m1), p(m2), p(m3), l3) }
write.table(data.frame(gene_id=base, gene=nmv[base], biotype=btv[base],
  in_M1=base%in%m1$sig, in_M2=base%in%m2$sig, in_M3=base%in%m3$sig, in_M4=base%in%m4$sig),
  file.path(OUT,sprintf("baseline_fate_%s.tsv",SIGNATURE)), sep="\t", quote=FALSE, row.names=FALSE)
saveRDS(list(m0=m0,m1=m1,m2=m2,m3=m3,m4=m4,meta=mt), file.path(OUT,sprintf("models_%s.rds",SIGNATURE)))
logf("\n=========== HOW TO READ THIS ===========")
logf(" M1 kills genes -> MDD signal MEDIATED by composition (as in ASD / SCZ-frontal).")
logf(" M1 does NOT kill them -> MDD signal is cell-intrinsic, confirming the arm's original conclusion")
logf("    on the same evidential footing as the other two arms.")
sink(file.path(OUT,"sessionInfo.txt")); print(sessionInfo()); sink()
logf("\noutputs -> %s", OUT)
