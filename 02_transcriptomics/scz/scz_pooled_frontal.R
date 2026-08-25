#!/usr/bin/env Rscript
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(library(DESeq2))
HAVE_SVA <- requireNamespace("sva", quietly = TRUE)
if (HAVE_SVA) suppressMessages(library(sva))
logf <- function(...) cat(sprintf(...), "\n"); CASE <- "SCZ"
NEURO <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
CS <- file.path(NEURO,"counts_scz"); SHEETS <- file.path(CS,"sheets")
OUT <- file.path(CS,"POOLED_frontal"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
ASD_DE <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
GTFMAP <- file.path(ASD_DE,"gene_id2name.tsv"); BIOTYPE <- file.path(ASD_DE,"gene_id2biotype.tsv")
KEEP_BIOTYPES <- c("protein_coding","lncRNA")
logf("=== SCZ FRONTAL POOLED DE (covariate route MAIN) ===")
read_counts <- function(path) {
  cm <- read.delim(path, comment.char="#", check.names=FALSE, stringsAsFactors=FALSE)
  gid <- cm[[1]]; cm[[1]] <- NULL
  cm <- cm[, !(colnames(cm) %in% c("Chr","Start","End","Strand","Length")), drop=FALSE]
  hit <- regmatches(colnames(cm), regexpr("SRR[0-9]+", colnames(cm)))
  cm <- cm[, nchar(hit)>0, drop=FALSE]
  colnames(cm) <- regmatches(colnames(cm), regexpr("SRR[0-9]+", colnames(cm)))
  rownames(cm) <- gid; m <- as.matrix(cm); m[is.na(m)] <- 0; storage.mode(m) <- "integer"; m }
datasets <- list(
  GSE80655=list(counts=file.path(NEURO,"PRJNA319583/counts/PRJNA319583_gene_counts_matrix.tsv"),
                sheet=file.path(SHEETS,"GSE80655_DLPFC_samplesheet.tsv")),
  GSE78936=list(counts=file.path(NEURO,"PRJNA314463/counts/PRJNA314463_gene_counts.txt"),
                sheet=file.path(SHEETS,"GSE78936_samplesheet.tsv")))
mats <- list(); metas <- list()
for (nm in names(datasets)) {
  d <- datasets[[nm]]; ss <- read.delim(d$sheet, stringsAsFactors=FALSE, check.names=FALSE)
  rownames(ss) <- ss$srr; cm <- read_counts(d$counts)
  common <- intersect(colnames(cm), rownames(ss)); cm <- cm[,common,drop=FALSE]; ss <- ss[common,,drop=FALSE]
  ss$study <- nm; ss <- ss[,c("srr","diagnosis","brain_region","study")]
  logf("%-10s : %3d samples | %s", nm, ncol(cm), paste(names(table(ss$brain_region)),table(ss$brain_region),sep="=",collapse=", "))
  mats[[nm]] <- cm; metas[[nm]] <- ss }
genes <- Reduce(intersect, lapply(mats,rownames))
bt <- read.delim(BIOTYPE, header=FALSE, col.names=c("gid","biotype")); btv <- setNames(bt$biotype,bt$gid)
genes <- genes[!is.na(btv[genes]) & btv[genes] %in% KEEP_BIOTYPES]
logf("common genes after biotype filter: %d", length(genes))
counts <- do.call(cbind, lapply(mats, function(m) m[genes,,drop=FALSE]))
meta <- do.call(rbind, metas); rownames(meta) <- meta$srr; meta <- meta[colnames(counts),,drop=FALSE]
stopifnot(identical(colnames(counts), rownames(meta)))
meta$diagnosis <- relevel(factor(meta$diagnosis), ref="Control")
meta$study <- factor(meta$study); meta$brain_region <- factor(meta$brain_region)
logf("\nPOOLED FRONTAL: %d samples x %d genes", ncol(counts), nrow(counts))
cat("\nstudy x diagnosis:\n"); print(table(meta$study, meta$diagnosis))
cat("\nbrain_region x study:\n"); print(table(meta$brain_region, meta$study))
sg <- min(table(meta$diagnosis)); keep <- rowSums(counts>=10)>=sg
logf("\nlow-count prefilter: %d -> %d genes (smallest group n=%d)", nrow(counts), sum(keep), sg)
counts <- counts[keep,,drop=FALSE]
mm <- model.matrix(~ study + brain_region + diagnosis, data=meta)
full_rank <- (qr(mm)$rank==ncol(mm))
if (full_rank) { MAIN_DESIGN <- ~ study + brain_region + diagnosis; logf("design: ~study+brain_region+diagnosis (full rank)")
} else { MAIN_DESIGN <- ~ study + diagnosis
  logf("WARN: ~study+brain_region+diagnosis rank-deficient; fallback ~study+diagnosis (region=limitation)") }
annot <- function(ids){ gn<-setNames(rep(NA_character_,length(ids)),ids)
  if(file.exists(GTFMAP)){m<-read.delim(GTFMAP,header=FALSE,col.names=c("gid","name"));gn[]<-setNames(m$name,m$gid)[ids]}
  list(name=gn, biotype=btv[ids]) }
run_de <- function(cnt, md, design, label, combat=FALSE) {
  if (combat) { logf("\n--- %s: ComBat-seq(batch=study)->~diagnosis ---", label)
    adj <- ComBat_seq(counts=cnt, batch=md$study, group=md$diagnosis)
    dds <- DESeqDataSetFromMatrix(adj, colData=md, design=~ diagnosis)
  } else { logf("\n--- %s: DESeq2 %s ---", label, deparse(design))
    dds <- DESeqDataSetFromMatrix(cnt, colData=md, design=design) }
  dds <- DESeq(dds)
  cf <- grep(sprintf("diagnosis_%s_vs_Control",CASE), resultsNames(dds), value=TRUE)
  res <- results(dds, contrast=c("diagnosis",CASE,"Control"), alpha=0.05)
  styp <- if(requireNamespace("apeglm",quietly=TRUE))"apeglm" else "normal"
  resL <- lfcShrink(dds, coef=cf, type=styp); a <- annot(rownames(res))
  o <- data.frame(gene_id=rownames(res), gene_name=a$name, biotype=a$biotype,
    baseMean=res$baseMean, log2FC_shrunk=resL$log2FoldChange, pvalue=res$pvalue, padj=res$padj, row.names=NULL)
  o[order(o$padj),] }
sig <- function(o,lfc=0.3) subset(o, !is.na(padj)&padj<0.05&abs(log2FC_shrunk)>lfc)
main <- run_de(counts, meta, MAIN_DESIGN, "MAIN (covariate route)")
write.table(main, file.path(OUT,"POOLED_frontal_MAIN_full_results.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
for (l in c(0,0.3,0.5,1)){s<-sig(main,l); logf("  MAIN padj05 absLFC>%.1f : %5d (up=%d down=%d)",l,nrow(s),sum(s$log2FC_shrunk>0),sum(s$log2FC_shrunk<0))}
write.table(sig(main,0.3), file.path(OUT,"POOLED_frontal_MAIN_sig_padj05_absLFC0.3.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
if (HAVE_SVA) { cc <- run_de(counts, meta, NULL, "CROSS-CHECK (ComBat-seq)", combat=TRUE)
  write.table(cc, file.path(OUT,"POOLED_frontal_ComBat_full_results.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  gM<-sig(main,0.3)$gene_id; gC<-sig(cc,0.3)$gene_id; ov<-intersect(gM,gC)
  logf("\n=== CONCORDANCE main vs ComBat (padj05,|LFC|>0.3) ===")
  logf("  main=%d combat=%d overlap=%d (%.1f%% of main)",length(gM),length(gC),length(ov),ifelse(length(gM)>0,100*length(ov)/length(gM),NA))
  writeLines(ov, file.path(OUT,"concordant_main_combat.txt"))
} else logf("\nWARN: sva missing - skip ComBat")
logf("\n=== LEAVE-ONE-DATASET-OUT (main model) ===")
core <- sig(main,0.3)$gene_id
for (drop in names(datasets)) {
  ok <- meta$study!=drop; md <- droplevels(meta[ok,]); cn <- counts[,ok,drop=FALSE]
  dsn <- if (nlevels(md$study)<2) {
           if (nlevels(droplevels(md$brain_region))>1) {
             mmb <- model.matrix(~ brain_region + diagnosis, md)
             if (qr(mmb)$rank==ncol(mmb)) ~ brain_region+diagnosis else ~ diagnosis
           } else ~ diagnosis
         } else MAIN_DESIGN
  cn <- cn[rowSums(cn>=10)>=min(table(md$diagnosis)),,drop=FALSE]
  r <- tryCatch(run_de(cn,md,dsn,sprintf("LODO drop %s",drop)), error=function(e){logf("  drop %s FAILED: %s",drop,conditionMessage(e));NULL})
  if(!is.null(r)){rem<-sig(r,0.3)$gene_id; kc<-intersect(core,rem)
    logf("  drop %-10s: %d sig | core retained %d/%d (%.0f%%)",drop,length(rem),length(kc),length(core),ifelse(length(core)>0,100*length(kc)/length(core),0))
    write.table(sig(r,0.3), file.path(OUT,sprintf("LODO_drop_%s_sig.tsv",drop)), sep="\t", quote=FALSE, row.names=FALSE)} }
logf("\n=== BETWEEN-DATASET DIRECTION CONCORDANCE ===")
per <- list()
for (nm in names(datasets)) { ok<-meta$study==nm; md<-droplevels(meta[ok,]); cn<-counts[,ok,drop=FALSE]
  cn <- cn[rowSums(cn>=10)>=min(table(md$diagnosis)),,drop=FALSE]
  dsn <- if (nlevels(droplevels(md$brain_region))>1) {
           mmb <- model.matrix(~ brain_region+diagnosis, md)
           if (qr(mmb)$rank==ncol(mmb)) ~ brain_region+diagnosis else ~ diagnosis
         } else ~ diagnosis
  per[[nm]] <- tryCatch(run_de(cn,md,dsn,sprintf("per %s",nm)), error=function(e) NULL) }
if (all(!sapply(per,is.null))) { a<-per[[1]]; b<-per[[2]]
  m<-merge(a[,c("gene_id","log2FC_shrunk","padj")], b[,c("gene_id","log2FC_shrunk","padj")], by="gene_id", suffixes=c(".a",".b"))
  r_all<-cor(m$log2FC_shrunk.a,m$log2FC_shrunk.b,use="complete.obs")
  msig<-subset(m,(!is.na(padj.a)&padj.a<0.05)|(!is.na(padj.b)&padj.b<0.05))
  r_sig<-if(nrow(msig)>2)cor(msig$log2FC_shrunk.a,msig$log2FC_shrunk.b,use="complete.obs") else NA
  same<-if(nrow(msig)>0)mean(sign(msig$log2FC_shrunk.a)==sign(msig$log2FC_shrunk.b)) else NA
  logf("  LFC corr (all genes)      : r=%.3f", r_all)
  logf("  LFC corr (sig either)     : r=%.3f (n=%d)", r_sig, nrow(msig))
  logf("  same-direction (sig either): %.1f%%", 100*same)
  write.table(m, file.path(OUT,"between_dataset_LFC.tsv"), sep="\t", quote=FALSE, row.names=FALSE) }
pdf(file.path(OUT,"POOLED_frontal_QC.pdf"), width=7, height=6)
hist(main$pvalue, breaks=50, col="grey", main="Frontal pool MAIN p-values", xlab="p"); dev.off()
sink(file.path(OUT,"sessionInfo.txt")); print(sessionInfo()); sink()
logf("\noutputs -> %s", OUT); logf("=== done ===")
