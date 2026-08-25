.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(DESeq2); library(ggplot2); library(pheatmap)
  library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)})
logf <- function(...) cat(sprintf(...),"\n")
NP  <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
CS  <- file.path(NP,"counts_scz"); SHD <- file.path(CS,"sheets")
OR  <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
OUT <- file.path(CS,"HIPP_final"); dir.create(OUT, showWarnings=FALSE)
sv <- function(p,w=7,h=6){ ggsave(file.path(OUT,paste0(p,".pdf")),width=w,height=h)
                           ggsave(file.path(OUT,paste0(p,".png")),width=w,height=h,dpi=150) }

## ---------- data ----------
bt <- read.delim(file.path(OR,"gene_id2biotype.tsv"),header=FALSE,col.names=c("gid","bt"))
nm <- read.delim(file.path(OR,"gene_id2name.tsv"),  header=FALSE,col.names=c("gid","sym"))
nmv <- setNames(nm$sym, nm$gid)
cm <- read.delim(file.path(NP,"PRJNA574470/counts/PRJNA574470_gene_counts_matrix.tsv"),
                 comment.char="#", check.names=FALSE)
g <- cm[[1]]; cm[[1]] <- NULL
cm <- cm[, !colnames(cm) %in% c("Chr","Start","End","Strand","Length"), drop=FALSE]
colnames(cm) <- regmatches(colnames(cm), regexpr("SRR[0-9]+", colnames(cm)))
cnt <- as.matrix(cm); rownames(cnt) <- g; storage.mode(cnt) <- "integer"
ss <- merge(read.delim(file.path(SHD,"GSE138082_samplesheet.tsv"), check.names=FALSE),
            read.delim(file.path(SHD,"GSE138082_prua.tsv")), by="srr")
rownames(ss) <- ss$srr
k <- intersect(colnames(cnt), rownames(ss)); cnt <- cnt[,k]; ss <- ss[k,]
cnt <- cnt[rownames(cnt) %in% bt$gid[bt$bt %in% c("protein_coding","lncRNA")],]
ss$diagnosis <- relevel(factor(ss$diagnosis),"Control")
for (v in c("age","pmi","rin","prua")) ss[[paste0(v,"_z")]] <- as.numeric(scale(as.numeric(ss[[v]])))
cnt <- cnt[rowSums(cnt>=10) >= min(table(ss$diagnosis)), ]
logf("hippocampus: n=%d | %d genes", nrow(ss), nrow(cnt))

dds <- suppressMessages(DESeq(DESeqDataSetFromMatrix(cnt, ss, ~ age_z+sex+pmi_z+rin_z+diagnosis), quiet=TRUE))
res <- as.data.frame(results(dds, name=tail(resultsNames(dds),1)))
res$gene <- rownames(res); res$sym <- ifelse(is.na(nmv[res$gene]), res$gene, nmv[res$gene])
write.csv(res, file.path(OUT,"hipp_R1_full.csv"), row.names=FALSE)
sig <- subset(res, !is.na(padj) & padj<0.05)
logf("R1: padj05 = %d | >=0.3 = %d | >=1 = %d", nrow(sig),
     sum(abs(sig$log2FoldChange)>=0.3), sum(abs(sig$log2FoldChange)>=1))

## ---------- STEP 7: exploratory ----------
vsd <- vst(dds, blind=TRUE); V <- assay(vsd)
pc <- prcomp(t(V), scale.=FALSE); pv <- summary(pc)$importance[2,]*100
pd <- data.frame(PC1=pc$x[,1], PC2=pc$x[,2], ss)
ggplot(pd, aes(PC1,PC2,colour=diagnosis)) + geom_point(size=2.6) + theme_bw() +
  labs(x=sprintf("PC1 (%.1f%%)",pv[1]), y=sprintf("PC2 (%.1f%%)",pv[2]),
       title="GSE138082 hippocampus - PCA (VST)"); sv("Fig_H1_PCA")

logf("\n--- PC vs covariate association (p-values) ---")
covs <- c("diagnosis","sex","age","pmi","rin","prua")
tab <- do.call(rbind, lapply(1:5, function(i){
  sapply(covs, function(v){
    x <- ss[[v]]
    if (is.character(x)||is.factor(x)) signif(kruskal.test(pc$x[,i], factor(x))$p.value,3)
    else signif(cor.test(pc$x[,i], as.numeric(x))$p.value,3) }) }))
rownames(tab) <- paste0("PC",1:5); print(tab)
write.csv(tab, file.path(OUT,"PC_covariate_pvalues.csv"))

d <- dist(t(V)); pheatmap(as.matrix(d), annotation_col=ss[,c("diagnosis","sex"),drop=FALSE],
  show_rownames=FALSE, show_colnames=FALSE, main="Sample-distance clustering",
  filename=file.path(OUT,"Fig_H2_distance.pdf"), width=8, height=7)
md <- colMeans(as.matrix(d)); out <- names(md)[md > mean(md)+3*sd(md)]
logf("outlier screen (>3 SD mean distance): %s", if(length(out)) paste(out,collapse=", ") else "none")

## ---------- STEP 9: visualisation ----------
res$tier <- ifelse(is.na(res$padj)|res$padj>=0.05,"NS",
            ifelse(res$log2FoldChange>0,"UP","DOWN"))
ggplot(res, aes(log2FoldChange, -log10(pvalue), colour=tier)) +
  geom_point(alpha=.5, size=1) +
  scale_colour_manual(values=c(NS="grey80", UP="#C0392B", DOWN="#2471A3")) +
  geom_vline(xintercept=c(-1,1), linetype=2, colour="grey40") + theme_bw() +
  labs(title=sprintf("SCZ hippocampus R1 - %d genes FDR<0.05", nrow(sig)))
sv("Fig_H3_volcano")
ggplot(res, aes(log10(baseMean+1), log2FoldChange, colour=tier)) +
  geom_point(alpha=.5,size=1) + geom_hline(yintercept=0) +
  scale_colour_manual(values=c(NS="grey80", UP="#C0392B", DOWN="#2471A3")) +
  theme_bw() + labs(title="MA plot - SCZ hippocampus"); sv("Fig_H4_MA")

top <- head(sig[order(sig$padj),], 50)
mat <- V[top$gene,]; rownames(mat) <- top$sym
pheatmap(t(scale(t(mat))), annotation_col=ss[,c("diagnosis"),drop=FALSE],
  show_colnames=FALSE, fontsize_row=7, main="Top 50 DE genes (VST, row-scaled)",
  filename=file.path(OUT,"Fig_H5_heatmap.pdf"), width=8, height=9)

## ---------- STEP 10: functional interpretation ----------
ent <- function(x) na.omit(AnnotationDbi::mapIds(org.Hs.eg.db, sub("\\..*","",x),
                  "ENTREZID","ENSEMBL", multiVals="first"))
uni <- ent(res$gene)
for (dir in c("ALL","UP","DOWN")) {
  s <- if (dir=="ALL") sig else sig[sig$log2FoldChange*ifelse(dir=="UP",1,-1) > 0,]
  if (nrow(s) < 5) { logf("ORA %s: too few genes", dir); next }
  ge <- ent(s$gene)
  for (ont in c("BP","MF","CC")) {
    e <- enrichGO(ge, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont=ont,
                  universe=uni, pAdjustMethod="BH", qvalueCutoff=0.2)
    n <- if (is.null(e)) 0 else sum(e@result$p.adjust < 0.05)
    logf("ORA GO-%s %-4s : %d terms", ont, dir, n)
    if (n) { write.csv(e@result, file.path(OUT,sprintf("ORA_GO%s_%s.csv",ont,dir)))
             print(dotplot(e, showCategory=15)); sv(sprintf("Fig_H6_GO%s_%s",ont,dir),8,7) }
  }
  kk <- enrichKEGG(ge, organism="hsa", universe=as.character(uni), pAdjustMethod="BH")
  n <- if (is.null(kk)) 0 else sum(kk@result$p.adjust < 0.05)
  logf("ORA KEGG %-4s : %d terms", dir, n)
  if (n) write.csv(kk@result, file.path(OUT,sprintf("ORA_KEGG_%s.csv",dir)))
}
rk <- res[!is.na(res$stat),]; rk$e <- ent(rk$gene)[match(sub("\\..*","",rk$gene), names(ent(rk$gene)))]
rl <- setNames(rk$stat, rk$e); rl <- sort(rl[!is.na(names(rl))], decreasing=TRUE)
rl <- rl[!duplicated(names(rl))]
gs <- gseGO(rl, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="BP", pAdjustMethod="BH", verbose=FALSE)
n <- if (is.null(gs)) 0 else sum(gs@result$p.adjust < 0.05)
logf("GSEA GO-BP : %d terms", n)
if (n) write.csv(gs@result, file.path(OUT,"GSEA_GOBP.csv"))
logf("\nOutputs -> %s", OUT)
