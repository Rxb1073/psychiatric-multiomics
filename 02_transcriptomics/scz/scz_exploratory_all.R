## STEP 7 - exploratory analysis for all three SCZ cohorts
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(DESeq2); library(ggplot2); library(pheatmap)})
logf <- function(...) cat(sprintf(...),"\n")
NP  <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
CS  <- file.path(NP,"counts_scz"); SHD <- file.path(CS,"sheets")
OR  <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
RT  <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/GSE80655/metadata/PRJNA319583_SraRunTable.csv"
OUT <- file.path(CS,"EXPLORATORY"); dir.create(OUT, showWarnings=FALSE)
bt  <- read.delim(file.path(OR,"gene_id2biotype.tsv"), header=FALSE, col.names=c("gid","bt"))
KEEP <- bt$gid[bt$bt %in% c("protein_coding","lncRNA")]

rc <- function(p){ cm <- read.delim(p, comment.char="#", check.names=FALSE)
  g <- cm[[1]]; cm[[1]] <- NULL
  cm <- cm[, !colnames(cm) %in% c("Chr","Start","End","Strand","Length"), drop=FALSE]
  colnames(cm) <- regmatches(colnames(cm), regexpr("SRR[0-9]+", colnames(cm)))
  m <- as.matrix(cm); rownames(m) <- g; storage.mode(m) <- "integer"; m }

explore <- function(cnt, ss, lab, covs) {
  ss$diagnosis <- droplevels(factor(ss$diagnosis))
  cnt <- cnt[rownames(cnt) %in% KEEP, , drop=FALSE]
  cnt <- cnt[rowSums(cnt>=10) >= min(table(ss$diagnosis)), , drop=FALSE]
  logf("\n===== %s : n=%d (%s) | %d genes =====", lab, nrow(ss),
       paste(names(table(ss$diagnosis)), table(ss$diagnosis), sep="=", collapse=" "), nrow(cnt))
  dds <- DESeqDataSetFromMatrix(cnt, ss, ~1)
  V <- assay(vst(dds, blind=TRUE))
  pc <- prcomp(t(V)); pv <- summary(pc)$importance[2,]*100
  pd <- data.frame(PC1=pc$x[,1], PC2=pc$x[,2], ss)
  p <- ggplot(pd, aes(PC1,PC2,colour=diagnosis)) + geom_point(size=2.6) + theme_bw() +
    labs(x=sprintf("PC1 (%.1f%%)",pv[1]), y=sprintf("PC2 (%.1f%%)",pv[2]), title=paste(lab,"- PCA (VST)"))
  if ("brain_region" %in% names(ss) && length(unique(ss$brain_region))>1)
    p <- p + aes(shape=brain_region)
  ggsave(file.path(OUT,paste0("PCA_",lab,".pdf")), p, width=7, height=6)
  ggsave(file.path(OUT,paste0("PCA_",lab,".png")), p, width=7, height=6, dpi=150)

  cv <- intersect(covs, names(ss))
  cv <- cv[sapply(cv, function(v){ x<-ss[[v]]
       if(is.character(x)||is.factor(x)) length(unique(x))>1 && length(unique(x))<nrow(ss)
       else sum(!is.na(suppressWarnings(as.numeric(as.character(x)))))>=10 })]
  tab <- do.call(rbind, lapply(1:min(5,ncol(pc$x)), function(i)
    sapply(cv, function(v){ x <- ss[[v]]
      if (is.character(x)||is.factor(x)) signif(kruskal.test(pc$x[,i], factor(x))$p.value,3)
      else { y <- suppressWarnings(as.numeric(as.character(x)))
             signif(cor.test(pc$x[,i], y)$p.value,3) } })))
  rownames(tab) <- paste0("PC",1:nrow(tab))
  logf("-- PC vs covariate association (p) --"); print(tab)
  write.csv(tab, file.path(OUT,paste0("PCcov_",lab,".csv")))
  logf("   variance explained PC1-5: %s", paste(sprintf("%.1f%%",pv[1:5]), collapse=" "))

  d <- dist(t(V))
  ann <- ss[, intersect(c("diagnosis","brain_region","sex"), names(ss)), drop=FALSE]
  ann <- ann[, sapply(ann, function(x) !all(is.na(x)) && length(unique(x[!is.na(x)])) > 1), drop=FALSE]
  for (j in seq_len(ncol(ann))) ann[[j]] <- factor(as.character(ann[[j]]))
  if (!ncol(ann)) ann <- NA
  pheatmap(as.matrix(d), annotation_col=ann, show_rownames=FALSE, show_colnames=FALSE,
    main=paste(lab,"- sample distances"),
    filename=file.path(OUT,paste0("dist_",lab,".pdf")), width=8, height=7)
  md <- colMeans(as.matrix(d)); o <- names(md)[md > mean(md)+3*sd(md)]
  logf("   outlier screen (>3 SD): %s", if(length(o)) paste(o,collapse=", ") else "none")
}

## ---- GSE80655, collapsed to 48 ----
ss <- read.delim(file.path(SHD,"GSE80655_DLPFC_samplesheet.tsv"), check.names=FALSE)
pr <- read.delim(file.path(SHD,"GSE80655_DLPFC_prua.tsv")); ss <- merge(ss,pr[,c("srr","prua")],by="srr")
rownames(ss) <- ss$srr
cnt <- rc(file.path(NP,"PRJNA319583/counts/PRJNA319583_gene_counts_matrix.tsv"))
k <- intersect(colnames(cnt), rownames(ss)); cnt <- cnt[,k]; ss <- ss[k,]
mp <- setNames(read.csv(RT, check.names=FALSE)$BioSample, read.csv(RT, check.names=FALSE)$Run)
bio <- mp[colnames(cnt)]
idx <- split(seq_len(ncol(cnt)), bio)
c2 <- vapply(idx, function(i) as.integer(rowSums(cnt[,i,drop=FALSE])), integer(nrow(cnt)))
rownames(c2) <- rownames(cnt)
ke <- !duplicated(bio); ss <- ss[ke,]; rownames(ss) <- bio[ke]; cnt <- c2[, rownames(ss)]
logf("GSE80655 collapsed: %d -> %d", length(bio), ncol(cnt))
explore(cnt, ss, "GSE80655_DLPFC_n48", c("diagnosis","sex","age","pmi","ph","prua"))

## ---- GSE78936, full 52 and BA11-only ----
ss <- read.delim(file.path(SHD,"GSE78936_samplesheet.tsv"), check.names=FALSE)
pr <- read.delim(file.path(SHD,"GSE78936_prua.tsv")); ss <- merge(ss,pr[,c("srr","prua")],by="srr")
rownames(ss) <- ss$srr
cnt <- rc(file.path(NP,"PRJNA314463/counts/PRJNA314463_gene_counts.txt"))
k <- intersect(colnames(cnt), rownames(ss)); cnt <- cnt[,k]; ss <- ss[k,]
explore(cnt, ss, "GSE78936_all52", c("diagnosis","brain_region","prua"))
b <- ss$brain_region=="Brodmann's Area 11"
explore(cnt[,rownames(ss)[b]], ss[b,], "GSE78936_BA11_n28", c("diagnosis","prua"))

## ---- GSE138082 hippocampus ----
ss <- read.delim(file.path(SHD,"GSE138082_samplesheet.tsv"), check.names=FALSE)
pr <- read.delim(file.path(SHD,"GSE138082_prua.tsv")); ss <- merge(ss,pr[,c("srr","prua")],by="srr")
rownames(ss) <- ss$srr
cnt <- rc(file.path(NP,"PRJNA574470/counts/PRJNA574470_gene_counts_matrix.tsv"))
k <- intersect(colnames(cnt), rownames(ss)); cnt <- cnt[,k]; ss <- ss[k,]
explore(cnt, ss, "GSE138082_hipp_n78", c("diagnosis","sex","age","pmi","rin","prua"))
logf("\nOutputs -> %s", OUT)
