## Frontal pool on UNIQUE DONORS only: GSE80655 (48, collapsed) + GSE78936 BA11 (28)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); suppressMessages(library(DESeq2))
logf <- function(...) cat(sprintf(...),"\n")
NP <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
CS <- file.path(NP,"counts_scz"); SHD <- file.path(CS,"sheets")
OR <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
RT <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/GSE80655/metadata/PRJNA319583_SraRunTable.csv"
OUT <- file.path(CS,"FRONTAL_pool_clean"); dir.create(OUT, showWarnings=FALSE)
bt <- read.delim(file.path(OR,"gene_id2biotype.tsv"), header=FALSE, col.names=c("gid","bt"))
KEEP <- bt$gid[bt$bt %in% c("protein_coding","lncRNA")]
rc <- function(p){ cm<-read.delim(p,comment.char="#",check.names=FALSE); g<-cm[[1]]; cm[[1]]<-NULL
  cm<-cm[,!colnames(cm) %in% c("Chr","Start","End","Strand","Length"),drop=FALSE]
  colnames(cm)<-regmatches(colnames(cm),regexpr("SRR[0-9]+",colnames(cm)))
  m<-as.matrix(cm); rownames(m)<-g; storage.mode(m)<-"integer"; m }

## GSE80655, collapsed to 48
s1 <- read.delim(file.path(SHD,"GSE80655_DLPFC_samplesheet.tsv"), check.names=FALSE)
s1 <- merge(s1, read.delim(file.path(SHD,"GSE80655_DLPFC_prua.tsv"))[,c("srr","prua")], by="srr")
rownames(s1) <- s1$srr
c1 <- rc(file.path(NP,"PRJNA319583/counts/PRJNA319583_gene_counts_matrix.tsv"))
k <- intersect(colnames(c1), rownames(s1)); c1 <- c1[,k]; s1 <- s1[k,]
rt <- read.csv(RT, check.names=FALSE); bio <- setNames(rt$BioSample, rt$Run)[colnames(c1)]
idx <- split(seq_len(ncol(c1)), bio)
c1 <- vapply(idx, function(i) as.integer(rowSums(c1[,i,drop=FALSE])), integer(nrow(c1)))
rownames(c1) <- rownames(rc(file.path(NP,"PRJNA319583/counts/PRJNA319583_gene_counts_matrix.tsv")))
ke <- !duplicated(bio); s1 <- s1[ke,]; rownames(s1) <- bio[ke]; c1 <- c1[, rownames(s1)]

## GSE78936 BA11 only
s2 <- read.delim(file.path(SHD,"GSE78936_samplesheet.tsv"), check.names=FALSE)
s2 <- merge(s2, read.delim(file.path(SHD,"GSE78936_prua.tsv"))[,c("srr","prua")], by="srr")
rownames(s2) <- s2$srr; s2 <- s2[s2$brain_region=="Brodmann's Area 11",]
c2 <- rc(file.path(NP,"PRJNA314463/counts/PRJNA314463_gene_counts.txt"))
k <- intersect(colnames(c2), rownames(s2)); c2 <- c2[,k]; s2 <- s2[k,]

g <- intersect(rownames(c1), rownames(c2)); g <- g[g %in% KEEP]
cnt <- cbind(c1[g,], c2[g,])
mt <- data.frame(diagnosis=c(s1$diagnosis, s2$diagnosis),
                 study=c(rep("GSE80655",ncol(c1)), rep("GSE78936",ncol(c2))),
                 prua=c(s1$prua, s2$prua), row.names=colnames(cnt))
mt$diagnosis <- relevel(factor(mt$diagnosis),"Control"); mt$study <- factor(mt$study)
mt$prua_z <- ave(mt$prua, mt$study, FUN=function(x) as.numeric(scale(x)))
cnt <- cnt[rowSums(cnt>=10) >= min(table(mt$diagnosis)),]
logf("CLEAN FRONTAL POOL: n=%d donors (%s) | %d genes", nrow(mt),
     paste(names(table(mt$diagnosis)), table(mt$diagnosis), sep="=", collapse=" "), nrow(cnt))
print(table(mt$study, mt$diagnosis))

fit <- function(form,lab,tag){
  dds <- suppressMessages(DESeq(DESeqDataSetFromMatrix(cnt, mt, form), quiet=TRUE))
  r <- as.data.frame(results(dds, name=tail(resultsNames(dds),1)))
  write.csv(r, file.path(OUT,paste0(tag,".csv")))
  f <- function(l) sum(!is.na(r$padj) & r$padj<0.05 & abs(r$log2FoldChange)>=l)
  logf("  %-32s padj05:%5d | >=0.3:%5d | >=0.5:%5d | >=1.0:%4d", lab, f(0), f(0.3), f(0.5), f(1))
  rownames(r)[!is.na(r$padj) & r$padj<0.05]
}
logf("\n=========== POOLED FRONTAL, UNIQUE DONORS ===========")
p0 <- fit(~ study + diagnosis,            "P0 ~study+diagnosis",   "pool_P0")
p1 <- fit(~ study + prua_z + diagnosis,   "P1 +PRUA",              "pool_P1")
logf("\nOutputs -> %s", OUT)
