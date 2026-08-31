.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(DESeq2); library(org.Hs.eg.db); library(AnnotationDbi)})
logf <- function(...) cat(sprintf(...),"\n")
AUT <- "/rds/projects/e/elhamsak-dissertation-autism"
AUT2<- "/rds/projects/e/elhamsak-dissertation-autism-2"
STG <- file.path(AUT,"RNA_seq/_pheno_staging")
DEC <- file.path(AUT2,"RNA_seq/Deconvolution/cell_type_proportions_MB.csv")
OR  <- file.path(AUT2,"RNA_seq/DE_analysis")
OUT <- file.path(OR,"DISSERTATION_final/Transcriptomics_ASD/tables/nocomposition")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
QC_FAIL <- c("SRR5938468","SRR5938447","SRR5938445","SRR1510129",
             "SRR5938444","SRR5938441","SRR1004284","SRR5938466")
cohorts <- list(
  GSE102741   = list(c=file.path(AUT ,"RNA_seq/GSE102741/counts/GSE102741_gene_counts.txt"),
                     s=file.path(STG,"GSE102741_samplesheet.tsv")),
  PRJNA434002 = list(c=file.path(AUT2,"RNA_seq/PRJNA434002_bulk/counts/PRJNA434002_150bp_gene_counts.txt"),
                     s=file.path(STG,"PRJNA434002_samplesheet.tsv")))
read_fc <- function(p){ f<-readLines(p,n=1); sk<-if(grepl("^# ",f))1 else 0
  m<-read.delim(p,header=TRUE,skip=sk,check.names=FALSE,row.names=1)
  d<-intersect(c("Chr","Start","End","Strand","Length"),colnames(m))
  if(length(d)) m<-m[,setdiff(colnames(m),d),drop=FALSE]
  colnames(m)<-regmatches(colnames(m),regexpr("[SED]RR[0-9]+",colnames(m))); as.matrix(m) }
ndx <- function(x){ x<-tolower(trimws(x))
  ifelse(grepl("asd|autism",x),"ASD", ifelse(grepl("control|ctrl|con\\b",x),"Control","OTHER")) }
prop <- read.csv(DEC, stringsAsFactors=FALSE); prop$srr <- gsub('"','',prop$srr)
mats <- list(); metas <- list()
for(nm in names(cohorts)){
  cm <- read_fc(cohorts[[nm]]$c)
  ss <- read.delim(cohorts[[nm]]$s, check.names=FALSE, stringsAsFactors=FALSE)
  idc<-grep("^srr$",names(ss),ignore.case=TRUE,value=TRUE)[1]
  dxc<-grep("diagnos",names(ss),ignore.case=TRUE,value=TRUE)[1]
  sxc<-grep("^sex$|gender",names(ss),ignore.case=TRUE,value=TRUE)[1]
  agc<-grep("^age$",names(ss),ignore.case=TRUE,value=TRUE)[1]
  ss$.id <- regmatches(ss[[idc]], regexpr("[SED]RR[0-9]+", ss[[idc]]))
  md <- data.frame(srr=ss$.id, series=nm, diagnosis=ndx(ss[[dxc]]),
                   sex=tolower(substr(trimws(ss[[sxc]]),1,1)),
                   age=suppressWarnings(as.numeric(as.character(ss[[agc]]))), stringsAsFactors=FALSE)
  common <- intersect(colnames(cm), md$srr)
  mats[[nm]] <- cm[,common]; metas[[nm]] <- md[match(common, md$srr),] }
genes  <- Reduce(intersect, lapply(mats, rownames))
counts <- do.call(cbind, lapply(mats, function(x) x[genes,,drop=FALSE]))
meta   <- do.call(rbind, metas); rownames(meta) <- meta$srr; meta <- meta[colnames(counts),]
meta$Astro <- prop$Astrocytes[match(meta$srr, prop$srr)]
meta$Neuro <- prop$Neurons[match(meta$srr, prop$srr)]
ok <- meta$diagnosis %in% c("ASD","Control") & meta$sex %in% c("m","f") &
      !is.na(meta$age) & !is.na(meta$Astro) & meta$age < 18 & !(meta$srr %in% QC_FAIL)
counts <- round(counts[,ok]); m <- meta[ok,]
pv <- function(mm) suppressWarnings(wilcox.test(age~diagnosis, mm)$p.value)
logf("pre-match n=%d age p=%.4g", nrow(m), pv(m))
while(pv(m) <= 0.05){
  medA<-median(m$age[m$diagnosis=="ASD"]); medC<-median(m$age[m$diagnosis=="Control"])
  if(sum(m$diagnosis=="ASD")<=8 || sum(m$diagnosis=="Control")<=8) break
  i <- if(medA>=medC) which(m$diagnosis=="ASD") else which(m$diagnosis=="Control")
  d <- if(medA>=medC) i[which.max(m$age[i])] else i[which.min(m$age[i])]
  m <- m[-d,,drop=FALSE] }
counts <- counts[, rownames(m)]
logf("post-match n=%d (ASD=%d Ctrl=%d) age p=%.4g", nrow(m), sum(m$diagnosis=="ASD"),
     sum(m$diagnosis=="Control"), pv(m))
m$diagnosis <- factor(as.character(m$diagnosis), levels=c("Control","ASD"))
m$series <- droplevels(factor(m$series)); m$sex <- droplevels(factor(m$sex))
for(v in c("Astro","Neuro")) m[[paste0(v,"_z")]] <- as.numeric(scale(m[[v]]))
bt <- read.delim(file.path(OR,"gene_id2biotype.tsv"), header=FALSE, col.names=c("gid","bt"))
KEEP <- bt$gid[bt$bt %in% c("protein_coding","lncRNA")]
counts <- counts[rownames(counts) %in% KEEP, , drop=FALSE]
counts <- counts[rowSums(counts>=10) >= min(table(m$diagnosis)), ]
logf("genes after filtering: %d", nrow(counts))
run <- function(form,label){
  dds <- DESeqDataSetFromMatrix(counts, m, as.formula(form)); dds <- DESeq(dds, quiet=TRUE)
  r <- results(dds, name=tail(resultsNames(dds),1)); r <- r[!is.na(r$padj),]
  logf("\n=== %s : %s ===", label, form)
  logf("  FDR<0.05          : %d", sum(r$padj<0.05))
  logf("  FDR<0.05 |LFC|>=1 : %d", sum(r$padj<0.05 & abs(r$log2FoldChange)>=1))
  if(sum(r$padj<0.05)){ top <- r[order(r$padj),][1:min(6,sum(r$padj<0.05)),]; print(data.frame(id=rownames(top), log2FC=round(top$log2FoldChange,3), padj=signif(top$padj,3)), row.names=FALSE)}; if(FALSE){
    sym <- AnnotationDbi::mapIds(org.Hs.eg.db, sub("\\..*","",rownames(top)),"SYMBOL","ENSEMBL",multiVals="first")
    print(data.frame(gene=ifelse(is.na(sym),rownames(top),sym),
                     log2FC=round(top$log2FoldChange,3), padj=signif(top$padj,3)), row.names=FALSE) }
  write.csv(as.data.frame(r), file.path(OUT, sprintf("ASD_pergene_%s.csv", label))); invisible(r) }
run("~ series + sex + Astro_z + Neuro_z + diagnosis", "with_composition")
run("~ series + sex + diagnosis",                      "without_composition")
logf("\n=== COMPOSITION TESTED SEPARATELY ===")
comp <- do.call(rbind, lapply(c("Astro","Neuro"), function(v){
  fit <- lm(as.formula(sprintf("%s ~ series + sex + diagnosis", v)), data=m)
  s <- summary(fit)$coefficients; r <- s[grep("^diagnosis", rownames(s)), ]
  data.frame(cell_type=v, mean_control=mean(m[[v]][m$diagnosis=="Control"]),
             mean_ASD=mean(m[[v]][m$diagnosis=="ASD"]), beta=unname(r[1]), p=unname(r[4])) }))
print(format(comp, digits=3), row.names=FALSE)
write.csv(comp, file.path(OUT,"ASD_composition_tested_separately.csv"), row.names=FALSE)
logf("\nwritten to %s", OUT)
