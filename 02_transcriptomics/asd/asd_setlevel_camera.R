## Set-level marker testing in the AUTISM arm, inter-gene correlation ESTIMATED.
## Loading block copied from asd_child_confoundfree.R, so the sample set matches
## the arm's primary analysis. Run with and without composition covariates,
## because autism is the only arm carrying them in the model.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(limma); library(edgeR); library(org.Hs.eg.db); library(AnnotationDbi)})
logf <- function(...) cat(sprintf(...),"\n")
AUT <- "/rds/projects/e/elhamsak-dissertation-autism"
AUT2<- "/rds/projects/e/elhamsak-dissertation-autism-2"
STG <- file.path(AUT,"RNA_seq/_pheno_staging")
DEC <- file.path(AUT2,"RNA_seq/Deconvolution/cell_type_proportions_MB.csv")
OR  <- file.path(AUT2,"RNA_seq/DE_analysis")
OUT <- file.path(OR,"DISSERTATION_final/Transcriptomics_ASD/tables/setlevel")
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
for(v in c("Astro","Neuro")) m[[v]] <- as.numeric(scale(m[[v]]))
bt <- read.delim(file.path(OR,"gene_id2biotype.tsv"), header=FALSE, col.names=c("gid","bt"))
KEEP <- bt$gid[bt$bt %in% c("protein_coding","lncRNA")]
counts <- counts[rownames(counts) %in% KEEP, , drop=FALSE]
counts <- counts[rowSums(counts>=10) >= min(table(m$diagnosis)), ]
logf("genes after biotype and low-count filtering: %d", nrow(counts))
sym <- AnnotationDbi::mapIds(org.Hs.eg.db, sub("\\..*","",rownames(counts)),
                             "SYMBOL","ENSEMBL", multiVals="first")
SETS <- list(
  GABAergic      = c("SST","PVALB","GAD1","GAD2","VIP","NPY","CCK","LHX6","RELN","CNR1"),
  IEG            = c("EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3","FOS","FOSB","JUN","JUNB",
                     "ARC","ATF3","DUSP1","DUSP6","NPAS4"),
  excitatory     = c("SLC17A7","SLC17A6","CAMK2A","SATB2","TBR1","NEUROD6"),
  astrocyte      = c("GFAP","ALDH1L1","AQP4","S100B","SLC1A2","SLC1A3","VIM"),
  oligodendrocyte= c("MBP","PLP1","MOG","MAG","CNP","SOX10"),
  CalbindinCR    = c("CALB1","CALB2"))
idx <- lapply(SETS, function(s) which(sym %in% s))
idx <- idx[sapply(idx, length) >= 2]
logf("marker sets testable: %s", paste(sprintf("%s(%d)", names(idx), sapply(idx,length)), collapse=" "))
y <- calcNormFactors(DGEList(counts))
models <- list(
  with_composition    = if(nlevels(m$series)>1) ~ series+sex+Astro+Neuro+diagnosis else ~ sex+Astro+Neuro+diagnosis,
  without_composition = if(nlevels(m$series)>1) ~ series+sex+diagnosis            else ~ sex+diagnosis)
res <- list()
for(nm in names(models)){
  des <- model.matrix(models[[nm]], m)
  stopifnot(tail(colnames(des),1) == "diagnosisASD")
  v   <- voom(y, des)
  cam <- camera(v, idx, des, contrast="diagnosisASD", inter.gene.cor=NA)
  logf("\n===== %s : %s =====", nm, paste(deparse(models[[nm]]), collapse=""))
  print(cam)
  write.csv(cam, file.path(OUT, sprintf("ASD_CAMERA_%s.csv", nm)))
  res[[nm]] <- data.frame(model=nm, set=rownames(cam),
                          cam[,c("NGenes","Correlation","Direction","PValue","FDR")])
}
S <- do.call(rbind, res); rownames(S) <- NULL
write.csv(S, file.path(OUT,"ASD_setlevel_summary.csv"), row.names=FALSE)
cat("\n===== IEG and GABAergic, all three arms =====\n")
cat("  SCZ (camera, igc estimated): IEG corr 0.501, FDR 0.112 | GABAergic FDR 0.0082\n")
cat("  MDD (camera, igc estimated): IEG corr 0.581, FDR 0.023\n")
print(S[S$set %in% c("IEG","GABAergic"), c("model","set","NGenes","Correlation","Direction","PValue","FDR")], row.names=FALSE)
logf("\nwritten to %s", OUT)
