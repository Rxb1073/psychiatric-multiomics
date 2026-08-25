#!/usr/bin/env Rscript
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(library(DESeq2))
logf<-function(...) cat(sprintf(...),"\n")
AUT<-"/rds/projects/e/elhamsak-dissertation-autism"; AUT2<-"/rds/projects/e/elhamsak-dissertation-autism-2"
STG<-file.path(AUT,"RNA_seq/_pheno_staging"); DEC<-file.path(AUT2,"RNA_seq/Deconvolution/cell_type_proportions_MB.csv")
cohorts<-list(
 GSE102741=list(c=file.path(AUT,"RNA_seq/GSE102741/counts/GSE102741_gene_counts.txt"),s=file.path(STG,"GSE102741_samplesheet.tsv")),
 GSE59288_GSE51264=list(c=file.path(AUT,"RNA_seq/GSE59288_GSE51264/counts/GSE59288_GSE51264_gene_counts.txt"),s=file.path(STG,"GSE59288_GSE51264_samplesheet.tsv")),
 PRJNA434002=list(c=file.path(AUT2,"RNA_seq/PRJNA434002_bulk/counts/PRJNA434002_150bp_gene_counts.txt"),s=file.path(STG,"PRJNA434002_samplesheet.tsv")))
read_fc<-function(p){f<-readLines(p,n=1);sk<-if(grepl("^# ",f))1 else 0;m<-read.delim(p,header=TRUE,skip=sk,check.names=FALSE,row.names=1)
 d<-intersect(c("Chr","Start","End","Strand","Length"),colnames(m));if(length(d))m<-m[,setdiff(colnames(m),d),drop=FALSE]
 colnames(m)<-regmatches(colnames(m),regexpr("[SED]RR[0-9]+",colnames(m)));as.matrix(m)}
ndx<-function(x){x<-tolower(trimws(x));ifelse(grepl("asd|autism",x),"ASD",ifelse(grepl("control|ctrl|con\\b",x),"Control","OTHER"))}
prop<-read.csv(DEC,stringsAsFactors=FALSE);prop$srr<-gsub('"','',prop$srr)
build<-function(use){
 mats<-list();metas<-list()
 for(nm in use){cm<-read_fc(cohorts[[nm]]$c);ss<-read.delim(cohorts[[nm]]$s,check.names=FALSE,stringsAsFactors=FALSE)
  idc<-grep("^srr$",names(ss),ignore.case=TRUE,value=TRUE)[1];dxc<-grep("diagnos",names(ss),ignore.case=TRUE,value=TRUE)[1]
  sxc<-grep("^sex$|gender",names(ss),ignore.case=TRUE,value=TRUE)[1];agc<-grep("^age$",names(ss),ignore.case=TRUE,value=TRUE)[1]
  ss$.id<-regmatches(ss[[idc]],regexpr("[SED]RR[0-9]+",ss[[idc]]))
  md<-data.frame(srr=ss$.id,series=nm,diagnosis=ndx(ss[[dxc]]),sex=tolower(substr(trimws(ss[[sxc]]),1,1)),
    age=suppressWarnings(as.numeric(as.character(ss[[agc]]))),stringsAsFactors=FALSE)
  common<-intersect(colnames(cm),md$srr);md<-md[match(common,md$srr),];metas[[nm]]<-md}
 meta<-do.call(rbind,metas);rownames(meta)<-meta$srr
 meta$Astro<-prop$Astrocytes[match(meta$srr,prop$srr)]
 ok<-meta$diagnosis%in%c("ASD","Control")&meta$sex%in%c("m","f")&!is.na(meta$age)&!is.na(meta$Astro)&meta$age<18
 m<-meta[ok,]
 pv<-function(mm) suppressWarnings(wilcox.test(age~diagnosis,mm)$p.value)
 p0<-pv(m); n0<-nrow(m)
 while(pv(m)<=0.05){medA<-median(m$age[m$diagnosis=="ASD"]);medC<-median(m$age[m$diagnosis=="Control"])
  if(sum(m$diagnosis=="ASD")<=12||sum(m$diagnosis=="Control")<=12)break
  if(medA>=medC){idx<-which(m$diagnosis=="ASD");d<-idx[which.max(m$age[idx])]}else{idx<-which(m$diagnosis=="Control");d<-idx[which.min(m$age[idx])]}
  m<-m[-d,,drop=FALSE]}
 list(m=m,p0=p0,n0=n0,p1=pv(m))}
report<-function(use,tag){
 o<-build(use);m<-o$m
 logf("\n=== %s ===",tag)
 logf("pre-match:  n=%d  Wilcoxon age p=%.4g",o$n0,o$p0)
 logf("post-match: n=%d (ASD=%d Ctrl=%d)  Wilcoxon age p=%.4g",nrow(m),sum(m$diagnosis=="ASD"),sum(m$diagnosis=="Control"),o$p1)
 logf("post-match age: ASD median=%.1f  Control median=%.1f",median(m$age[m$diagnosis=="ASD"]),median(m$age[m$diagnosis=="Control"]))}
report(names(cohorts),"POOLED 3-cohort")
report("GSE59288_GSE51264","GSE59288/51264 only")
report(c("GSE59288_GSE51264","PRJNA434002"),"2-cohort (original)")
