#!/usr/bin/env Rscript
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(DESeq2); library(org.Hs.eg.db)})
logf<-function(...) cat(sprintf(...),"\n")
AUT<-"/rds/projects/e/elhamsak-dissertation-autism"; AUT2<-"/rds/projects/e/elhamsak-dissertation-autism-2"
STG<-file.path(AUT,"RNA_seq/_pheno_staging"); DEC<-file.path(AUT2,"RNA_seq/Deconvolution/cell_type_proportions_MB.csv")
OUT<-file.path(AUT2,"RNA_seq/DE_analysis/POOLED_cortex")
# RIN available only for these two:
cohorts<-list(
 GSE102741=list(c=file.path(AUT,"RNA_seq/GSE102741/counts/GSE102741_gene_counts.txt"),s=file.path(STG,"GSE102741_samplesheet.tsv")),
 GSE59288_GSE51264=list(c=file.path(AUT,"RNA_seq/GSE59288_GSE51264/counts/GSE59288_GSE51264_gene_counts.txt"),s=file.path(STG,"GSE59288_GSE51264_samplesheet.tsv")))
read_fc<-function(p){f<-readLines(p,n=1);sk<-if(grepl("^# ",f))1 else 0;m<-read.delim(p,header=TRUE,skip=sk,check.names=FALSE,row.names=1)
 d<-intersect(c("Chr","Start","End","Strand","Length"),colnames(m));if(length(d))m<-m[,setdiff(colnames(m),d),drop=FALSE]
 colnames(m)<-regmatches(colnames(m),regexpr("[SED]RR[0-9]+",colnames(m)));as.matrix(m)}
ndx<-function(x){x<-tolower(trimws(x));ifelse(grepl("asd|autism",x),"ASD",ifelse(grepl("control|ctrl|con\\b",x),"Control","OTHER"))}
prop<-read.csv(DEC,stringsAsFactors=FALSE);prop$srr<-gsub('"','',prop$srr)

mats<-list();metas<-list()
for(nm in names(cohorts)){
 cm<-read_fc(cohorts[[nm]]$c);ss<-read.delim(cohorts[[nm]]$s,check.names=FALSE,stringsAsFactors=FALSE)
 idc<-grep("^srr$",names(ss),ignore.case=TRUE,value=TRUE)[1];dxc<-grep("diagnos",names(ss),ignore.case=TRUE,value=TRUE)[1]
 sxc<-grep("^sex$|gender",names(ss),ignore.case=TRUE,value=TRUE)[1];agc<-grep("^age$",names(ss),ignore.case=TRUE,value=TRUE)[1]
 rnc<-grep("^rin$|rin",names(ss),ignore.case=TRUE,value=TRUE)[1]
 logf("%-18s RIN column detected: %s",nm,ifelse(is.na(rnc),"NONE",rnc))
 ss$.id<-regmatches(ss[[idc]],regexpr("[SED]RR[0-9]+",ss[[idc]]))
 md<-data.frame(srr=ss$.id,series=nm,diagnosis=ndx(ss[[dxc]]),sex=tolower(substr(trimws(ss[[sxc]]),1,1)),
   age=suppressWarnings(as.numeric(as.character(ss[[agc]]))),
   rin=if(is.na(rnc)) NA_real_ else suppressWarnings(as.numeric(as.character(ss[[rnc]]))),stringsAsFactors=FALSE)
 common<-intersect(colnames(cm),md$srr);cm<-cm[,common];md<-md[match(common,md$srr),];mats[[nm]]<-cm;metas[[nm]]<-md}
genes<-Reduce(intersect,lapply(mats,rownames));counts<-do.call(cbind,lapply(mats,function(x)x[genes,,drop=FALSE]))
meta<-do.call(rbind,metas);rownames(meta)<-meta$srr;meta<-meta[colnames(counts),]
meta$Astro<-prop$Astrocytes[match(meta$srr,prop$srr)];meta$Neuro<-prop$Neurons[match(meta$srr,prop$srr)]

# child + require RIN present
ok<-meta$diagnosis%in%c("ASD","Control")&meta$sex%in%c("m","f")&!is.na(meta$age)&!is.na(meta$Astro)&!is.na(meta$rin)&meta$age<18
counts<-round(counts[,ok]);m<-meta[ok,]
cat("\n===== RIN by group (is RIN confounded with diagnosis?) =====\n")
for(gp in c("ASD","Control")){x<-m$rin[m$diagnosis==gp];logf("%-8s n=%d  RIN median=%.2f mean=%.2f range %.1f-%.1f",gp,length(x),median(x),mean(x),min(x),max(x))}
logf("Wilcoxon RIN ASD vs Control: p=%.3g",suppressWarnings(wilcox.test(rin~diagnosis,m)$p.value))

# age-match (same rule)
pv<-function(mm) suppressWarnings(wilcox.test(age~diagnosis,mm)$p.value)
while(pv(m)<=0.05){medA<-median(m$age[m$diagnosis=="ASD"]);medC<-median(m$age[m$diagnosis=="Control"])
 if(sum(m$diagnosis=="ASD")<=12||sum(m$diagnosis=="Control")<=12)break
 if(medA>=medC){idx<-which(m$diagnosis=="ASD");d<-idx[which.max(m$age[idx])]}else{idx<-which(m$diagnosis=="Control");d<-idx[which.min(m$age[idx])]}
 m<-m[-d,,drop=FALSE]}
counts<-counts[,rownames(m)]
m$diagnosis<-factor(as.character(m$diagnosis),levels=c("Control","ASD"));m$series<-droplevels(factor(m$series));m$sex<-droplevels(factor(m$sex))
for(v in c("Astro","Neuro","rin"))m[[v]]<-as.numeric(scale(m[[v]]))
counts<-counts[rowSums(counts>=10)>=min(table(m$diagnosis)),]
logf("\nmatched pool with RIN: n=%d (ASD=%d Ctrl=%d) genes=%d",nrow(m),sum(m$diagnosis=="ASD"),sum(m$diagnosis=="Control"),nrow(counts))

run<-function(des,tag){
 dds<-DESeq(DESeqDataSetFromMatrix(counts,m,des),quiet=TRUE)
 r<-as.data.frame(results(dds,contrast=c("diagnosis","ASD","Control"),alpha=0.05));r$ens<-sub("\\..*","",rownames(r))
 logf("[%s] %-45s  FDR<0.05=%d | FDR<0.05&|LFC|>=1=%d",tag,deparse(des),sum(r$padj<0.05,na.rm=TRUE),sum(r$padj<0.05&abs(r$log2FoldChange)>=1,na.rm=TRUE))
 r}
sr<-if(nlevels(m$series)>1)"series+" else ""
cat("\n===== MODEL COMPARISON: without vs with RIN =====\n")
r_no <-run(as.formula(paste0("~ ",sr,"sex+Astro+Neuro+diagnosis")),"noRIN")
r_rin<-run(as.formula(paste0("~ ",sr,"sex+rin+Astro+Neuro+diagnosis")),"withRIN")

key<-c("DLX1","DLX2","ARX","LHX6","SYNGAP1","DCX","SOX11","DRAXIN","DPYSL3","ST8SIA2")
e<-mapIds(org.Hs.eg.db,key,"ENSEMBL","SYMBOL")
cat("\n===== KEY GENES: padj without vs with RIN =====\n")
cat(sprintf("%-9s %12s %12s\n","gene","padj_noRIN","padj_RIN"))
for(g in key){i1<-match(e[g],r_no$ens);i2<-match(e[g],r_rin$ens)
 cat(sprintf("%-9s %12.3g %12.3g\n",g,if(is.na(i1))NA else r_no$padj[i1],if(is.na(i2))NA else r_rin$padj[i2]))}
