#!/usr/bin/env Rscript
# Pooled 3-cohort childhood ASD DE, age-balanced + composition-adjusted.
# Produces the n=89 primary pooled result (985 FDR<0.05; 53 with |LFC|>=1).
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(DESeq2); library(org.Hs.eg.db)})
QC_FAIL<-c("SRR5938468","SRR5938447","SRR5938445","SRR1510129","SRR5938444","SRR5938441","SRR1004284","SRR5938466")
logf<-function(...) cat(sprintf(...),"\n")
AUT<-"/rds/projects/e/elhamsak-dissertation-autism"; AUT2<-"/rds/projects/e/elhamsak-dissertation-autism-2"
STG<-file.path(AUT,"RNA_seq/_pheno_staging"); DEC<-file.path(AUT2,"RNA_seq/Deconvolution/cell_type_proportions_MB.csv")
OUT<-file.path(AUT2,"RNA_seq/DE_analysis/POOLED_cortex")
cohorts<-list(
 GSE102741=list(c=file.path(AUT,"RNA_seq/GSE102741/counts/GSE102741_gene_counts.txt"),s=file.path(STG,"GSE102741_samplesheet.tsv")),
 GSE59288_GSE51264=list(c=file.path(AUT,"RNA_seq/GSE59288_GSE51264/counts/GSE59288_GSE51264_gene_counts.txt"),s=file.path(STG,"GSE59288_GSE51264_samplesheet.tsv")),
 PRJNA434002=list(c=file.path(AUT2,"RNA_seq/PRJNA434002_bulk/counts/PRJNA434002_150bp_gene_counts.txt"),s=file.path(STG,"PRJNA434002_samplesheet.tsv")))
for(nm in names(cohorts)) for(k in c("c","s")) if(!file.exists(cohorts[[nm]][[k]])) stop(sprintf("MISSING %s (%s)",nm,k))
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
 ss$.id<-regmatches(ss[[idc]],regexpr("[SED]RR[0-9]+",ss[[idc]]))
 md<-data.frame(srr=ss$.id,series=nm,diagnosis=ndx(ss[[dxc]]),sex=tolower(substr(trimws(ss[[sxc]]),1,1)),
   age=suppressWarnings(as.numeric(as.character(ss[[agc]]))),stringsAsFactors=FALSE)
 common<-intersect(colnames(cm),md$srr);cm<-cm[,common];md<-md[match(common,md$srr),];mats[[nm]]<-cm;metas[[nm]]<-md}
genes<-Reduce(intersect,lapply(mats,rownames));counts<-do.call(cbind,lapply(mats,function(x)x[genes,,drop=FALSE]))
meta<-do.call(rbind,metas);rownames(meta)<-meta$srr;meta<-meta[colnames(counts),]
meta$Astro<-prop$Astrocytes[match(meta$srr,prop$srr)];meta$Neuro<-prop$Neurons[match(meta$srr,prop$srr)]
ok<-meta$diagnosis%in%c("ASD","Control")&meta$sex%in%c("m","f")&!is.na(meta$age)&!is.na(meta$Astro)&meta$age<18&!(meta$srr%in%QC_FAIL)
counts<-round(counts[,ok]);m<-meta[ok,]
cat("===== PER-COHORT n AFTER CHILD FILTER =====\n");print(table(series=m$series,diagnosis=m$diagnosis))
pv<-function(mm) suppressWarnings(wilcox.test(age~diagnosis,mm)$p.value)
logf("pre-match age p=%.4g (n=%d)",pv(m),nrow(m))
while(pv(m)<=0.05){medA<-median(m$age[m$diagnosis=="ASD"]);medC<-median(m$age[m$diagnosis=="Control"])
 if(sum(m$diagnosis=="ASD")<=12||sum(m$diagnosis=="Control")<=12){logf("floor reached");break}
 if(medA>=medC){i<-which(m$diagnosis=="ASD");d<-i[which.max(m$age[i])]}else{i<-which(m$diagnosis=="Control");d<-i[which.min(m$age[i])]}
 m<-m[-d,,drop=FALSE]}
counts<-counts[,rownames(m)]
logf("post-match n=%d (ASD=%d Ctrl=%d) age p=%.4g",nrow(m),sum(m$diagnosis=="ASD"),sum(m$diagnosis=="Control"),pv(m))
m$diagnosis<-factor(as.character(m$diagnosis),levels=c("Control","ASD"))
m$series<-droplevels(factor(m$series));m$sex<-droplevels(factor(m$sex))
for(v in c("Astro","Neuro"))m[[v]]<-as.numeric(scale(m[[v]]))
counts<-counts[rowSums(counts>=10)>=min(table(m$diagnosis)),]
des<-if(nlevels(m$series)>1) ~ series+sex+Astro+Neuro+diagnosis else ~ sex+Astro+Neuro+diagnosis
logf("design: %s | genes=%d",deparse(des),nrow(counts))
dds<-DESeq(DESeqDataSetFromMatrix(counts,m,des),quiet=TRUE)
r<-as.data.frame(results(dds,contrast=c("diagnosis","ASD","Control"),alpha=0.05))
r$ens<-sub("\\..*","",rownames(r));r$sym<-mapIds(org.Hs.eg.db,r$ens,"SYMBOL","ENSEMBL");r$sym[is.na(r$sym)]<-r$ens[is.na(r$sym)]
logf("FDR<0.05 = %d | FDR<0.05 & |LFC|>=1 = %d",sum(r$padj<0.05,na.rm=TRUE),sum(r$padj<0.05&abs(r$log2FoldChange)>=1,na.rm=TRUE))
write.csv(r[order(r$padj),],file.path(OUT,"DE_ASD_child_3dataset_full.csv"),row.names=FALSE)
st<-r[which(r$padj<0.05&abs(r$log2FoldChange)>=1),];write.csv(st[order(st$padj),],file.path(OUT,"DE_ASD_child_3dataset_strict.csv"),row.names=FALSE)
key<-c("DLX1","DLX2","ARX","LHX6","SYNGAP1","DCX","SOX11","DRAXIN","DPYSL3","ST8SIA2","PCDH19","PCDHGC3")
e<-mapIds(org.Hs.eg.db,key,"ENSEMBL","SYMBOL")
cat("\n===== KEY GENES =====\n")
for(g in key){i<-match(e[g],r$ens);if(!is.na(i))cat(sprintf("%-9s LFC=%6.2f padj=%9.3g\n",g,r$log2FoldChange[i],r$padj[i]))}
