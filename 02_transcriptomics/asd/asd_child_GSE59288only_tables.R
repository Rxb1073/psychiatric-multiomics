#!/usr/bin/env Rscript
# Single-cohort (GSE59288/GSE51264) childhood ASD DE, age-matched + composition-adjusted.
# Reproduces the n=36 result (352 FDR<0.05; 186 with |LFC|>=1).
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(DESeq2); library(org.Hs.eg.db)})
logf<-function(...) cat(sprintf(...),"\n")
AUT<-"/rds/projects/e/elhamsak-dissertation-autism"; AUT2<-"/rds/projects/e/elhamsak-dissertation-autism-2"
STG<-file.path(AUT,"RNA_seq/_pheno_staging"); DEC<-file.path(AUT2,"RNA_seq/Deconvolution/cell_type_proportions_MB.csv")
OUT<-file.path(AUT2,"RNA_seq/DE_analysis/POOLED_cortex")
CNT<-file.path(AUT,"RNA_seq/GSE59288_GSE51264/counts/GSE59288_GSE51264_gene_counts.txt")
SHT<-file.path(STG,"GSE59288_GSE51264_samplesheet.tsv")
for(p in c(CNT,SHT,DEC)) if(!file.exists(p)) stop(sprintf("MISSING: %s",p))
read_fc<-function(p){f<-readLines(p,n=1);sk<-if(grepl("^# ",f))1 else 0;m<-read.delim(p,header=TRUE,skip=sk,check.names=FALSE,row.names=1)
 d<-intersect(c("Chr","Start","End","Strand","Length"),colnames(m));if(length(d))m<-m[,setdiff(colnames(m),d),drop=FALSE]
 colnames(m)<-regmatches(colnames(m),regexpr("[SED]RR[0-9]+",colnames(m)));as.matrix(m)}
ndx<-function(x){x<-tolower(trimws(x));ifelse(grepl("asd|autism",x),"ASD",ifelse(grepl("control|ctrl|con\\b",x),"Control","OTHER"))}
prop<-read.csv(DEC,stringsAsFactors=FALSE);prop$srr<-gsub('"','',prop$srr)
cm<-read_fc(CNT);ss<-read.delim(SHT,check.names=FALSE,stringsAsFactors=FALSE)
idc<-grep("^srr$",names(ss),ignore.case=TRUE,value=TRUE)[1];dxc<-grep("diagnos",names(ss),ignore.case=TRUE,value=TRUE)[1]
sxc<-grep("^sex$|gender",names(ss),ignore.case=TRUE,value=TRUE)[1];agc<-grep("^age$",names(ss),ignore.case=TRUE,value=TRUE)[1]
ss$.id<-regmatches(ss[[idc]],regexpr("[SED]RR[0-9]+",ss[[idc]]))
md<-data.frame(srr=ss$.id,diagnosis=ndx(ss[[dxc]]),sex=tolower(substr(trimws(ss[[sxc]]),1,1)),
  age=suppressWarnings(as.numeric(as.character(ss[[agc]]))),stringsAsFactors=FALSE)
common<-intersect(colnames(cm),md$srr);cm<-cm[,common];md<-md[match(common,md$srr),]
md$Astro<-prop$Astrocytes[match(md$srr,prop$srr)];md$Neuro<-prop$Neurons[match(md$srr,prop$srr)];rownames(md)<-md$srr
ok<-md$diagnosis%in%c("ASD","Control")&md$sex%in%c("m","f")&!is.na(md$age)&!is.na(md$Astro)&md$age<18
counts<-round(cm[,ok]);m<-md[ok,]
pv<-function(mm) suppressWarnings(wilcox.test(age~diagnosis,mm)$p.value)
logf("pre-match n=%d age p=%.4g",nrow(m),pv(m))
while(pv(m)<=0.05){medA<-median(m$age[m$diagnosis=="ASD"]);medC<-median(m$age[m$diagnosis=="Control"])
 if(sum(m$diagnosis=="ASD")<=12||sum(m$diagnosis=="Control")<=12){logf("floor reached");break}
 if(medA>=medC){i<-which(m$diagnosis=="ASD");d<-i[which.max(m$age[i])]}else{i<-which(m$diagnosis=="Control");d<-i[which.min(m$age[i])]}
 m<-m[-d,,drop=FALSE]}
counts<-counts[,rownames(m)]
logf("post-match n=%d (ASD=%d Ctrl=%d) age p=%.4g | med ASD=%.1f Ctrl=%.1f",nrow(m),sum(m$diagnosis=="ASD"),
  sum(m$diagnosis=="Control"),pv(m),median(m$age[m$diagnosis=="ASD"]),median(m$age[m$diagnosis=="Control"]))
m$diagnosis<-factor(as.character(m$diagnosis),levels=c("Control","ASD"));m$sex<-droplevels(factor(m$sex))
for(v in c("Astro","Neuro"))m[[v]]<-as.numeric(scale(m[[v]]))
counts<-counts[rowSums(counts>=10)>=min(table(m$diagnosis)),]
des<-if(nlevels(m$sex)>1) ~ sex+Astro+Neuro+diagnosis else ~ Astro+Neuro+diagnosis
logf("design: %s | genes=%d",deparse(des),nrow(counts))
dds<-DESeq(DESeqDataSetFromMatrix(counts,m,des),quiet=TRUE)
r<-as.data.frame(results(dds,contrast=c("diagnosis","ASD","Control"),alpha=0.05))
r$ens<-sub("\\..*","",rownames(r));r$sym<-mapIds(org.Hs.eg.db,r$ens,"SYMBOL","ENSEMBL");r$sym[is.na(r$sym)]<-r$ens[is.na(r$sym)]
logf("FDR<0.05 = %d | FDR<0.05 & |LFC|>=1 = %d",sum(r$padj<0.05,na.rm=TRUE),sum(r$padj<0.05&abs(r$log2FoldChange)>=1,na.rm=TRUE))
write.csv(r[order(r$padj),],file.path(OUT,"DE_ASD_child_GSE59288only_full.csv"),row.names=FALSE)
st<-r[which(r$padj<0.05&abs(r$log2FoldChange)>=1),];write.csv(st[order(st$padj),],file.path(OUT,"DE_ASD_child_GSE59288only_strict.csv"),row.names=FALSE)
key<-c("DLX1","DLX2","ARX","LHX6","SYNGAP1","DCX","SOX11","DRAXIN","DPYSL3","ST8SIA2","PCDH19","PCDHGC3")
e<-mapIds(org.Hs.eg.db,key,"ENSEMBL","SYMBOL")
cat("\n===== KEY GENES =====\n")
for(g in key){i<-match(e[g],r$ens);if(!is.na(i))cat(sprintf("%-9s LFC=%6.2f padj=%9.3g\n",g,r$log2FoldChange[i],r$padj[i]))}
