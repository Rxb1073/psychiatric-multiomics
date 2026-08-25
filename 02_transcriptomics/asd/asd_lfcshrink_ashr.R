#!/usr/bin/env Rscript
# lfcShrink (ashr) applied to all three ASD childhood analyses.
# Recomputes each dds, then writes shrunken results alongside the MLE versions.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(DESeq2); library(org.Hs.eg.db)})
logf<-function(...) cat(sprintf(...),"\n")
AUT<-"/rds/projects/e/elhamsak-dissertation-autism"; AUT2<-"/rds/projects/e/elhamsak-dissertation-autism-2"
STG<-file.path(AUT,"RNA_seq/_pheno_staging"); DEC<-file.path(AUT2,"RNA_seq/Deconvolution/cell_type_proportions_MB.csv")
OUT<-file.path(AUT2,"RNA_seq/DE_analysis/POOLED_cortex")
FDR<-0.05; LFC<-1
QC_FAIL<-c("SRR5938468","SRR5938447","SRR5938445","SRR1510129","SRR5938444","SRR5938441","SRR1004284","SRR5938466")

allc<-list(
 GSE102741=list(c=file.path(AUT,"RNA_seq/GSE102741/counts/GSE102741_gene_counts.txt"),s=file.path(STG,"GSE102741_samplesheet.tsv")),
 GSE59288_GSE51264=list(c=file.path(AUT,"RNA_seq/GSE59288_GSE51264/counts/GSE59288_GSE51264_gene_counts.txt"),s=file.path(STG,"GSE59288_GSE51264_samplesheet.tsv")),
 PRJNA434002=list(c=file.path(AUT2,"RNA_seq/PRJNA434002_bulk/counts/PRJNA434002_150bp_gene_counts.txt"),s=file.path(STG,"PRJNA434002_samplesheet.tsv")))
read_fc<-function(p){f<-readLines(p,n=1);sk<-if(grepl("^# ",f))1 else 0;m<-read.delim(p,header=TRUE,skip=sk,check.names=FALSE,row.names=1)
 d<-intersect(c("Chr","Start","End","Strand","Length"),colnames(m));if(length(d))m<-m[,setdiff(colnames(m),d),drop=FALSE]
 colnames(m)<-regmatches(colnames(m),regexpr("[SED]RR[0-9]+",colnames(m)));as.matrix(m)}
ndx<-function(x){x<-tolower(trimws(x));ifelse(grepl("asd|autism",x),"ASD",ifelse(grepl("control|ctrl|con\\b",x),"Control","OTHER"))}
prop<-read.csv(DEC,stringsAsFactors=FALSE);prop$srr<-gsub('"','',prop$srr)

build<-function(use,qc){
 mats<-list();metas<-list()
 for(nm in use){cm<-read_fc(allc[[nm]]$c);ss<-read.delim(allc[[nm]]$s,check.names=FALSE,stringsAsFactors=FALSE)
  idc<-grep("^srr$",names(ss),ignore.case=TRUE,value=TRUE)[1];dxc<-grep("diagnos",names(ss),ignore.case=TRUE,value=TRUE)[1]
  sxc<-grep("^sex$|gender",names(ss),ignore.case=TRUE,value=TRUE)[1];agc<-grep("^age$",names(ss),ignore.case=TRUE,value=TRUE)[1]
  ss$.id<-regmatches(ss[[idc]],regexpr("[SED]RR[0-9]+",ss[[idc]]))
  md<-data.frame(srr=ss$.id,series=nm,diagnosis=ndx(ss[[dxc]]),sex=tolower(substr(trimws(ss[[sxc]]),1,1)),
    age=suppressWarnings(as.numeric(as.character(ss[[agc]]))),stringsAsFactors=FALSE)
  common<-intersect(colnames(cm),md$srr);cm<-cm[,common];md<-md[match(common,md$srr),];mats[[nm]]<-cm;metas[[nm]]<-md}
 genes<-Reduce(intersect,lapply(mats,rownames));counts<-do.call(cbind,lapply(mats,function(x)x[genes,,drop=FALSE]))
 meta<-do.call(rbind,metas);rownames(meta)<-meta$srr;meta<-meta[colnames(counts),]
 meta$Astro<-prop$Astrocytes[match(meta$srr,prop$srr)];meta$Neuro<-prop$Neurons[match(meta$srr,prop$srr)]
 qcv<-if(qc) QC_FAIL else character(0)
 ok<-meta$diagnosis%in%c("ASD","Control")&meta$sex%in%c("m","f")&!is.na(meta$age)&!is.na(meta$Astro)&
     meta$age<18&!(meta$srr%in%qcv)
 counts<-round(counts[,ok]);m<-meta[ok,]
 pv<-function(mm) suppressWarnings(wilcox.test(age~diagnosis,mm)$p.value)
 flr<-if(length(use)==1) 12 else 8
 while(pv(m)<=0.05){medA<-median(m$age[m$diagnosis=="ASD"]);medC<-median(m$age[m$diagnosis=="Control"])
  if(sum(m$diagnosis=="ASD")<=flr||sum(m$diagnosis=="Control")<=flr)break
  if(medA>=medC){i<-which(m$diagnosis=="ASD");d<-i[which.max(m$age[i])]}else{i<-which(m$diagnosis=="Control");d<-i[which.min(m$age[i])]}
  m<-m[-d,,drop=FALSE]}
 counts<-counts[,rownames(m)]
 m$diagnosis<-factor(as.character(m$diagnosis),levels=c("Control","ASD"))
 m$series<-droplevels(factor(m$series));m$sex<-droplevels(factor(m$sex))
 for(v in c("Astro","Neuro"))m[[v]]<-as.numeric(scale(m[[v]]))
 counts<-counts[rowSums(counts>=10)>=min(table(m$diagnosis)),]
 des<-if(nlevels(m$series)>1) ~ series+sex+Astro+Neuro+diagnosis else ~ sex+Astro+Neuro+diagnosis
 list(dds=DESeq(DESeqDataSetFromMatrix(counts,m,des),quiet=TRUE),n=nrow(m))}

run<-function(use,tag,qc){
 o<-build(use,qc); dds<-o$dds
 cn<-resultsNames(dds); coef<-grep("diagnosis_ASD_vs_Control",cn,value=TRUE)[1]
 if(is.na(coef)) stop("diagnosis coefficient not found: ",paste(cn,collapse=", "))
 mle<-as.data.frame(results(dds,contrast=c("diagnosis","ASD","Control"),alpha=FDR))
 shr<-as.data.frame(lfcShrink(dds,contrast=c("diagnosis","ASD","Control"),type="ashr"))
 d<-data.frame(ens=sub("\\..*","",rownames(mle)),
               lfc_MLE=mle$log2FoldChange, lfc_shrunk=shr$log2FoldChange,
               padj=mle$padj, baseMean=mle$baseMean, stringsAsFactors=FALSE)
 d$sym<-mapIds(org.Hs.eg.db,d$ens,"SYMBOL","ENSEMBL"); d$sym[is.na(d$sym)]<-d$ens[is.na(d$sym)]
 sM<-sum(d$padj<FDR & abs(d$lfc_MLE)>=LFC,na.rm=TRUE)
 sS<-sum(d$padj<FDR & abs(d$lfc_shrunk)>=LFC,na.rm=TRUE)
 logf("[%-16s] n=%3d | FDR<0.05=%4d | strict MLE=%4d | strict SHRUNK=%4d",
      tag,o$n,sum(d$padj<FDR,na.rm=TRUE),sM,sS)
 write.csv(d[order(d$padj),],file.path(OUT,sprintf("DE_%s_lfcShrink.csv",tag)),row.names=FALSE)
 key<-c("DLX1","DLX2","ARX","LHX6","SYNGAP1","DCX","SOX11","DRAXIN","DPYSL3","ST8SIA2")
 k<-d[d$sym%in%key,c("sym","lfc_MLE","lfc_shrunk","padj")]
 if(nrow(k)) {cat("   key genes (MLE -> shrunk):\n"); print(k[order(k$padj),],row.names=FALSE)}
 invisible(d)}

cat("===== lfcShrink (ashr) across ASD childhood analyses =====\n")
run(names(allc),"ASD_child_3dataset",FALSE)
run("GSE59288_GSE51264","ASD_child_GSE59288only",FALSE)
run(c("GSE102741","PRJNA434002"),"ASD_child_confoundfree",TRUE)
logf("\nOutputs -> %s",OUT)
