#!/usr/bin/env Rscript
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(DESeq2); library(org.Hs.eg.db)})
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

run_cohort<-function(nm){
 cm<-read_fc(cohorts[[nm]]$c);ss<-read.delim(cohorts[[nm]]$s,check.names=FALSE,stringsAsFactors=FALSE)
 idc<-grep("^srr$",names(ss),ignore.case=TRUE,value=TRUE)[1];dxc<-grep("diagnos",names(ss),ignore.case=TRUE,value=TRUE)[1]
 sxc<-grep("^sex$|gender",names(ss),ignore.case=TRUE,value=TRUE)[1];agc<-grep("^age$",names(ss),ignore.case=TRUE,value=TRUE)[1]
 ss$.id<-regmatches(ss[[idc]],regexpr("[SED]RR[0-9]+",ss[[idc]]))
 md<-data.frame(srr=ss$.id,diagnosis=ndx(ss[[dxc]]),sex=tolower(substr(trimws(ss[[sxc]]),1,1)),
   age=suppressWarnings(as.numeric(as.character(ss[[agc]]))),stringsAsFactors=FALSE)
 common<-intersect(colnames(cm),md$srr);cm<-cm[,common];md<-md[match(common,md$srr),]
 md$Astro<-prop$Astrocytes[match(md$srr,prop$srr)];md$Neuro<-prop$Neurons[match(md$srr,prop$srr)];rownames(md)<-md$srr
 ok<-md$diagnosis%in%c("ASD","Control")&md$sex%in%c("m","f")&!is.na(md$age)&!is.na(md$Astro)&md$age<18
 cm<-round(cm[,ok]);md<-md[ok,]
 logf("%-18s child n=%d  (ASD=%d Ctrl=%d)",nm,nrow(md),sum(md$diagnosis=="ASD"),sum(md$diagnosis=="Control"))
 if(sum(md$diagnosis=="ASD")<3||sum(md$diagnosis=="Control")<3){logf("   -> too few in one group; skipping DE for %s",nm);return(NULL)}
 md$diagnosis<-factor(as.character(md$diagnosis),levels=c("Control","ASD"));md$sex<-droplevels(factor(md$sex))
 for(v in c("Astro","Neuro"))md[[v]]<-as.numeric(scale(md[[v]]))
 cm<-cm[rowSums(cm>=10)>=min(table(md$diagnosis)),]
 des<-if(nlevels(md$sex)>1) ~ sex+Astro+Neuro+diagnosis else ~ Astro+Neuro+diagnosis
 dds<-DESeq(DESeqDataSetFromMatrix(cm,md,des),quiet=TRUE)
 r<-as.data.frame(results(dds,contrast=c("diagnosis","ASD","Control"),alpha=0.05));r$ens<-sub("\\..*","",rownames(r))
 logf("   composition-adjusted significant padj<0.05 = %d",sum(r$padj<0.05,na.rm=TRUE))
 r}

cat("===== PER-COHORT child composition-adjusted DE =====\n")
res<-lapply(names(cohorts),run_cohort);names(res)<-names(cohorts)

key<-c("DLX1","DLX2","ARX","LHX6","SYNGAP1","PCDH19","DCX","PCDHGC3","SOX11","DRAXIN","DPYSL3","ST8SIA2")
e<-mapIds(org.Hs.eg.db,key,"ENSEMBL","SYMBOL")
cat("\n===== KEY GENES: LFC (padj) per cohort =====\n")
cat(sprintf("%-9s | %-20s | %-20s | %-20s\n","gene","GSE102741","GSE59288_51264","PRJNA434002"))
fmt<-function(r, e){if(is.null(r))return("   --      ");i<-match(e,r$ens);if(is.na(i))return("   NA      ");sprintf("%6.2f (%7.3g)",r$log2FoldChange[i],r$padj[i])}
for(g in key) cat(sprintf("%-9s | %-20s | %-20s | %-20s\n",g,fmt(res[[1]],e[g]),fmt(res[[2]],e[g]),fmt(res[[3]],e[g])))

# direction concordance across the cohorts that ran
cat("\n===== DIRECTION across cohorts (key genes) =====\n")
ok_res<-res[!sapply(res,is.null)]
if(length(ok_res)>=2){
 signs<-sapply(ok_res,function(r){sapply(e,function(x){i<-match(x,r$ens);if(is.na(i))NA else sign(r$log2FoldChange[i])})})
 allneg<-apply(signs,1,function(s){s<-s[!is.na(s)];length(s)>0&&all(s<0)})
 logf("key genes downregulated in ALL cohorts that ran: %d/%d",sum(allneg,na.rm=TRUE),length(allneg))
}
