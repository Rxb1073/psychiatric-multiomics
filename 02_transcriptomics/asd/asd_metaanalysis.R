#!/usr/bin/env Rscript
# META-ANALYSIS of the two confound-free cohorts (GSE102741, PRJNA434002).
# Per-cohort DESeq2, then inverse-variance fixed- and random-effects
# combination of log2FC using lfcSE. Reports I^2 heterogeneity.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(DESeq2); library(org.Hs.eg.db)})
logf<-function(...) cat(sprintf(...),"\n")
AUT<-"/rds/projects/e/elhamsak-dissertation-autism"; AUT2<-"/rds/projects/e/elhamsak-dissertation-autism-2"
STG<-file.path(AUT,"RNA_seq/_pheno_staging"); DEC<-file.path(AUT2,"RNA_seq/Deconvolution/cell_type_proportions_MB.csv")
OUT<-file.path(AUT2,"RNA_seq/DE_analysis/POOLED_cortex")
QC_FAIL<-c("SRR5938468","SRR5938447","SRR5938445","SRR1510129","SRR5938444","SRR5938441","SRR1004284","SRR5938466")

cohorts<-list(
 GSE102741=list(c=file.path(AUT,"RNA_seq/GSE102741/counts/GSE102741_gene_counts.txt"),s=file.path(STG,"GSE102741_samplesheet.tsv")),
 PRJNA434002=list(c=file.path(AUT2,"RNA_seq/PRJNA434002_bulk/counts/PRJNA434002_150bp_gene_counts.txt"),s=file.path(STG,"PRJNA434002_samplesheet.tsv")))
read_fc<-function(p){f<-readLines(p,n=1);sk<-if(grepl("^# ",f))1 else 0;m<-read.delim(p,header=TRUE,skip=sk,check.names=FALSE,row.names=1)
 d<-intersect(c("Chr","Start","End","Strand","Length"),colnames(m));if(length(d))m<-m[,setdiff(colnames(m),d),drop=FALSE]
 colnames(m)<-regmatches(colnames(m),regexpr("[SED]RR[0-9]+",colnames(m)));as.matrix(m)}
ndx<-function(x){x<-tolower(trimws(x));ifelse(grepl("asd|autism",x),"ASD",ifelse(grepl("control|ctrl|con\\b",x),"Control","OTHER"))}
prop<-read.csv(DEC,stringsAsFactors=FALSE);prop$srr<-gsub('"','',prop$srr)

per_cohort<-function(nm){
 cm<-read_fc(cohorts[[nm]]$c);ss<-read.delim(cohorts[[nm]]$s,check.names=FALSE,stringsAsFactors=FALSE)
 idc<-grep("^srr$",names(ss),ignore.case=TRUE,value=TRUE)[1];dxc<-grep("diagnos",names(ss),ignore.case=TRUE,value=TRUE)[1]
 sxc<-grep("^sex$|gender",names(ss),ignore.case=TRUE,value=TRUE)[1];agc<-grep("^age$",names(ss),ignore.case=TRUE,value=TRUE)[1]
 ss$.id<-regmatches(ss[[idc]],regexpr("[SED]RR[0-9]+",ss[[idc]]))
 md<-data.frame(srr=ss$.id,diagnosis=ndx(ss[[dxc]]),sex=tolower(substr(trimws(ss[[sxc]]),1,1)),
   age=suppressWarnings(as.numeric(as.character(ss[[agc]]))),stringsAsFactors=FALSE)
 common<-intersect(colnames(cm),md$srr);cm<-cm[,common];md<-md[match(common,md$srr),]
 md$Astro<-prop$Astrocytes[match(md$srr,prop$srr)];md$Neuro<-prop$Neurons[match(md$srr,prop$srr)];rownames(md)<-md$srr
 ok<-md$diagnosis%in%c("ASD","Control")&md$sex%in%c("m","f")&!is.na(md$age)&!is.na(md$Astro)&
     md$age<18&!(md$srr%in%QC_FAIL)
 cm<-round(cm[,ok]);md<-md[ok,]
 md$diagnosis<-factor(as.character(md$diagnosis),levels=c("Control","ASD"));md$sex<-droplevels(factor(md$sex))
 for(v in c("Astro","Neuro"))md[[v]]<-as.numeric(scale(md[[v]]))
 logf("%-14s n=%2d (ASD=%2d Ctrl=%2d)",nm,nrow(md),sum(md$diagnosis=="ASD"),sum(md$diagnosis=="Control"))
 cm<-cm[rowSums(cm>=10)>=min(table(md$diagnosis)),]
 des<-if(nlevels(md$sex)>1) ~ sex+Astro+Neuro+diagnosis else ~ Astro+Neuro+diagnosis
 dds<-DESeq(DESeqDataSetFromMatrix(cm,md,des),quiet=TRUE)
 r<-as.data.frame(results(dds,contrast=c("diagnosis","ASD","Control")))
 data.frame(ens=sub("\\..*","",rownames(r)),beta=r$log2FoldChange,se=r$lfcSE,
            padj_within=r$padj,baseMean=r$baseMean,stringsAsFactors=FALSE)}

cat("===== per-cohort DESeq2 (confound-free cohorts) =====\n")
a<-per_cohort("GSE102741"); b<-per_cohort("PRJNA434002")
m<-merge(a,b,by="ens",suffixes=c("_1","_2"))
m<-m[is.finite(m$beta_1)&is.finite(m$beta_2)&is.finite(m$se_1)&is.finite(m$se_2)&m$se_1>0&m$se_2>0,]
logf("genes tested in BOTH cohorts: %d",nrow(m))

# ---- fixed-effect inverse variance ----
w1<-1/m$se_1^2; w2<-1/m$se_2^2
m$beta_FE<-(w1*m$beta_1+w2*m$beta_2)/(w1+w2)
m$se_FE  <-sqrt(1/(w1+w2))
m$z_FE   <-m$beta_FE/m$se_FE
m$p_FE   <-2*pnorm(-abs(m$z_FE))
m$padj_FE<-p.adjust(m$p_FE,"BH")

# ---- heterogeneity (Q, I^2) and DerSimonian-Laird random effects ----
Q  <- w1*(m$beta_1-m$beta_FE)^2 + w2*(m$beta_2-m$beta_FE)^2
df <- 1                                    # k-1 with k=2
C  <- (w1+w2) - (w1^2+w2^2)/(w1+w2)
tau2 <- pmax(0,(Q-df)/C)
m$I2 <- pmax(0,(Q-df)/Q)*100
wr1<-1/(m$se_1^2+tau2); wr2<-1/(m$se_2^2+tau2)
m$beta_RE<-(wr1*m$beta_1+wr2*m$beta_2)/(wr1+wr2)
m$se_RE  <-sqrt(1/(wr1+wr2))
m$p_RE   <-2*pnorm(-abs(m$beta_RE/m$se_RE))
m$padj_RE<-p.adjust(m$p_RE,"BH")

m$sym<-mapIds(org.Hs.eg.db,m$ens,"SYMBOL","ENSEMBL");m$sym[is.na(m$sym)]<-m$ens[is.na(m$sym)]
m$concordant<-sign(m$beta_1)==sign(m$beta_2)

logf("\n*** FIXED-EFFECT  : FDR<0.05 = %d | + |beta|>=1 = %d",
     sum(m$padj_FE<0.05,na.rm=TRUE), sum(m$padj_FE<0.05&abs(m$beta_FE)>=1,na.rm=TRUE))
logf("*** RANDOM-EFFECT : FDR<0.05 = %d | + |beta|>=1 = %d",
     sum(m$padj_RE<0.05,na.rm=TRUE), sum(m$padj_RE<0.05&abs(m$beta_RE)>=1,na.rm=TRUE))
logf("direction-concordant genes: %d/%d (%.1f%%)",sum(m$concordant),nrow(m),100*mean(m$concordant))
logf("median I^2: %.1f%% | genes with I^2>50%%: %d (%.1f%%)",
     median(m$I2,na.rm=TRUE), sum(m$I2>50,na.rm=TRUE), 100*mean(m$I2>50,na.rm=TRUE))

o<-m[order(m$p_FE),c("sym","ens","beta_1","se_1","beta_2","se_2","beta_FE","se_FE","p_FE","padj_FE",
                     "beta_RE","padj_RE","I2","concordant","baseMean_1","baseMean_2")]
write.csv(o,file.path(OUT,"DE_ASD_child_metaanalysis_confoundfree.csv"),row.names=FALSE)
cat("\n===== TOP 20 by fixed-effect p =====\n")
print(head(o[,c("sym","beta_1","beta_2","beta_FE","p_FE","padj_FE","I2","concordant")],20),row.names=FALSE)

key<-c("DLX1","DLX2","ARX","LHX6","SYNGAP1","DCX","SOX11","DRAXIN","DPYSL3","ST8SIA2")
cat("\n===== KEY GENES =====\n")
k<-o[o$sym%in%key,c("sym","beta_1","beta_2","beta_FE","se_FE","p_FE","padj_FE","I2","concordant")]
print(k[order(k$p_FE),],row.names=FALSE)
