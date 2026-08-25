#!/usr/bin/env Rscript
# Step 7 (exploratory) + step 9 (MA plots) for the PRIMARY confound-free pool.
# PCA, sample clustering, outlier check on VST expression; MA plot per cohort.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
set.seed(42)
suppressMessages({library(DESeq2); library(ggplot2); library(ggrepel); library(pheatmap)})
logf<-function(...) cat(sprintf(...),"\n")
AUT<-"/rds/projects/e/elhamsak-dissertation-autism"; AUT2<-"/rds/projects/e/elhamsak-dissertation-autism-2"
STG<-file.path(AUT,"RNA_seq/_pheno_staging"); DEC<-file.path(AUT2,"RNA_seq/Deconvolution/cell_type_proportions_MB.csv")
PLT<-file.path(AUT2,"RNA_seq/DE_analysis/DISSERTATION_final/figures"); dir.create(PLT,showWarnings=FALSE,recursive=TRUE)
QC_FAIL<-c("SRR5938468","SRR5938447","SRR5938445","SRR1510129","SRR5938444","SRR5938441","SRR1004284","SRR5938466")
cohorts<-list(
 GSE102741=list(c=file.path(AUT,"RNA_seq/GSE102741/counts/GSE102741_gene_counts.txt"),s=file.path(STG,"GSE102741_samplesheet.tsv")),
 PRJNA434002=list(c=file.path(AUT2,"RNA_seq/PRJNA434002_bulk/counts/PRJNA434002_150bp_gene_counts.txt"),s=file.path(STG,"PRJNA434002_samplesheet.tsv")))
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
ok<-meta$diagnosis%in%c("ASD","Control")&meta$sex%in%c("m","f")&!is.na(meta$age)&!is.na(meta$Astro)&
    meta$age<18&!(meta$srr%in%QC_FAIL)
counts<-round(counts[,ok]);m<-meta[ok,]
m$diagnosis<-factor(m$diagnosis,levels=c("Control","ASD"));m$series<-factor(m$series);m$sex<-factor(m$sex)
counts<-counts[rowSums(counts>=10)>=min(table(m$diagnosis)),]
logf("primary pool: n=%d (ASD=%d Ctrl=%d), genes=%d",nrow(m),sum(m$diagnosis=="ASD"),sum(m$diagnosis=="Control"),nrow(counts))

dds<-DESeqDataSetFromMatrix(counts,m,~ series+sex+diagnosis)
dds<-estimateSizeFactors(dds)
vsd<-vst(dds,blind=TRUE)          # blind=TRUE for unsupervised QC
v<-assay(vsd)

# ---------- PCA ----------
pv<-rowVars(v); top<-order(pv,decreasing=TRUE)[1:min(2000,length(pv))]
pca<-prcomp(t(v[top,]),scale.=FALSE)
pct<-round(100*pca$sdev^2/sum(pca$sdev^2),1)
d<-data.frame(PC1=pca$x[,1],PC2=pca$x[,2],PC3=pca$x[,3],
              Diagnosis=m$diagnosis,Cohort=m$series,Age=m$age,Sex=m$sex,srr=rownames(m))
logf("PC1 %.1f%% | PC2 %.1f%% | PC3 %.1f%% of variance",pct[1],pct[2],pct[3])
# association of PCs with covariates
for(pc in c("PC1","PC2","PC3")){
  a<-summary(aov(d[[pc]]~d$Cohort))[[1]][["Pr(>F)"]][1]
  b<-summary(aov(d[[pc]]~d$Diagnosis))[[1]][["Pr(>F)"]][1]
  cc<-cor.test(d[[pc]],d$Age)$p.value
  logf("%s ~ cohort p=%.3g | diagnosis p=%.3g | age p=%.3g",pc,a,b,cc)}
p<-ggplot(d,aes(PC1,PC2,color=Cohort,shape=Diagnosis))+geom_point(size=3,alpha=.85)+
 labs(title="PCA - primary confound-free pool (VST, top 2000 variable genes)",
      subtitle=sprintf("n=%d; PC1 %.1f%%, PC2 %.1f%% of variance",nrow(m),pct[1],pct[2]),
      x=sprintf("PC1 (%.1f%%)",pct[1]),y=sprintf("PC2 (%.1f%%)",pct[2]))+
 theme_bw(base_size=12)
ggsave(file.path(PLT,"PCA_primary_confoundfree.png"),p,width=8,height=6.5,dpi=200)
ggsave(file.path(PLT,"PCA_primary_confoundfree.pdf"),p,width=8,height=6.5)

# ---------- sample-sample distance clustering ----------
sd<-as.matrix(dist(t(v[top,])))
ann<-data.frame(Diagnosis=m$diagnosis,Cohort=m$series,row.names=rownames(m))
png(file.path(PLT,"sample_distance_primary.png"),width=1100,height=1000,res=140,type="cairo")
pheatmap(sd,annotation_col=ann,annotation_row=ann,show_rownames=FALSE,show_colnames=FALSE,
  main="Sample-to-sample distance (VST) - primary pool",fontsize=9)
dev.off()

# ---------- outlier check: median correlation to other samples ----------
cm2<-cor(v[top,],method="spearman")
mc<-apply(cm2,1,function(x) median(x[-which.max(x)]))
out<-data.frame(srr=names(mc),median_cor=round(mc,3),Cohort=m$series,Diagnosis=m$diagnosis)
out<-out[order(out$median_cor),]
thr<-median(mc)-3*mad(mc)
logf("\noutlier screen: median inter-sample rho = %.3f; flag threshold (med-3*MAD) = %.3f",median(mc),thr)
logf("samples below threshold: %d",sum(mc<thr))
print(head(out,8),row.names=FALSE)
write.csv(out,file.path(PLT,"..","tables","sample_outlier_screen_primary.csv"),row.names=FALSE)

# ---------- MA plots (per cohort) ----------
for(nm in levels(m$series)){
  k<-m$series==nm; cc<-counts[,k]; mm<-droplevels(m[k,])
  for(v2 in c("Astro","Neuro")) mm[[v2]]<-as.numeric(scale(mm[[v2]]))
  cc<-cc[rowSums(cc>=10)>=min(table(mm$diagnosis)),]
  des<-if(nlevels(mm$sex)>1) ~ sex+Astro+Neuro+diagnosis else ~ Astro+Neuro+diagnosis
  d2<-DESeq(DESeqDataSetFromMatrix(cc,mm,des),quiet=TRUE)
  r<-as.data.frame(results(d2,contrast=c("diagnosis","ASD","Control")))
  r$sig<-!is.na(r$padj)&r$padj<0.05
  p<-ggplot(r,aes(baseMean,log2FoldChange,color=sig))+geom_point(size=.9,alpha=.5)+
   scale_x_log10()+scale_color_manual(values=c("FALSE"="grey72","TRUE"="#d7191c"),name="FDR<0.05")+
   geom_hline(yintercept=0,linetype=2,color="grey40")+
   labs(title=sprintf("MA plot - %s (childhood, composition-adjusted)",nm),
        x="mean of normalised counts",y="log2 fold change (ASD / Control)")+
   theme_bw(base_size=12)
  ggsave(file.path(PLT,sprintf("MA_%s.png",nm)),p,width=7.5,height=5.5,dpi=200)
  ggsave(file.path(PLT,sprintf("MA_%s.pdf",nm)),p,width=7.5,height=5.5)
  logf("wrote MA_%s",nm)}
logf("\ndone -> %s",PLT)
