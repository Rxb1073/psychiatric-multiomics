#!/usr/bin/env Rscript
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(DESeq2); library(org.Hs.eg.db); library(ggplot2); library(ggrepel)})
logf<-function(...) cat(sprintf(...),"\n")
AUT<-"/rds/projects/e/elhamsak-dissertation-autism"; AUT2<-"/rds/projects/e/elhamsak-dissertation-autism-2"
STG<-file.path(AUT,"RNA_seq/_pheno_staging"); DEC<-file.path(AUT2,"RNA_seq/Deconvolution/cell_type_proportions_MB.csv")
OUT<-file.path(AUT2,"RNA_seq/DE_analysis/POOLED_cortex"); PLT<-file.path(OUT,"plots"); dir.create(PLT,showWarnings=FALSE)
FDR<-0.05; LFC<-1
cohorts<-list(
 GSE102741=list(c=file.path(AUT,"RNA_seq/GSE102741/counts/GSE102741_gene_counts.txt"),s=file.path(STG,"GSE102741_samplesheet.tsv")),
 GSE59288_GSE51264=list(c=file.path(AUT,"RNA_seq/GSE59288_GSE51264/counts/GSE59288_GSE51264_gene_counts.txt"),s=file.path(STG,"GSE59288_GSE51264_samplesheet.tsv")),
 PRJNA434002=list(c=file.path(AUT2,"RNA_seq/PRJNA434002_bulk/counts/PRJNA434002_150bp_gene_counts.txt"),s=file.path(STG,"PRJNA434002_samplesheet.tsv")))
read_fc<-function(p){f<-readLines(p,n=1);sk<-if(grepl("^# ",f))1 else 0;m<-read.delim(p,header=TRUE,skip=sk,check.names=FALSE,row.names=1)
 d<-intersect(c("Chr","Start","End","Strand","Length"),colnames(m));if(length(d))m<-m[,setdiff(colnames(m),d),drop=FALSE]
 colnames(m)<-regmatches(colnames(m),regexpr("[SED]RR[0-9]+",colnames(m)));as.matrix(m)}
ndx<-function(x){x<-tolower(trimws(x));ifelse(grepl("asd|autism",x),"ASD",ifelse(grepl("control|ctrl|con\\b",x),"Control","OTHER"))}
prop<-read.csv(DEC,stringsAsFactors=FALSE);prop$srr<-gsub('"','',prop$srr)
key<-c("DLX1","DLX2","ARX","LHX6","SYNGAP1","DCX","SOX11","DRAXIN","DPYSL3","ST8SIA2","PCDH19","PCDHGC3")
build<-function(use){
 mats<-list();metas<-list()
 for(nm in use){cm<-read_fc(cohorts[[nm]]$c);ss<-read.delim(cohorts[[nm]]$s,check.names=FALSE,stringsAsFactors=FALSE)
  idc<-grep("^srr$",names(ss),ignore.case=TRUE,value=TRUE)[1];dxc<-grep("diagnos",names(ss),ignore.case=TRUE,value=TRUE)[1]
  sxc<-grep("^sex$|gender",names(ss),ignore.case=TRUE,value=TRUE)[1];agc<-grep("^age$",names(ss),ignore.case=TRUE,value=TRUE)[1]
  ss$.id<-regmatches(ss[[idc]],regexpr("[SED]RR[0-9]+",ss[[idc]]))
  md<-data.frame(srr=ss$.id,series=nm,diagnosis=ndx(ss[[dxc]]),sex=tolower(substr(trimws(ss[[sxc]]),1,1)),
    age=suppressWarnings(as.numeric(as.character(ss[[agc]]))),stringsAsFactors=FALSE)
  common<-intersect(colnames(cm),md$srr);cm<-cm[,common];md<-md[match(common,md$srr),];mats[[nm]]<-cm;metas[[nm]]<-md}
 genes<-Reduce(intersect,lapply(mats,rownames));counts<-do.call(cbind,lapply(mats,function(x)x[genes,,drop=FALSE]))
 meta<-do.call(rbind,metas);rownames(meta)<-meta$srr;meta<-meta[colnames(counts),]
 meta$Astro<-prop$Astrocytes[match(meta$srr,prop$srr)];meta$Neuro<-prop$Neurons[match(meta$srr,prop$srr)]
 ok<-meta$diagnosis%in%c("ASD","Control")&meta$sex%in%c("m","f")&!is.na(meta$age)&!is.na(meta$Astro)&meta$age<18
 counts<-round(counts[,ok]);m<-meta[ok,]
 pv<-function(mm) suppressWarnings(wilcox.test(age~diagnosis,mm)$p.value)
 while(pv(m)<=0.05){medA<-median(m$age[m$diagnosis=="ASD"]);medC<-median(m$age[m$diagnosis=="Control"])
  if(sum(m$diagnosis=="ASD")<=12||sum(m$diagnosis=="Control")<=12)break
  if(medA>=medC){idx<-which(m$diagnosis=="ASD");d<-idx[which.max(m$age[idx])]}else{idx<-which(m$diagnosis=="Control");d<-idx[which.min(m$age[idx])]}
  m<-m[-d,,drop=FALSE]}
 counts<-counts[,rownames(m)]
 m$diagnosis<-factor(as.character(m$diagnosis),levels=c("Control","ASD"));m$series<-droplevels(factor(m$series));m$sex<-droplevels(factor(m$sex))
 for(v in c("Astro","Neuro"))m[[v]]<-as.numeric(scale(m[[v]]))
 counts<-counts[rowSums(counts>=10)>=min(table(m$diagnosis)),]
 des<-if(nlevels(m$series)>1) ~ series+sex+Astro+Neuro+diagnosis else ~ sex+Astro+Neuro+diagnosis
 dds<-DESeq(DESeqDataSetFromMatrix(counts,m,des),quiet=TRUE)
 r<-as.data.frame(results(dds,contrast=c("diagnosis","ASD","Control"),alpha=FDR));r$ens<-sub("\\..*","",rownames(r))
 r$sym<-mapIds(org.Hs.eg.db,r$ens,"SYMBOL","ENSEMBL");r$sym[is.na(r$sym)]<-r$ens[is.na(r$sym)]
 list(m=m,r=r)}
volc<-function(o,tag,sub){
 r<-o$r; v<-r[!is.na(r$padj),]
 # STANDARD scheme: colour encodes DIRECTION among significant+|LFC| genes
 v$diffexpressed<-"NO"
 v$diffexpressed[v$padj<FDR & v$log2FoldChange>=LFC]<-"UP"
 v$diffexpressed[v$padj<FDR & v$log2FoldChange<=-LFC]<-"DOWN"
 v$diffexpressed<-factor(v$diffexpressed,levels=c("DOWN","UP","NO"))
 nUp<-sum(v$diffexpressed=="UP"); nDn<-sum(v$diffexpressed=="DOWN")
 logf("[%s] DOWN=%d UP=%d",tag,nDn,nUp)
 lab<-v[v$sym%in%key & v$diffexpressed!="NO",]; lab<-unique(rbind(lab,head(v[order(v$padj),],12)))
 p<-ggplot(v,aes(log2FoldChange,-log10(padj),color=diffexpressed))+geom_point(size=1.1,alpha=.6)+
  scale_color_manual(values=c(DOWN="#2c7bb6",UP="#d7191c",NO="grey72"),name="diffexpressed")+
  geom_vline(xintercept=c(-LFC,LFC),linetype=2,color="grey45")+geom_hline(yintercept=-log10(FDR),linetype=2,color="grey45")+
  geom_text_repel(data=lab,aes(label=sym),size=3,max.overlaps=25,show.legend=FALSE,seed=1)+
  labs(title="ASD vs Control - childhood cortex (age-matched, composition-adjusted)",
       subtitle=sprintf("%s\nDOWN in ASD: %d   UP in ASD: %d   (FDR<0.05 & |log2FC|>=1)",sub,nDn,nUp),
       x="log2 fold change (ASD / Control)",y="-log10 adjusted p")+
  theme_bw(base_size=12)+theme(legend.position="right",plot.title=element_text(size=12.5),
                               plot.subtitle=element_text(size=10.5,lineheight=1.05))
 ggsave(file.path(PLT,sprintf("volcano_%s_updown.png",tag)),p,width=8.3,height=6.8,dpi=200,type="cairo")
 ggsave(file.path(PLT,sprintf("volcano_%s_updown.pdf",tag)),p,width=8.3,height=6.8)}
cat("=== POOLED ===\n"); o3<-build(names(cohorts)); volc(o3,"ASD_child_3dataset",sprintf("pooled 3 cohorts, n=%d; NOTE signal is one-cohort-driven (see per-cohort analysis)",nrow(o3$m)))
cat("=== GSE59288 ONLY ===\n"); o1<-build("GSE59288_GSE51264"); volc(o1,"ASD_child_GSE59288only",sprintf("GSE59288/51264 only, n=%d (cohort carrying the signal)",nrow(o1$m)))
logf("done -> %s",PLT)
