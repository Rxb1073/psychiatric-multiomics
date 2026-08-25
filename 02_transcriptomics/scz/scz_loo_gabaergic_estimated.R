## Leave-one-gene-out for the GABAergic CAMERA result (Dr Hosny's request)
## Re-runs CAMERA 10 times, each dropping one marker, to confirm no single
## gene drives FDR = 8e-11. Reads the saved voom fit; recomputes on the pool.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(limma); library(edgeR); library(AnnotationDbi); library(org.Hs.eg.db)})
logf <- function(...) cat(sprintf(...),"\n")
NP <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
CS <- file.path(NP,"counts_scz"); SHD <- file.path(CS,"sheets")
OR <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
RT <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/GSE80655/metadata/PRJNA319583_SraRunTable.csv"
OUT <- file.path(CS,"CELLTYPE_SETS")
bt <- read.delim(file.path(OR,"gene_id2biotype.tsv"),header=FALSE,col.names=c("gid","bt"))
KEEP <- bt$gid[bt$bt %in% c("protein_coding","lncRNA")]
rc <- function(p){cm<-read.delim(p,comment.char="#",check.names=FALSE);g<-cm[[1]];cm[[1]]<-NULL
  cm<-cm[,!colnames(cm)%in%c("Chr","Start","End","Strand","Length"),drop=FALSE]
  colnames(cm)<-regmatches(colnames(cm),regexpr("SRR[0-9]+",colnames(cm)))
  m<-as.matrix(cm);rownames(m)<-g;storage.mode(m)<-"integer";m}
s1<-read.delim(file.path(SHD,"GSE80655_DLPFC_samplesheet.tsv"),check.names=FALSE)
s1<-merge(s1,read.delim(file.path(SHD,"GSE80655_DLPFC_prua.tsv"))[,c("srr","prua")],by="srr");rownames(s1)<-s1$srr
c1<-rc(file.path(NP,"PRJNA319583/counts/PRJNA319583_gene_counts_matrix.tsv"))
k<-intersect(colnames(c1),rownames(s1));c1<-c1[,k];s1<-s1[k,]
rt<-read.csv(RT,check.names=FALSE);bio<-setNames(rt$BioSample,rt$Run)[colnames(c1)]
ix<-split(seq_len(ncol(c1)),bio)
c1s<-vapply(ix,function(i)as.integer(rowSums(c1[,i,drop=FALSE])),integer(nrow(c1)))
rownames(c1s)<-rownames(c1);ke<-!duplicated(bio);s1<-s1[ke,];rownames(s1)<-bio[ke];c1<-c1s[,rownames(s1)]
s2<-read.delim(file.path(SHD,"GSE78936_samplesheet.tsv"),check.names=FALSE)
s2<-merge(s2,read.delim(file.path(SHD,"GSE78936_prua.tsv"))[,c("srr","prua")],by="srr");rownames(s2)<-s2$srr
s2<-s2[s2$brain_region=="Brodmann's Area 11",]
c2<-rc(file.path(NP,"PRJNA314463/counts/PRJNA314463_gene_counts.txt"))
k<-intersect(colnames(c2),rownames(s2));c2<-c2[,k];s2<-s2[k,]
g<-intersect(rownames(c1),rownames(c2));g<-g[g%in%KEEP]
cnt<-cbind(c1[g,],c2[g,])
mt<-data.frame(dx=factor(c(s1$diagnosis,s2$diagnosis),levels=c("Control","SCZ")),
  study=factor(c(rep("A",ncol(c1)),rep("B",ncol(c2)))),prua=scale(c(s1$prua,s2$prua))[,1],row.names=colnames(cnt))
cnt<-cnt[rowSums(cnt>=10)>=min(table(mt$dx)),]
y<-calcNormFactors(DGEList(cnt));des<-model.matrix(~study+prua+dx,mt);v<-voom(y,des)
sym<-AnnotationDbi::mapIds(org.Hs.eg.db,sub("\\..*","",rownames(v)),"SYMBOL","ENSEMBL",multiVals="first")
GABA<-c("SST","PVALB","GAD1","GAD2","VIP","NPY","CCK","LHX6","RELN","CNR1")
full<-which(sym %in% GABA)
r0<-camera(v, list(GABAergic=full), des, contrast="dxSCZ", inter.gene.cor=NA)
logf("full set (n=%d): FDR = %.2e  %s", length(full), r0$PValue, r0$Direction)
logf("\n---- leave-one-gene-out ----")
res<-data.frame()
for (drop in GABA){
  idx<-which(sym %in% setdiff(GABA,drop))
  r<-camera(v, list(s=idx), des, contrast="dxSCZ", inter.gene.cor=NA)
  res<-rbind(res,data.frame(dropped=drop,n=length(idx),FDR=r$PValue,dir=r$Direction))
  logf("  drop %-6s (n=%d): p = %.2e  %s", drop, length(idx), r$PValue, r$Direction)
}
logf("\nworst case (largest p after dropping any one gene): %.2e  [dropping %s]",
     max(res$FDR), res$dropped[which.max(res$FDR)])
logf("all ten remain Down: %s", all(res$dir=="Down"))
write.csv(res, file.path(OUT,"GABAergic_leave_one_out_estimated.csv"), row.names=FALSE)
logf("saved -> %s", file.path(OUT,"GABAergic_leave_one_out_estimated.csv"))
