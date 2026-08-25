## Cell-type marker set testing, 76 unique frontal donors. PRIMARY RESULT.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(limma); library(edgeR); library(AnnotationDbi); library(org.Hs.eg.db)})
logf <- function(...) cat(sprintf(...),"\n")
NP <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
CS <- file.path(NP,"counts_scz"); SHD <- file.path(CS,"sheets")
OR <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
RT <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/GSE80655/metadata/PRJNA319583_SraRunTable.csv"
OUT <- file.path(CS,"CELLTYPE_SETS"); dir.create(OUT, showWarnings=FALSE)
bt <- read.delim(file.path(OR,"gene_id2biotype.tsv"),header=FALSE,col.names=c("gid","bt"))
KEEP <- bt$gid[bt$bt %in% c("protein_coding","lncRNA")]
rc <- function(p){cm<-read.delim(p,comment.char="#",check.names=FALSE);g<-cm[[1]];cm[[1]]<-NULL
  cm<-cm[,!colnames(cm)%in%c("Chr","Start","End","Strand","Length"),drop=FALSE]
  colnames(cm)<-regmatches(colnames(cm),regexpr("SRR[0-9]+",colnames(cm)))
  m<-as.matrix(cm);rownames(m)<-g;storage.mode(m)<-"integer";m}
s1<-read.delim(file.path(SHD,"GSE80655_DLPFC_samplesheet.tsv"),check.names=FALSE)
s1<-merge(s1,read.delim(file.path(SHD,"GSE80655_DLPFC_prua.tsv"))[,c("srr","prua")],by="srr")
rownames(s1)<-s1$srr
c1<-rc(file.path(NP,"PRJNA319583/counts/PRJNA319583_gene_counts_matrix.tsv"))
k<-intersect(colnames(c1),rownames(s1));c1<-c1[,k];s1<-s1[k,]
rt<-read.csv(RT,check.names=FALSE);bio<-setNames(rt$BioSample,rt$Run)[colnames(c1)]
ix<-split(seq_len(ncol(c1)),bio)
c1s<-vapply(ix,function(i)as.integer(rowSums(c1[,i,drop=FALSE])),integer(nrow(c1)))
rownames(c1s)<-rownames(c1);ke<-!duplicated(bio)
s1<-s1[ke,];rownames(s1)<-bio[ke];c1<-c1s[,rownames(s1)]
s2<-read.delim(file.path(SHD,"GSE78936_samplesheet.tsv"),check.names=FALSE)
s2<-merge(s2,read.delim(file.path(SHD,"GSE78936_prua.tsv"))[,c("srr","prua")],by="srr")
rownames(s2)<-s2$srr;s2<-s2[s2$brain_region=="Brodmann's Area 11",]
c2<-rc(file.path(NP,"PRJNA314463/counts/PRJNA314463_gene_counts.txt"))
k<-intersect(colnames(c2),rownames(s2));c2<-c2[,k];s2<-s2[k,]
g<-intersect(rownames(c1),rownames(c2));g<-g[g%in%KEEP]
cnt<-cbind(c1[g,],c2[g,])
mt<-data.frame(dx=factor(c(s1$diagnosis,s2$diagnosis),levels=c("Control","SCZ")),
  study=factor(c(rep("GSE80655",ncol(c1)),rep("GSE78936",ncol(c2)))),
  prua=scale(c(s1$prua,s2$prua))[,1],row.names=colnames(cnt))
cnt<-cnt[rowSums(cnt>=10)>=min(table(mt$dx)),]
logf("frontal pool: %d donors (%s) | %d genes", nrow(mt),
     paste(names(table(mt$dx)),table(mt$dx),sep="=",collapse=" "), nrow(cnt))
y<-calcNormFactors(DGEList(cnt)); des<-model.matrix(~study+prua+dx,mt); v<-voom(y,des)
sym<-AnnotationDbi::mapIds(org.Hs.eg.db,sub("\\..*","",rownames(v)),"SYMBOL","ENSEMBL",multiVals="first")
sets<-list(
  GABAergic=c("SST","PVALB","GAD1","GAD2","VIP","NPY","CCK","LHX6","RELN","CNR1"),
  IEG=c("EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3","FOS","FOSB","JUN","JUNB","ARC","ATF3","DUSP1","DUSP6","NPAS4"),
  excitatory=c("SLC17A7","SLC17A6","CAMK2A","SATB2","TBR1","NEUROD6"),
  astrocyte=c("GFAP","ALDH1L1","AQP4","S100B","SLC1A2","SLC1A3","VIM"),
  oligodendrocyte=c("MBP","PLP1","MOG","MAG","CNP","SOX10"),
  CalbindinCR=c("CALB1","CALB2"))
idx<-lapply(sets,function(s) which(sym %in% s))
idx<-idx[sapply(idx,length)>=2]
cam<-camera(v, idx, des, contrast="dxSCZ")
IGC <- file.path(OUT,"camera_igc"); dir.create(IGC, showWarnings=FALSE)
cam_est  <- camera(v, idx, des, contrast="dxSCZ", inter.gene.cor=NA)
cam_zero <- camera(v, idx, des, contrast="dxSCZ", inter.gene.cor=0)
logf("\n===== CAMERA, correlation ESTIMATED (inter.gene.cor=NA) -- matches MDD =====")
print(cam_est);  write.csv(cam_est,  file.path(IGC,"CAMERA_estimated.csv"))
logf("\n===== CAMERA, inter.gene.cor=0 (true uncorrected) =====")
print(cam_zero); write.csv(cam_zero, file.path(IGC,"CAMERA_igc0.csv"))
logf("\n===== as archived, default inter.gene.cor=0.01 =====")
print(cam)
quit(save="no")
logf("\n=========== CAMERA (corrects inter-gene correlation) ===========")
print(cam); write.csv(cam, file.path(OUT,"CAMERA_celltype_sets.csv"))
logf("\n=========== uncorrected competitive test (for comparison) ===========")
fit<-eBayes(lmFit(v,des)); tt<-topTable(fit,coef="dxSCZ",number=Inf,sort.by="none")
tt$symbol<-sym[rownames(tt)]
st<-data.frame(set=names(idx), n=sapply(idx,length),
  median_t=sapply(idx,function(i) median(fit$t[i,"dxSCZ"])),
  p_uncorrected=sapply(idx,function(i)
    geneSetTest(i, fit$t[,"dxSCZ"], alternative="either", type="t", ranks.only=TRUE)))
print(st); write.csv(st, file.path(OUT,"uncorrected_setTest.csv"), row.names=FALSE)
write.csv(tt, file.path(OUT,"frontal_pool_voom_full.csv"))
logf("\n=========== PER-GENE, BY SET ===========")
allm<-do.call(rbind,lapply(names(sets),function(n){
  s<-tt[tt$symbol %in% sets[[n]],c("symbol","logFC","P.Value","adj.P.Val")]
  if(nrow(s)) cbind(set=n,s[order(s$P.Value),]) else NULL}))
print(allm); write.csv(allm, file.path(OUT,"marker_genes_per_set.csv"), row.names=FALSE)
gab<-tt[tt$symbol %in% sets$GABAergic,]
write.csv(gab[order(gab$P.Value),], file.path(OUT,"GABAergic_module.csv"))
logf("\nGABAergic module: %d genes, %d negative", nrow(gab), sum(gab$logFC<0))
logf("Outputs -> %s", OUT)
