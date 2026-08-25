## Composition-adjusted CAMERA, 76 unique frontal donors.
## Derived from scz_camera_igc_check.R; data-loading block unchanged.
## Asks whether the GABAergic downshift survives adjustment for cell-type
## proportion. All models use inter.gene.cor=NA (estimated), matching MDD.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(limma); library(edgeR); library(AnnotationDbi); library(org.Hs.eg.db)})
logf <- function(...) cat(sprintf(...),"\n")
NP <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
CS <- file.path(NP,"counts_scz"); SHD <- file.path(CS,"sheets")
OR <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
RT <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/GSE80655/metadata/PRJNA319583_SraRunTable.csv"
COMP <- file.path(OR,"DISSERTATION_final/Transcriptomics_SCZ/integration/SCZ_composition_donorcorrect_wide.csv")
OUT <- file.path(OR,"DISSERTATION_final/Transcriptomics_SCZ/tables/camera_composition")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

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

cp <- read.csv(COMP, check.names=FALSE, row.names=1)
miss_v <- setdiff(rownames(mt), rownames(cp))
if(length(miss_v)) stop("donor IDs do not match: ", paste(head(miss_v,5),collapse=", "))
cp <- cp[rownames(mt),,drop=FALSE]
logf("joined %d donors, both directions clean", nrow(cp))
mt$MB_neuronal <- scale(cp$MB_Neurons)[,1]
mt$CA_neuronal <- scale(cp$CA_Excitatory + cp$CA_Inhibitory)[,1]
mt$CA_inhib    <- scale(cp$CA_Inhibitory)[,1]

y <- calcNormFactors(DGEList(cnt))
sym <- AnnotationDbi::mapIds(org.Hs.eg.db, sub("\\..*","",rownames(y)), "SYMBOL","ENSEMBL", multiVals="first")
sets<-list(
  GABAergic=c("SST","PVALB","GAD1","GAD2","VIP","NPY","CCK","LHX6","RELN","CNR1"),
  IEG=c("EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3","FOS","FOSB","JUN","JUNB","ARC","ATF3","DUSP1","DUSP6","NPAS4"),
  excitatory=c("SLC17A7","SLC17A6","CAMK2A","SATB2","TBR1","NEUROD6"),
  astrocyte=c("GFAP","ALDH1L1","AQP4","S100B","SLC1A2","SLC1A3","VIM"),
  oligodendrocyte=c("MBP","PLP1","MOG","MAG","CNP","SOX10"),
  CalbindinCR=c("CALB1","CALB2"))
idx <- lapply(sets, function(s) which(sym %in% s))
idx <- idx[sapply(idx,length)>=2]

models <- list(
  baseline      = ~study+prua+dx,
  MB_neuronal   = ~study+prua+MB_neuronal+dx,
  CA_neuronal   = ~study+prua+CA_neuronal+dx,
  CA_inhibitory = ~study+prua+CA_inhib+dx)

res <- list()
for(nm in names(models)){
  des <- model.matrix(models[[nm]], mt)
  stopifnot(tail(colnames(des),1)=="dxSCZ")
  v <- voom(y, des)
  cam <- camera(v, idx, des, contrast="dxSCZ", inter.gene.cor=NA)
  logf("\n===== %s : %s =====", nm, paste(deparse(models[[nm]]),collapse=""))
  print(cam)
  write.csv(cam, file.path(OUT, sprintf("CAMERA_%s.csv", nm)))
  res[[nm]] <- data.frame(model=nm, set=rownames(cam), cam[,c("NGenes","Correlation","Direction","PValue","FDR")])
}
comb <- do.call(rbind, res); rownames(comb) <- NULL
write.csv(comb, file.path(OUT,"CAMERA_composition_comparison.csv"), row.names=FALSE)
logf("\n===== GABAergic across the four models =====")
print(comb[comb$set=="GABAergic", c("model","Correlation","Direction","PValue","FDR")])
logf("\n===== IEG across the four models =====")
print(comb[comb$set=="IEG", c("model","Direction","PValue","FDR")])
logf("\nwritten to %s", OUT)
