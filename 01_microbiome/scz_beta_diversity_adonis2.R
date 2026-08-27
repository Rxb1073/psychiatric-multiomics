.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(library(vegan)); set.seed(1)
S   <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA880407"
CM  <- file.path(S,"qiime2_v4","core-metrics-results")
MD  <- file.path(S,"metadata","metadata_qiime2.tsv")
GEN <- file.path(S,"genus_raw.tsv")
OUT <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final/Microbiome_SCZ/tables"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

md <- read.delim(MD, check.names=FALSE, stringsAsFactors=FALSE)
if (grepl("^#", md[1,1])) md <- md[-1,,drop=FALSE]
rownames(md) <- md[[1]]
cat("metadata:", nrow(md), "samples\ncolumns:", paste(names(md),collapse=", "), "\n")
dxcol <- grep("diagnos|disease|group|condition|phenotype", names(md), ignore.case=TRUE, value=TRUE)[1]
cat("using column:", dxcol, "\n"); md$dx <- factor(md[[dxcol]]); print(table(md$dx))

read_qza_dist <- function(qza){ tmp<-tempfile(); dir.create(tmp); unzip(qza, exdir=tmp)
  f <- list.files(tmp, pattern="distance-matrix\\.tsv$", recursive=TRUE, full.names=TRUE)[1]
  m <- read.delim(f, row.names=1, check.names=FALSE); unlink(tmp, recursive=TRUE)
  as.dist(as.matrix(m)) }

run_one <- function(d,label){ ids<-labels(d); keep<-intersect(ids, rownames(md))
  if(length(keep)<length(ids)) cat(sprintf("  [%s] %d of %d matched\n",label,length(keep),length(ids)))
  dd <- as.dist(as.matrix(d)[keep,keep]); g <- droplevels(md[keep,"dx"])
  pn <- adonis2(dd ~ g, permutations=999); bd <- betadisper(dd,g)
  pd <- permutest(bd, permutations=999)
  data.frame(metric=label, n=length(keep), R2=round(pn$R2[1],4), pseudo_F=round(pn$F[1],3),
             P_permanova=pn$`Pr(>F)`[1], P_permdisp=pd$tab$`Pr(>F)`[1], stringsAsFactors=FALSE) }

mats <- c(bray_curtis="bray_curtis_distance_matrix.qza", jaccard="jaccard_distance_matrix.qza",
          unweighted_unifrac="unweighted_unifrac_distance_matrix.qza",
          weighted_unifrac="weighted_unifrac_distance_matrix.qza")
res <- list()
for(nm in names(mats)){ f<-file.path(CM,mats[[nm]])
  if(file.exists(f)) res[[nm]] <- run_one(read_qza_dist(f), nm) else cat("missing:",f,"\n") }

if(file.exists(GEN)){ g<-read.delim(GEN,row.names=1,check.names=FALSE,skip=1,comment.char="")
  g<-g[,colnames(g) %in% rownames(md),drop=FALSE]; g<-g[rowSums(g)>0,,drop=FALSE]
  x<-t(as.matrix(g))+0.5; clr<-t(apply(x,1,function(r) log(r)-mean(log(r))))
  res[["aitchison"]] <- run_one(dist(clr,method="euclidean"),"aitchison") }

S2 <- do.call(rbind,res); rownames(S2)<-NULL
cat("\n=== SCZ BETA DIVERSITY, adonis2 ===\n"); print(S2, row.names=FALSE)
write.csv(S2, file.path(OUT,"SCZ_beta_diversity_adonis2.csv"), row.names=FALSE)
cat("\nwritten to", OUT, "\n")
