.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(library(vegan)); set.seed(1)
P   <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/phylogenetic"
CM  <- file.path(P,"core_metrics_phylo")
MD  <- file.path(P,"phylo_metadata.tsv")
OUT <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final/Microbiome_MDD/tables"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

md <- read.delim(MD, check.names=FALSE, stringsAsFactors=FALSE)
if (grepl("^#", md[1,1])) md <- md[-1,,drop=FALSE]
rownames(md) <- md[[1]]
cat("metadata:", nrow(md), "samples\ncolumns:", paste(names(md),collapse=", "), "\n")
dxc <- grep("disease|diagnos|group|condition", names(md), ignore.case=TRUE, value=TRUE)[1]
stc <- grep("study|cohort|bioproject|project", names(md), ignore.case=TRUE, value=TRUE)[1]
cat("diagnosis column:", dxc, " | study column:", stc, "\n")
md$dx <- factor(md[[dxc]]); md$st <- factor(md[[stc]])
print(table(md$st, md$dx))

rd <- function(q){ t<-tempfile(); dir.create(t); unzip(q, exdir=t)
  f<-list.files(t,pattern="distance-matrix\\.tsv$",recursive=TRUE,full.names=TRUE)[1]
  m<-read.delim(f,row.names=1,check.names=FALSE); unlink(t,recursive=TRUE)
  as.dist(as.matrix(m)) }

one <- function(d,lab){ keep<-intersect(labels(d), rownames(md))
  cat(sprintf("  [%s] n = %d\n", lab, length(keep)))
  dd<-as.dist(as.matrix(d)[keep,keep]); sub<-md[keep,]
  pn<-adonis2(dd ~ st + dx, data=sub, permutations=999, by="margin")
  bd<-betadisper(dd, droplevels(sub$dx)); pd<-permutest(bd, permutations=999)
  data.frame(metric=lab, n=length(keep),
             R2=round(pn["dx","R2"],4), pseudo_F=round(pn["dx","F"],3),
             P_permanova=pn["dx","Pr(>F)"], P_permdisp=pd$tab$`Pr(>F)`[1]) }

mats <- c(bray_curtis="bray_curtis_distance_matrix.qza", jaccard="jaccard_distance_matrix.qza",
          unweighted_unifrac="unweighted_unifrac_distance_matrix.qza",
          weighted_unifrac="weighted_unifrac_distance_matrix.qza")
res <- list()
for(nm in names(mats)){ f<-file.path(CM,mats[[nm]])
  if(file.exists(f)) res[[nm]]<-one(rd(f),nm) else cat("missing:",f,"\n") }
S <- do.call(rbind,res); rownames(S)<-NULL
cat("\n=== MDD BETA DIVERSITY, adonis2 (~ study + disease, marginal) ===\n")
print(S, row.names=FALSE)
write.csv(S, file.path(OUT,"MDD_beta_diversity_adonis2.csv"), row.names=FALSE)
cat("\nwritten to", OUT, "\n")
