.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
if (isTRUE(capabilities("cairo"))) options(bitmapType="cairo")
D   <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final"
OUT <- file.path(D,"Transcriptomics_ASD","figures"); dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
get <- function(path, pcol, sep=",") {
  if(!file.exists(path)){cat("missing:",path,"\n");return(NULL)}
  d <- read.delim(path, sep=sep, check.names=FALSE)
  if(!(pcol %in% names(d))){cat(basename(path),"cols:",paste(names(d),collapse=", "),"\n");return(NULL)}
  p <- as.numeric(d[[pcol]]); p[!is.na(p) & p>0] }
arms <- list(
  Autism        = get(file.path(D,"Transcriptomics_ASD/tables/DE_ASD_child_metaanalysis_confoundfree.csv"),"p_FE"),
  Schizophrenia = get(file.path(D,"Transcriptomics_SCZ/tables/pool_P1.csv"),"pvalue"),
  Depression    = get(file.path(D,"Transcriptomics_MDD/tables/01_PRIMARY_pooled_DE.tsv"),"pvalue","\t"))
arms <- arms[!sapply(arms,is.null)]
lam <- function(p) round(median(qchisq(1-p,1))/qchisq(0.5,1),3)
png(file.path(OUT,"FigSx_QQ_all_arms.png"), width=1500, height=520, res=150)
par(mfrow=c(1,max(1,length(arms))), mar=c(4.2,4.2,3,1))
for(nm in names(arms)){
  p <- sort(arms[[nm]]); n <- length(p)
  plot(-log10(ppoints(n)), -log10(p), pch=20, cex=.35, col="grey30",
       xlab="Expected -log10(P)", ylab="Observed -log10(P)",
       main=sprintf("%s  (n = %d)", nm, n))
  abline(0,1,col="firebrick",lwd=1.4)
  legend("topleft", bty="n", legend=paste("lambda =", lam(arms[[nm]])))
  cat(sprintf("%-14s n=%6d  lambda=%.3f  p<0.05: %d (%.2fx expected)\n",
      nm, n, lam(arms[[nm]]), sum(arms[[nm]]<0.05), sum(arms[[nm]]<0.05)/(0.05*n)))
}
dev.off(); cat("\nwritten to", OUT, "\n")
