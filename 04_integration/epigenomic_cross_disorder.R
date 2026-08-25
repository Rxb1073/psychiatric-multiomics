## Cross-disorder comparison at the EPIGENOMIC layer: do the autism and
## schizophrenia epigenomic gene sets overlap, and does either arm's set move in
## the other arm's transcriptome?
##
## This is the third cross-disorder axis and the only one not previously tested.
## The two epigenomic layers are different assay types -- CHD8 binding sites in
## autism, differentially methylated positions in schizophrenia -- so they cannot
## be compared as effect sizes. Two questions are available instead:
##
##   1. Do the two gene sets share more genes than chance? A shared gene would
##      be one both bound by CHD8 and differentially methylated in schizophrenia.
##   2. Does each arm's epigenomic set move in the OTHER arm's transcriptome?
##      This is the cross-disorder form of the within-arm tests already reported.
##
## Both use the expression-matched permutation of the integration chapter, since
## both sets are assay-defined and therefore biased toward abundant transcripts.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(IlluminaHumanMethylation450kanno.ilmn12.hg19)})
set.seed(1)

D   <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final"
M   <- "/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-methylation"
INT <- file.path(D,"INTEGRATION_inputs")
OUT <- file.path(INT,"13_epigenomic_cross_disorder"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
NPERM <- 10000

## ---- set 1: CHD8 targets (autism) ------------------------------------------
p2g <- read.csv(file.path(INT,"07_chd8_axis","CHD8_peak_to_gene.csv"), stringsAsFactors=FALSE)
cat(sprintf("CHD8 targets: %d genes (Ensembl IDs)\n", length(unique(p2g$ens))))

## ---- set 2: schizophrenia DMP genes ----------------------------------------
dmp <- read.csv(file.path(M,"GSE74193/results/DMCpGs_SCZ_postnatal_autosomal_213.csv"),
                row.names=1, check.names=FALSE)
ann <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
a <- ann[rownames(dmp), c("UCSC_RefGene_Name","UCSC_RefGene_Group")]
dmp_genes <- unique(unlist(strsplit(a$UCSC_RefGene_Name, ";")))
dmp_genes <- dmp_genes[dmp_genes != "" & !is.na(dmp_genes)]
cat(sprintf("SCZ DMP genes: %d (symbols)\n", length(dmp_genes)))

## ---- transcriptomes --------------------------------------------------------
asd <- read.csv(file.path(INT,"05_genelevel","ASD_meta.csv"), check.names=FALSE)
asd <- data.frame(ens=sub("\\..*$","",asd$ens), sym=asd$sym, beta=asd$beta_FE,
                  abund=rowMeans(asd[,c("baseMean_1","baseMean_2")], na.rm=TRUE),
                  stringsAsFactors=FALSE)
scz <- read.csv(file.path(INT,"05_genelevel","SCZ_frontal.csv"), check.names=FALSE)
names(scz)[1] <- "ens"
scz <- data.frame(ens=sub("\\..*$","",scz$ens), sym=scz$symbol, beta=scz$logFC,
                  abund=scz$AveExpr, stringsAsFactors=FALSE)
for(nm in c("asd","scz")){
  x <- get(nm); x <- x[!is.na(x$beta) & !is.na(x$abund) & !duplicated(x$ens), ]
  x$stratum <- cut(rank(x$abund, ties.method="first"),
                   breaks=quantile(rank(x$abund, ties.method="first"), probs=seq(0,1,0.05)),
                   include.lowest=TRUE, labels=FALSE)
  assign(nm, x) }
cat(sprintf("ASD transcriptome %d genes | SCZ transcriptome %d genes\n", nrow(asd), nrow(scz)))

## ---- question 1: do the two epigenomic sets overlap? -----------------------
cat("\n", strrep("=",70), "\nOVERLAP BETWEEN THE TWO EPIGENOMIC SETS\n", strrep("=",70), "\n", sep="")
chd8_sym <- unique(asd$sym[asd$ens %in% p2g$ens])
chd8_sym <- chd8_sym[!is.na(chd8_sym)]
pool <- intersect(asd$sym, scz$sym); pool <- pool[!is.na(pool)]
A <- intersect(chd8_sym, pool); B <- intersect(dmp_genes, pool)
ov <- intersect(A, B)
ph <- phyper(length(ov)-1, length(A), length(pool)-length(A), length(B), lower.tail=FALSE)
cat(sprintf("shared gene universe: %d\nCHD8 targets in universe: %d\nSCZ DMP genes in universe: %d\n",
            length(pool), length(A), length(B)))
cat(sprintf("overlap: %d (%s)\nhypergeometric p = %.4f | %.1f-fold\n",
            length(ov), paste(head(ov,20), collapse=", "), ph,
            (length(ov)/length(B)) / (length(A)/length(pool))))

## ---- question 2: does each set move in the other arm's transcriptome? ------
matched_test <- function(x, genes, label){
  idx <- which(x$sym %in% genes | x$ens %in% genes)
  if(length(idx) < 20){ cat(sprintf("\n%s: %d genes -- too few\n", label, length(idx))); return(NULL) }
  mt <- median(x$beta[idx])
  tab <- table(x$stratum[idx]); rest <- setdiff(seq_len(nrow(x)), idx)
  pool2 <- split(rest, x$stratum[rest])
  null <- replicate(NPERM, {
    s <- unlist(lapply(names(tab), function(k)
      if(!is.null(pool2[[k]])) sample(pool2[[k]], min(tab[[k]], length(pool2[[k]]))) else integer(0)))
    median(x$beta[s]) })
  p <- (sum(abs(null - median(null)) >= abs(mt - median(null))) + 1)/(NPERM + 1)
  cat(sprintf("\n%s\n  %d genes | median %+.4f | matched null %+.4f (SD %.4f) | p = %.4f | %s\n",
              label, length(idx), mt, median(null), sd(null), p,
              if(mt > median(null)) "UP" else "DOWN"))
  data.frame(test=label, n=length(idx), median=mt, null=median(null), p=p)
}

cat("\n", strrep("=",70), "\nEACH SET IN THE OTHER ARM'S TRANSCRIPTOME\n", strrep("=",70), "\n", sep="")
R <- list()
R$a <- matched_test(scz, chd8_sym,  "CHD8 targets (ASD layer) in the SCHIZOPHRENIA transcriptome")
R$b <- matched_test(asd, dmp_genes, "SCZ DMP genes (SCZ layer) in the AUTISM transcriptome")
cat("\nfor reference, the within-arm results already reported:\n")
cat("  CHD8 targets in the autism transcriptome      : displaced UP, p < 0.0001\n")
cat("  SCZ DMP genes in the schizophrenia transcriptome: no coherent relationship\n")

S <- do.call(rbind, R)
if(!is.null(S)) write.csv(S, file.path(OUT,"epigenomic_cross_disorder.csv"), row.names=FALSE)
write.csv(data.frame(gene=ov), file.path(OUT,"shared_epigenomic_genes.csv"), row.names=FALSE)
cat(sprintf("\nwritten to %s\n", OUT))
