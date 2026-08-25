## SCZ cross-omics link: does prefrontal DNA methylation align with prefrontal
## gene expression in schizophrenia?
##
## LAYERS. Methylation: 213 differentially methylated positions (FDR<0.05,
## |delta beta| >= 0.05) from GSE74193, postnatal autosomal analysis. Expression:
## pooled frontal limma-voom topTable, 76 donors, GSE80655 + GSE78936. Different
## cohorts, no shared donors -- the layers meet at the gene, not the sample.
##
## DIRECTION IS THE POINT. An overlap test would only ask whether the same genes
## appear twice. Methylation carries a sign, so a specific prediction is available
## and is fixed BEFORE looking at the result:
##   promoter (TSS200/TSS1500/5'UTR/1stExon) hypermethylation -> DOWNregulation
##   gene body hypermethylation              -> UPregulation (transcription-coupled)
## These are opposite, so promoter and body probes are tested SEPARATELY. Pooling
## them would cancel a real effect. 87% of the 213 positions are hypomethylated,
## so the promoter set predicts upregulation and the body set downregulation.
##
## TEST. Expression-matched permutation, as used for the CHD8 axis: methylation
## arrays interrogate CpG islands at active promoters, so DMP genes are not a
## random draw with respect to expression level. Comparing them against ALL other
## genes would confound methylation with abundance. cameraPR is reported alongside
## for comparison, across assumed inter-gene correlations.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(limma); library(IlluminaHumanMethylation450kanno.ilmn12.hg19)})
set.seed(1)

M   <- "/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-methylation"
INT <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final/INTEGRATION_inputs"
OUT <- file.path(INT, "08_methylation_axis"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
NPERM <- 10000

## ---- DMPs and their 450k annotation ----------------------------------------
dmp <- read.csv(file.path(M,"GSE74193/results/DMCpGs_SCZ_postnatal_autosomal_213.csv"),
                row.names=1, check.names=FALSE)
cat(sprintf("DMPs: %d | hypomethylated %d (%.1f%%)\n", nrow(dmp),
            sum(dmp$delta_beta < 0), 100*mean(dmp$delta_beta < 0)))

ann <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
a <- ann[rownames(dmp), c("chr","pos","UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_Island")]
dmp <- cbind(dmp, a)

## a probe may list several genes/regions, semicolon separated -> one row per pair
expand <- do.call(rbind, lapply(seq_len(nrow(dmp)), function(i){
  g <- strsplit(dmp$UCSC_RefGene_Name[i], ";")[[1]]
  r <- strsplit(dmp$UCSC_RefGene_Group[i], ";")[[1]]
  if(!length(g)) return(NULL)
  if(length(r) != length(g)) r <- rep(r[1], length(g))
  u <- !duplicated(paste(g, r))
  data.frame(cpg=rownames(dmp)[i], gene=g[u], region=r[u],
             delta_beta=dmp$delta_beta[i], logFC_meth=dmp$logFC[i],
             island=dmp$Relation_to_Island[i], stringsAsFactors=FALSE)
}))
cat(sprintf("annotated: %d of %d probes carry a gene | %d unique genes\n",
            length(unique(expand$cpg)), nrow(dmp), length(unique(expand$gene))))
cat("\nprobe positions relative to gene:\n"); print(sort(table(expand$region), decreasing=TRUE))

PROMOTER <- c("TSS200","TSS1500","5'UTR","1stExon")
expand$part <- ifelse(expand$region %in% PROMOTER, "promoter",
               ifelse(expand$region == "Body", "body", "other"))

## ---- expression ------------------------------------------------------------
ex <- read.csv(file.path(INT,"05_genelevel","SCZ_frontal.csv"), check.names=FALSE)
names(ex)[1] <- "ens"
ex <- ex[!is.na(ex$symbol) & !is.na(ex$logFC) & !duplicated(ex$symbol), ]
ex$stat <- sign(ex$logFC) * -log10(pmax(ex$P.Value, .Machine$double.xmin))
ex$abund <- ex$AveExpr                       # voom log-CPM: the abundance covariate
cat(sprintf("\nexpression: %d genes with a symbol\n", nrow(ex)))
ex$stratum <- cut(rank(ex$abund, ties.method="first"),
                  breaks=quantile(rank(ex$abund, ties.method="first"), probs=seq(0,1,0.05)),
                  include.lowest=TRUE, labels=FALSE)

test_set <- function(genes, label, predicted){
  idx <- which(ex$symbol %in% genes)
  if(length(idx) < 15){ cat(sprintf("\n%-28s %d genes -- too few to test\n", label, length(idx))); return(NULL) }
  mt <- median(ex$logFC[idx])
  tab  <- table(ex$stratum[idx])
  rest <- setdiff(seq_len(nrow(ex)), idx)
  pool <- split(rest, ex$stratum[rest])
  null <- replicate(NPERM, {
    s <- unlist(lapply(names(tab), function(k)
      if(!is.null(pool[[k]])) sample(pool[[k]], min(tab[[k]], length(pool[[k]]))) else integer(0)))
    median(ex$logFC[s]) })
  p <- (sum(abs(null - median(null)) >= abs(mt - median(null))) + 1)/(NPERM + 1)
  obs <- if(mt > median(null)) "UP" else "DOWN"
  cat(sprintf("\n--- %s: %d genes ---\n", label, length(idx)))
  cat(sprintf("  median logFC: set %+.4f | expression-matched null %+.4f (SD %.4f)\n",
              mt, median(null), sd(null)))
  cat(sprintf("  permutation p = %.4f | observed %s | predicted %s -> %s\n",
              p, obs, predicted, if(obs==predicted) "CONSISTENT" else "OPPOSITE"))
  cp <- do.call(rbind, lapply(c(0.01,0.05,0.10), function(r){
    x <- cameraPR(setNames(ex$stat, ex$symbol), list(S=idx), inter.gene.cor=r)
    data.frame(inter_gene_cor=r, Direction=x$Direction, PValue=x$PValue) }))
  print(cp, row.names=FALSE)
  data.frame(set=label, n=length(idx), median_logFC=mt, null_median=median(null),
             p_matched=p, observed=obs, predicted=predicted,
             cameraPR_igc0.05=cp$PValue[cp$inter_gene_cor==0.05])
}

res <- list()
## hypomethylated dominate, so promoter hypo -> predicted UP, body hypo -> predicted DOWN
res$prom_hypo <- test_set(unique(expand$gene[expand$part=="promoter" & expand$delta_beta<0]),
                          "promoter, hypomethylated", "UP")
res$prom_hyper<- test_set(unique(expand$gene[expand$part=="promoter" & expand$delta_beta>0]),
                          "promoter, hypermethylated", "DOWN")
res$body_hypo <- test_set(unique(expand$gene[expand$part=="body" & expand$delta_beta<0]),
                          "gene body, hypomethylated", "DOWN")
res$body_hyper<- test_set(unique(expand$gene[expand$part=="body" & expand$delta_beta>0]),
                          "gene body, hypermethylated", "UP")
res$all       <- test_set(unique(expand$gene), "all DMP genes (no direction)", "NA")

## per-gene table for the discussion
mg <- merge(expand, ex[,c("symbol","logFC","P.Value","adj.P.Val")],
            by.x="gene", by.y="symbol")
mg <- mg[order(mg$P.Value), ]
write.csv(mg, file.path(OUT,"SCZ_DMP_genes_with_expression.csv"), row.names=FALSE)
cat(sprintf("\n%d DMP gene-probe pairs matched to expression | top by expression p:\n", nrow(mg)))
print(head(mg[,c("gene","cpg","region","delta_beta","logFC","P.Value","adj.P.Val")], 12), row.names=FALSE)

S <- do.call(rbind, res)
if(!is.null(S)){ cat("\n================ summary ================\n"); print(S, row.names=FALSE)
  write.csv(S, file.path(OUT,"SCZ_methylation_axis_summary.csv"), row.names=FALSE) }
cat(sprintf("\nwritten to %s\n", OUT))
