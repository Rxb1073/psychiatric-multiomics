## ASD cross-omics link, v2: CHD8 binding targets against the ASD meta-analysis.
##
## Three changes from v1, each fixing something that run exposed.
##
## 1. PEAK RANKING. All 5,507 consensus sites annotate to ~6,300 promoters, a
##    quarter of the tested transcriptome. A competitive test on a set that large
##    is degenerate: the "rest of the transcriptome" it compares against contains
##    the same signal, which is why v1's p-value ran from 4e-58 to 0.42 depending
##    on an assumed inter-gene correlation. Sites are therefore ranked by MACS2
##    fold enrichment -- taking the MINIMUM across CHD8 antibody preparations, so
##    a site scores highly only if every preparation supports it -- and tested at
##    several set sizes.
##
## 2. ONE GENE PER PEAK. v1 took every promoter overlapping a peak, so gene count
##    exceeded peak count. Each peak is now assigned to its single nearest TSS.
##
## 3. EXPRESSION-MATCHED NULL. CHD8 binds active promoters, so its targets are
##    more highly expressed than average, and expression level is itself related
##    to effect size in a meta-analysis. Comparing targets against ALL other genes
##    therefore confounds binding with abundance. The primary test compares the
##    observed median effect against null sets drawn to match the targets'
##    baseMean distribution, so the question becomes: do CHD8 targets shift more
##    than equally-expressed non-targets?
##
## Directional prediction: CHD8 haploinsufficiency lowers expression of bound
## genes, so targets are predicted to shift DOWN. v1 found them shifted UP.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(limma); library(GenomicRanges)})
set.seed(1)

EPI <- "/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-CHD8/GSE61492/peaks_corrected"
GTF <- "/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf"
INT <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final/INTEGRATION_inputs"
OUT <- file.path(INT, "07_chd8_axis"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
FLANK <- 2000; NPERM <- 2000

## ---- TSS from the shared Ensembl GTF ---------------------------------------
cat("reading gene records from the shared GTF...\n")
gt <- read.delim(GTF, header=FALSE, comment.char="#", quote="",
                 col.names=c("chr","src","feature","start","end","score","strand","frame","attr"),
                 colClasses=c("character","NULL","character","integer","integer",
                              "NULL","character","NULL","character"))
gt <- gt[gt$feature=="gene", ]
gt$ens <- sub('.*gene_id "([^"]+)".*', "\\1", gt$attr)
gt$bt  <- sub('.*gene_biotype "([^"]+)".*', "\\1", gt$attr)
gt <- gt[gt$bt %in% c("protein_coding","lncRNA"), ]
gt$tss <- ifelse(gt$strand=="+", gt$start, gt$end)
tss <- GRanges(paste0("chr", gt$chr), IRanges(gt$tss, gt$tss), ens=gt$ens)
cat(sprintf("  %d protein-coding/lncRNA TSSs\n", length(tss)))

## ---- consensus peaks, scored by minimum fold enrichment --------------------
bed <- read.delim(file.path(EPI,"CHD8_specific_consensus.bed"), header=FALSE)[,1:3]
names(bed) <- c("chr","start","end")
pk <- GRanges(bed$chr, IRanges(bed$start+1, bed$end))
cat(sprintf("\nconsensus: %d sites, mean width %.0f bp\n", length(pk), mean(width(pk))))

npf <- list.files(EPI, pattern="^GSE61492_CHD8_.*_peaks\\.narrowPeak$", full.names=TRUE)
cat(sprintf("scoring against %d CHD8 preparations:\n  %s\n", length(npf),
            paste(basename(npf), collapse="\n  ")))
sc <- matrix(NA_real_, nrow=length(pk), ncol=length(npf),
             dimnames=list(NULL, basename(npf)))
for(j in seq_along(npf)){
  np <- read.delim(npf[j], header=FALSE)
  g  <- GRanges(np$V1, IRanges(np$V2+1, np$V3), fe=np$V7)     # V7 = fold enrichment
  ov <- findOverlaps(pk, g)
  best <- tapply(g$fe[subjectHits(ov)], queryHits(ov), max)
  sc[as.integer(names(best)), j] <- best
}
cat(sprintf("  sites with a score in all preparations: %d of %d\n",
            sum(complete.cases(sc)), nrow(sc)))
## minimum across preparations: conservative, rewards agreement
bed$score <- apply(sc, 1, function(r) if(all(is.na(r))) NA else min(r, na.rm=TRUE))
cat(sprintf("  fold enrichment (min across preps): median %.2f, range %.2f-%.2f\n",
            median(bed$score, na.rm=TRUE), min(bed$score, na.rm=TRUE), max(bed$score, na.rm=TRUE)))

## ---- one gene per peak: nearest TSS within FLANK ---------------------------
nr <- distanceToNearest(pk, tss)
ok <- !is.na(mcols(nr)$distance) & mcols(nr)$distance <= FLANK
peak2gene <- data.frame(peak = queryHits(nr)[ok],
                        ens  = tss$ens[subjectHits(nr)[ok]],
                        dist = mcols(nr)$distance[ok],
                        score= bed$score[queryHits(nr)[ok]],
                        stringsAsFactors=FALSE)
cat(sprintf("\n%d of %d peaks within %d bp of a TSS (%.1f%%) -> %d unique genes\n",
            nrow(peak2gene), length(pk), FLANK, 100*nrow(peak2gene)/length(pk),
            length(unique(peak2gene$ens))))
## one row per gene, keeping its best-scoring peak
peak2gene <- peak2gene[order(-peak2gene$score), ]
peak2gene <- peak2gene[!duplicated(peak2gene$ens), ]

## ---- ASD ranked list -------------------------------------------------------
asd <- read.csv(file.path(INT,"05_genelevel","ASD_meta.csv"), check.names=FALSE)
asd$ens <- sub("\\..*$", "", asd$ens)
asd <- asd[!is.na(asd$beta_FE) & !is.na(asd$p_FE) & !duplicated(asd$ens), ]
asd$p_FE[asd$p_FE <= 0] <- .Machine$double.xmin
asd$stat <- sign(asd$beta_FE) * -log10(asd$p_FE)
asd$baseMean <- rowMeans(asd[, c("baseMean_1","baseMean_2")], na.rm=TRUE)
cat(sprintf("ASD meta: %d genes usable\n", nrow(asd)))

## expression strata for the matched null
asd$stratum <- cut(rank(asd$baseMean, ties.method="first"),
                   breaks=quantile(rank(asd$baseMean, ties.method="first"), probs=seq(0,1,0.05)),
                   include.lowest=TRUE, labels=FALSE)

test_set <- function(genes, label){
  idx <- which(asd$ens %in% genes)
  if(length(idx) < 20){ cat(sprintf("\n%s: only %d genes, skipped\n", label, length(idx))); return(NULL) }
  mt <- median(asd$beta_FE[idx]); mo <- median(asd$beta_FE[-idx])

  ## abundance check: are targets more highly expressed?
  em <- median(asd$baseMean[idx]); eo <- median(asd$baseMean[-idx])

  ## expression-matched null: resample non-targets to match stratum composition
  tab <- table(asd$stratum[idx])
  pool <- split(setdiff(seq_len(nrow(asd)), idx), asd$stratum[setdiff(seq_len(nrow(asd)), idx)])
  null <- replicate(NPERM, {
    s <- unlist(lapply(names(tab), function(k)
      if(!is.null(pool[[k]])) sample(pool[[k]], min(tab[[k]], length(pool[[k]]))) else integer(0)))
    median(asd$beta_FE[s]) })
  p_matched <- (sum(abs(null - median(null)) >= abs(mt - median(null))) + 1) / (NPERM + 1)

  cat(sprintf("\n--- %s: %d genes ---\n", label, length(idx)))
  cat(sprintf("  median beta_FE  targets %+.4f | all non-targets %+.4f\n", mt, mo))
  cat(sprintf("  median baseMean targets %8.1f | all non-targets %8.1f  (%.2fx)\n", em, eo, em/eo))
  cat(sprintf("  expression-matched null: median %+.4f (SD %.4f) -> p = %.4f\n",
              median(null), sd(null), p_matched))
  cat(sprintf("  direction: %s\n", if(mt < median(null)) "DOWN vs matched null" else "UP vs matched null"))
  cp <- do.call(rbind, lapply(c(0.01, 0.05, 0.10), function(r){
    x <- cameraPR(setNames(asd$stat, asd$ens), list(S=idx), inter.gene.cor=r)
    data.frame(inter_gene_cor=r, Direction=x$Direction, PValue=x$PValue) }))
  print(cp, row.names=FALSE)
  data.frame(set=label, n=length(idx), median_target=mt, median_matched_null=median(null),
             p_matched=p_matched, abundance_ratio=em/eo,
             cameraPR_p_igc0.05=cp$PValue[cp$inter_gene_cor==0.05])
}

res <- list()
for(N in c(1000, 2000, 3000, nrow(peak2gene))){
  g <- head(peak2gene$ens[order(-peak2gene$score)], N)
  res[[as.character(N)]] <- test_set(g, sprintf("top %d peaks by fold enrichment", N))
}
S <- do.call(rbind, res)
cat("\n================ summary ================\n"); print(S, row.names=FALSE)
write.csv(S, file.path(OUT,"CHD8_ranked_summary.csv"), row.names=FALSE)
write.csv(peak2gene, file.path(OUT,"CHD8_peak_to_gene.csv"), row.names=FALSE)
cat(sprintf("\nwritten to %s\n", OUT))
