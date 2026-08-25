suppressPackageStartupMessages({library(csaw); library(edgeR); library(GenomicRanges)})
S <- "/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-ChIP-seq/GSE215991"
e <- new.env(); load(file.path(S, "csaw/GSE215991_csaw_v2.RData"), envir = e)
filtered <- get("filtered", envir=e); y <- get("y", envir=e); merged <- get("merged", envir=e)

srr <- sub("\\.filtered\\.bam$", "", basename(colData(filtered)$bam.files))
fem <- c("SRR21948649","SRR21948657","SRR21948663","SRR21948666","SRR21948673","SRR21948675")
sex <- factor(ifelse(srr %in% fem, "F", "M"), levels=c("M","F"))
print(table(sex, colData(filtered)$group))

des <- model.matrix(~ sex)
y2  <- estimateDisp(y, des)
res <- glmQLFTest(glmQLFit(y2, des, robust=TRUE), coef=2)
tab <- combineTests(merged$ids, res$table)

sig <- which(tab$FDR < 0.05)
cat("\n=== SEX CONTRAST ===\nregions:", nrow(tab), " FDR<0.05:", length(sig), "\n")
if (length(sig)) {
  ch <- as.character(seqnames(merged$regions)[sig])
  cat("on chrX/chrY:", sum(ch %in% c("chrX","chrY")),
      sprintf(" (%.1f%%)\n", 100*mean(ch %in% c("chrX","chrY"))))
  o <- order(tab$FDR[sig])
  print(head(data.frame(chr=ch, start=start(merged$regions)[sig],
                        FDR=signif(tab$FDR[sig],3), dir=tab$direction[sig])[o, ], 10))
}
saveRDS(list(tab=tab, regions=merged$regions, sex=sex),
        file.path(S, "csaw", "sexcontrol_GSE215991.rds"))
