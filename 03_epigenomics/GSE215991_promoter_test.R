suppressPackageStartupMessages({library(csaw); library(edgeR)
  library(GenomicRanges); library(rtracklayer)})
S   <- "/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-ChIP-seq/GSE215991"
GTF <- "/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf"

e <- new.env(); load(file.path(S,"csaw/GSE215991_csaw_v2.RData"), envir=e)
filtered <- get("filtered", envir=e); y <- get("y", envir=e)
rr <- rowRanges(filtered)

gtf <- import(GTF, feature.type="gene")
seqlevels(gtf) <- paste0("chr", seqlevels(gtf))          # Ensembl -> UCSC
gtf <- keepSeqlevels(gtf, intersect(seqlevels(gtf), seqlevels(rr)), pruning.mode="coarse")
cat("genes:", length(gtf), " shared seqlevels:", length(seqlevels(gtf)), "\n")

prom <- reduce(promoters(gtf, upstream=2000, downstream=2000))
cat("promoter regions:", length(prom), "\n")

hit  <- findOverlaps(rr, prom)
keep <- unique(queryHits(hit))
cat("windows retained:", length(keep), "of", length(rr), "\n")
stopifnot(length(keep) > 1000)

ids <- subjectHits(hit)[match(keep, queryHits(hit))]
yp  <- y[keep, ]
grp <- factor(colData(filtered)$group, levels=c("NPC","SCZ"))
des <- model.matrix(~ grp)
yp  <- estimateDisp(yp, des)
res <- glmQLFTest(glmQLFit(yp, des, robust=TRUE), coef=2)
tab <- combineTests(ids, res$table)

cat("\n=== PROMOTER-RESTRICTED, SCZ vs NPC ===\n")
cat("promoters tested:", nrow(tab), " FDR<0.05:", sum(tab$FDR<0.05, na.rm=TRUE), "\n")
cat("nominal p<0.05:", sum(tab$PValue<0.05, na.rm=TRUE),
    " (chance", round(nrow(tab)*0.05), ")\n")
print(head(tab[order(tab$PValue), ], 10))

saveRDS(list(tab=tab, prom=prom[sort(unique(ids))]),
        file.path(S,"csaw","promoter_test_GSE215991.rds"))
