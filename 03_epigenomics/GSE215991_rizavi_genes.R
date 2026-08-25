suppressPackageStartupMessages({library(csaw); library(edgeR)
  library(GenomicRanges); library(rtracklayer)})
S   <- "/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-ChIP-seq/GSE215991"
GTF <- "/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf"

e <- new.env(); load(file.path(S,"csaw/GSE215991_csaw_v2.RData"), envir=e)
filtered <- get("filtered", envir=e); y <- get("y", envir=e)
rr <- rowRanges(filtered)

g <- import(GTF, feature.type="gene")
seqlevels(g) <- paste0("chr", seqlevels(g))
g <- keepSeqlevels(g, intersect(seqlevels(g), seqlevels(rr)), pruning.mode="coarse")

# Rizavi et al. promoter definition: 10 kb centred on TSS
prom <- promoters(g, upstream=5000, downstream=5000)
mcols(prom)$sym <- mcols(g)$gene_name
cat("promoters (+/-5kb):", length(prom), "\n")

hit  <- findOverlaps(rr, prom)
keep <- unique(queryHits(hit))
cat("windows retained:", length(keep), "of", length(rr), "\n")

grp <- factor(colData(filtered)$group, levels=c("NPC","SCZ"))
des <- model.matrix(~ grp)
yp  <- estimateDisp(y[keep,], des)
res <- glmQLFTest(glmQLFit(yp, des, robust=TRUE), coef=2)

idx <- match(keep, queryHits(hit))
ids <- subjectHits(hit)[idx]
tab <- combineTests(ids, res$table)
sym <- mcols(prom)$sym[as.integer(names(table(ids)))]

cat("\n=== +/-5kb PROMOTERS, SCZ vs NPC ===\n")
cat("tested:", nrow(tab), " FDR<0.05:", sum(tab$FDR<0.05, na.rm=TRUE),
    " raw p<0.01:", sum(tab$PValue<0.01, na.rm=TRUE),
    " (chance", round(nrow(tab)*0.01), ")  [they report 159]\n")

tgt <- c("GPHN","CDHR2","GPX6","OXT","FURIN","NRGN","AKT3","CDH20","LIPJ","C4B")
cat("\n=== THEIR 10 SELECTED GENES IN YOUR DATA ===\n")
out <- do.call(rbind, lapply(tgt, function(s){
  i <- which(sym == s)
  if(!length(i)) data.frame(gene=s, logFC=NA, PValue=NA, FDR=NA, dir=NA)
  else data.frame(gene=s, logFC=round(tab$rep.logFC[i[1]],3),
                  PValue=signif(tab$PValue[i[1]],3),
                  FDR=signif(tab$FDR[i[1]],3), dir=tab$direction[i[1]])
}))
print(out, row.names=FALSE)
saveRDS(list(tab=tab, sym=sym), file.path(S,"csaw","rizavi_promoter_test.rds"))
