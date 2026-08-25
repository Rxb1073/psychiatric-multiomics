.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)})
logf <- function(...) cat(sprintf(...),"\n")
D <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/counts_scz/HIPP_final"
r <- read.csv(file.path(D,"hipp_R1_full.csv"))
ecm <- read.csv(file.path(D,"ORA_GOCC_DOWN.csv"))
row <- ecm[grepl("^extracellular matrix$", ecm$Description),][1,]
gs  <- strsplit(as.character(row$geneID), "/")[[1]]
sym <- AnnotationDbi::mapIds(org.Hs.eg.db, gs, "SYMBOL","ENTREZID", multiVals="first")
ens <- AnnotationDbi::mapIds(org.Hs.eg.db, gs, "ENSEMBL","ENTREZID", multiVals="first")
logf("ECM term genes: %d  ->  %s", length(gs), paste(na.omit(sym), collapse=", "))
r$base <- sub("\\..*","",r$gene)
sub <- r[r$base %in% na.omit(ens),]
logf("matched in R1 results: %d", nrow(sub))
print(sub[order(sub$padj), c("sym","log2FoldChange","padj")])
write.csv(sub, file.path(D,"ECM_genes_R1.csv"), row.names=FALSE)
