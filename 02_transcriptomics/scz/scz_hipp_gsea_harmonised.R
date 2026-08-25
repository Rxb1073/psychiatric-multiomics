## GSEA harmonised to ASD ranking: sign(LFC) * -log10(p), GO BP/MF/CC + KEGG
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)})
logf <- function(...) cat(sprintf(...),"\n")
D <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/counts_scz/HIPP_final"
r <- read.csv(file.path(D,"hipp_R1_full.csv"))
r <- r[!is.na(r$pvalue) & !is.na(r$log2FoldChange), ]
r$entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, sub("\\..*","",r$gene),
                                  "ENTREZID","ENSEMBL", multiVals="first")
r <- r[!is.na(r$entrez) & !duplicated(r$entrez), ]
r$rankstat <- sign(r$log2FoldChange) * -log10(pmax(r$pvalue, .Machine$double.xmin))
gl <- sort(setNames(r$rankstat, r$entrez), decreasing=TRUE)
logf("ranked list: %d genes (ASD-harmonised rankstat)", length(gl))
for (ont in c("BP","MF","CC")) {
  g <- tryCatch(gseGO(geneList=gl, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont=ont,
                pAdjustMethod="BH", pvalueCutoff=0.05, verbose=FALSE), error=function(e) NULL)
  n <- if (is.null(g)) 0 else nrow(g@result)
  logf("GSEA GO-%s : %d terms", ont, n)
  if (n) {
    write.csv(g@result, file.path(D, sprintf("GSEA_harm_GO%s.csv", ont)), row.names=FALSE)
    x <- g@result[order(g@result$p.adjust), c("Description","NES","p.adjust")]
    print(head(x, 10))
  }
}
k <- tryCatch(gseKEGG(geneList=gl, organism="hsa", pAdjustMethod="BH",
              pvalueCutoff=0.05, verbose=FALSE), error=function(e) NULL)
n <- if (is.null(k)) 0 else nrow(k@result)
logf("GSEA KEGG : %d terms", n)
if (n) {
  write.csv(k@result, file.path(D,"GSEA_harm_KEGG.csv"), row.names=FALSE)
  print(head(k@result[order(k@result$p.adjust), c("Description","NES","p.adjust")], 10))
}
