#!/usr/bin/env Rscript
# GO/KEGG enrichment on the ASD childhood DE results.
#  (a) ORA  on the strict gene lists (FDR<0.05 & |LFC|>=1)
#  (b) GSEA on the full ranked gene list (better powered, uses all genes)
# Universe = genes actually tested by DESeq2 (NOT all annotated genes).
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
ok <- suppressMessages(requireNamespace("clusterProfiler", quietly=TRUE))
if(!ok) stop("clusterProfiler not installed. Run:\n  install.packages('BiocManager'); BiocManager::install(c('clusterProfiler','org.Hs.eg.db','enrichplot'))")
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db); library(ggplot2)})
logf <- function(...) cat(sprintf(...), "\n")

BASE <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
IN   <- file.path(BASE, "POOLED_cortex")
OUT  <- file.path(BASE, "DISSERTATION_final", "enrichment"); dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
FDR <- 0.05; LFC <- 1

analyses <- list(
  pooled_3dataset = "DE_ASD_child_3dataset_full.csv",
  GSE59288only    = "DE_ASD_child_GSE59288only_full.csv")

run_one <- function(tag, fn){
  f <- file.path(IN, fn); if(!file.exists(f)){logf("MISSING %s - skip", f); return(invisible())}
  r <- read.csv(f, stringsAsFactors=FALSE)
  r <- r[!is.na(r$padj) & !is.na(r$ens), ]

  # ---- map to ENTREZ (clusterProfiler's native id) ----
  r$entrez <- suppressMessages(mapIds(org.Hs.eg.db, r$ens, "ENTREZID", "ENSEMBL"))
  r <- r[!is.na(r$entrez), ]
  r <- r[!duplicated(r$entrez), ]

  universe <- r$entrez                                    # all tested genes
  strict   <- r$entrez[r$padj < FDR & abs(r$log2FoldChange) >= LFC]
  down     <- r$entrez[r$padj < FDR & r$log2FoldChange <= -LFC]
  up       <- r$entrez[r$padj < FDR & r$log2FoldChange >=  LFC]

  logf("\n===== %s =====", tag)
  logf("universe=%d | strict=%d (down=%d up=%d)", length(universe), length(strict), length(down), length(up))

  # ---------- (a) ORA ----------
  sets <- list(strict=strict, down=down, up=up)
  for(sn in names(sets)){
    g <- sets[[sn]]; if(length(g) < 5){logf("  [%s] only %d genes - too few for ORA, skipped", sn, length(g)); next}
    for(ont in c("BP","MF","CC")){
      e <- tryCatch(enrichGO(gene=g, universe=universe, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
                             ont=ont, pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.10,
                             readable=TRUE), error=function(e) NULL)
      n <- if(is.null(e)) 0 else nrow(as.data.frame(e))
      logf("  ORA GO:%s [%s] -> %d significant terms", ont, sn, n)
      if(n > 0){
        df <- as.data.frame(e)
        write.csv(df, file.path(OUT, sprintf("ORA_GO%s_%s_%s.csv", ont, sn, tag)), row.names=FALSE)
        cat(sprintf("      top: %s\n", paste(head(df$Description,5), collapse=" | ")))
        if(n >= 3){
          p <- dotplot(e, showCategory=min(15,n)) + ggtitle(sprintf("%s - GO:%s (%s)", tag, ont, sn)) +
               theme(axis.text.y=element_text(size=8))
          ggsave(file.path(OUT, sprintf("ORA_GO%s_%s_%s.png", ont, sn, tag)), p,
                 width=8, height=6.5, dpi=200, type="cairo")
        }
      }
    }
    k <- tryCatch(enrichKEGG(gene=g, universe=universe, organism="hsa",
                             pAdjustMethod="BH", pvalueCutoff=0.05), error=function(e) NULL)
    nk <- if(is.null(k)) 0 else nrow(as.data.frame(k))
    logf("  ORA KEGG [%s] -> %d significant pathways", sn, nk)
    if(nk > 0) write.csv(as.data.frame(k), file.path(OUT, sprintf("ORA_KEGG_%s_%s.csv", sn, tag)), row.names=FALSE)
  }

  # ---------- (b) GSEA on full ranked list ----------
  # rank by signed -log10(p): keeps direction, uses every tested gene
  r$rankstat <- sign(r$log2FoldChange) * -log10(pmax(r$pvalue, .Machine$double.xmin))
  gl <- sort(setNames(r$rankstat, r$entrez), decreasing=TRUE)
  gl <- gl[is.finite(gl)]
  for(ont in c("BP","CC")){
    g <- tryCatch(gseGO(geneList=gl, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont=ont,
                        pAdjustMethod="BH", pvalueCutoff=0.05, minGSSize=15, maxGSSize=500,
                        verbose=FALSE, seed=TRUE), error=function(e) NULL)
    n <- if(is.null(g)) 0 else nrow(as.data.frame(g))
    logf("  GSEA GO:%s -> %d significant terms", ont, n)
    if(n > 0){
      df <- as.data.frame(g)
      df$direction <- ifelse(df$NES > 0, "up in ASD", "down in ASD")
      write.csv(df, file.path(OUT, sprintf("GSEA_GO%s_%s.csv", ont, tag)), row.names=FALSE)
      cat(sprintf("      top: %s\n", paste(head(df$Description,5), collapse=" | ")))
      if(n >= 3){
        p <- dotplot(g, showCategory=min(15,n), split=".sign") + facet_grid(.~.sign) +
             ggtitle(sprintf("%s - GSEA GO:%s", tag, ont)) + theme(axis.text.y=element_text(size=8))
        ggsave(file.path(OUT, sprintf("GSEA_GO%s_%s.png", ont, tag)), p,
               width=9.5, height=6.5, dpi=200, type="cairo")
      }
    }
  }
}

for(tag in names(analyses)) run_one(tag, analyses[[tag]])
logf("\nAll enrichment output -> %s", OUT)
