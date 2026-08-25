## Rebuild the host KEGG GSEA WITHOUT the significance filter, so the cross-omics
## correlation can run on all tested pathways rather than the significant subset.
##
## The archived tables (GSEA_KEGG_primary.csv: 70 rows SCZ, 126 MDD) were written
## with pvalueCutoff=0.05, and the script that produced them is not on disk. This
## rebuilds the ranked list from the gene-level tables per the chapter's methods
## -- sign(log2FC) x -log10(p), gene sets 15-500, BH -- and runs gseKEGG with
## pvalueCutoff=1.
##
## VERIFICATION IS THE POINT. The significant subset of the new output is compared
## against the archived table. If IDs and NES reproduce, the ranking was rebuilt
## faithfully and the full list can be trusted. If not, the archived tables were
## built from something other than what the methods describe, and that must be
## resolved before either version is used.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)})
set.seed(1)

D   <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final"
INT <- file.path(D, "INTEGRATION_inputs")
OUT <- file.path(INT, "06_kegg_axis", "host_unfiltered"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

cat("clusterProfiler", as.character(packageVersion("clusterProfiler")),
    "| org.Hs.eg.db", as.character(packageVersion("org.Hs.eg.db")), "\n")

pick <- function(df, cands, what){
  hit <- cands[cands %in% names(df)]
  if(!length(hit)) stop(sprintf("no %s column found; available: %s", what, paste(names(df), collapse=", ")))
  if(length(hit)>1) cat(sprintf("  NOTE: several %s columns (%s); using %s\n", what, paste(hit,collapse=", "), hit[1]))
  hit[1]
}

build_rank <- function(df, arm){
  pcol <- pick(df, c("P.Value","pvalue","p_FE","p_fe","PValue","pval","q_fe"), "p-value")
  lcol <- pick(df, c("logFC","log2FoldChange","beta_fe","beta_FE","lfc"), "effect size")
  scol <- if("symbol" %in% names(df)) "symbol" else if("sym" %in% names(df)) "sym" else NA
  cat(sprintf("  %s: ranking on sign(%s) x -log10(%s)\n", arm, lcol, pcol))
  if(pcol %in% c("q_fe")) cat("  *** WARNING: only an ADJUSTED p-value is available in this file.\n",
                              "      The ranking is therefore not exactly as the methods describe.\n")
  d <- data.frame(sym = if(is.na(scol)) NA else df[[scol]],
                  lfc = suppressWarnings(as.numeric(df[[lcol]])),
                  p   = suppressWarnings(as.numeric(df[[pcol]])),
                  stringsAsFactors = FALSE)
  ## gene identifier: symbol if present, else the row name / first column
  if(is.na(scol)){
    id <- if("gene" %in% names(df)) df$gene else rownames(df)
    ens <- sub("\\..*$","", id)
    d$sym <- AnnotationDbi::mapIds(org.Hs.eg.db, ens, "SYMBOL","ENSEMBL", multiVals="first")
  }
  d <- d[!is.na(d$sym) & !is.na(d$lfc) & !is.na(d$p), ]
  d$p[d$p <= 0] <- .Machine$double.xmin
  d$stat <- sign(d$lfc) * -log10(d$p)
  d$entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, d$sym, "ENTREZID","SYMBOL", multiVals="first")
  d <- d[!is.na(d$entrez) & !duplicated(d$entrez), ]
  cat(sprintf("  %s: %d genes with Entrez ID\n", arm, nrow(d)))
  gl <- setNames(d$stat, d$entrez)
  sort(gl, decreasing=TRUE)
}

run_arm <- function(arm, path, archived){
  cat(sprintf("\n================ %s ================\n", arm))
  df <- if(grepl("\\.tsv$", path)) read.delim(path, check.names=FALSE)
        else read.csv(path, check.names=FALSE, row.names=1)
  cat(sprintf("  input: %s (%d rows)\n  columns: %s\n", basename(path), nrow(df),
              paste(names(df), collapse=", ")))
  gl <- build_rank(df, arm)

  k <- gseKEGG(geneList=gl, organism="hsa", minGSSize=15, maxGSSize=500,
               pAdjustMethod="BH", pvalueCutoff=1, verbose=FALSE, seed=TRUE)
  r <- as.data.frame(k)
  cat(sprintf("\n  UNFILTERED: %d pathways tested | %d at p.adjust<0.05\n",
              nrow(r), sum(r$p.adjust < 0.05)))
  write.csv(r, file.path(OUT, sprintf("GSEA_KEGG_%s_unfiltered.csv", arm)), row.names=FALSE)

  ## ---- reproduction check against the archived table ----
  if(file.exists(archived)){
    old <- read.csv(archived, check.names=FALSE)
    new_sig <- r[r$p.adjust < 0.05, ]
    cat(sprintf("\n  --- reproduction check ---\n  archived: %d rows | new significant: %d rows\n",
                nrow(old), nrow(new_sig)))
    shared <- intersect(old$ID, new_sig$ID)
    cat(sprintf("  IDs in common: %d (%.0f%% of archived)\n", length(shared), 100*length(shared)/nrow(old)))
    if(length(shared)){
      m <- merge(old[,c("ID","NES")], new_sig[,c("ID","NES")], by="ID", suffixes=c("_old","_new"))
      cat(sprintf("  NES correlation on shared IDs: r = %.4f | max abs difference = %.4f\n",
                  cor(m$NES_old, m$NES_new), max(abs(m$NES_old - m$NES_new))))
      if(cor(m$NES_old, m$NES_new) > 0.999 && max(abs(m$NES_old-m$NES_new)) < 0.01)
        cat("  => REPRODUCED. The unfiltered list is faithful to the archived analysis.\n")
      else
        cat("  => DOES NOT REPRODUCE. The archived table came from a different ranking;\n",
            "     resolve before using either version.\n")
    }
    missing <- setdiff(old$ID, new_sig$ID)
    if(length(missing)) cat(sprintf("  archived IDs absent from new significant set: %d (%s%s)\n",
        length(missing), paste(head(missing,5), collapse=", "), if(length(missing)>5) ", ..." else ""))
  } else cat(sprintf("\n  (no archived table at %s to check against)\n", archived))
  invisible(r)
}

r_scz <- run_arm("SCZ", file.path(INT,"05_genelevel","SCZ_frontal.csv"),
                 file.path(D,"Transcriptomics_SCZ/enrichment/GSEA_KEGG_primary.csv"))
r_mdd <- run_arm("MDD", file.path(INT,"05_genelevel","MDD_meta.tsv"),
                 file.path(D,"Transcriptomics_MDD/enrichment/GSEA_KEGG_primary.csv"))

cat(sprintf("\nwritten to %s\n", OUT))
