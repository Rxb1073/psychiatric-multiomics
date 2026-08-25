## Third integration axis: microbial pathways against epigenomically-implicated genes.
##
## The two axes already reported both use the transcriptome as the common layer.
## This one does not, so the layers must be brought into a shared feature space by
## mapping the epigenomic gene set onto KEGG pathways and asking whether the
## pathways it implicates overlap those altered in the gut.
##
## WHAT THIS CAN AND CANNOT SHOW. A shared KEGG identifier means a bacterial map
## and a human map carry the same five-digit number, not that the same enzymes are
## involved. Agreement would indicate that both layers implicate overlapping
## metabolic categories; it would not indicate a shared molecular process. The test
## is also weaker than the gene-level tests already reported, because pathway
## membership is a coarser unit than a ranked gene list.
##
## POWER. Autism has 4,351 CHD8 target genes against 14 microbial pathways.
## Schizophrenia has roughly 160 methylation genes against 38 pathways, which will
## map to few pathways and is expected to be underpowered; it is run and reported
## regardless, so the null is stated rather than the arm quietly omitted.
## Depression retained no microbial pathway after class filtering, so no test.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)
                  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)})
set.seed(1)

D   <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final"
M   <- "/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-methylation"
INT <- file.path(D, "INTEGRATION_inputs")
AX  <- file.path(INT, "06_kegg_axis")
OUT <- file.path(INT, "14_microbiome_epigenomics"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

## ---- epigenomic gene sets --------------------------------------------------
## Autism: CHD8 target genes (Ensembl), converted to Entrez for KEGG
p2g  <- read.csv(file.path(INT,"07_chd8_axis","CHD8_peak_to_gene.csv"), stringsAsFactors=FALSE)
asd_ens <- unique(sub("\\..*$","",p2g$ens))
asd_ent <- na.omit(unique(AnnotationDbi::mapIds(org.Hs.eg.db, asd_ens, "ENTREZID","ENSEMBL", multiVals="first")))
cat(sprintf("ASD  CHD8 targets: %d genes -> %d with Entrez ID\n", length(asd_ens), length(asd_ent)))

## Schizophrenia: genes annotated to the 213 differentially methylated positions
dmp <- read.csv(file.path(M,"GSE74193/results/DMCpGs_SCZ_postnatal_autosomal_213.csv"),
                row.names=1, check.names=FALSE)
ann <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
sym <- unique(unlist(strsplit(ann[rownames(dmp),"UCSC_RefGene_Name"], ";")))
sym <- sym[sym != "" & !is.na(sym)]
scz_ent <- na.omit(unique(AnnotationDbi::mapIds(org.Hs.eg.db, sym, "ENTREZID","SYMBOL", multiVals="first")))
cat(sprintf("SCZ  DMP genes:     %d symbols -> %d with Entrez ID\n", length(sym), length(scz_ent)))

## ---- map each epigenomic set onto human KEGG pathways ----------------------
## Over-representation, not GSEA: these are unordered sets, not ranked lists.
kegg_of <- function(entrez, label){
  e <- enrichKEGG(gene=entrez, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
  if(is.null(e) || !nrow(as.data.frame(e))){ cat(sprintf("  %s: no KEGG result\n", label)); return(NULL) }
  r <- as.data.frame(e); r$suffix <- sub("^hsa","",r$ID)
  cat(sprintf("  %s: %d KEGG pathways tested, %d at q<0.05\n", label, nrow(r), sum(r$qvalue<0.05, na.rm=TRUE)))
  r
}
cat("\nmapping epigenomic sets to human KEGG pathways\n")
k_asd <- kegg_of(asd_ent, "ASD CHD8 targets")
k_scz <- kegg_of(scz_ent, "SCZ DMP genes")

## ---- microbial pathways ----------------------------------------------------
mb <- function(a){
  f <- file.path(AX, sprintf("KO_pathway_DA_%s_filtered.tsv", a))
  d <- read.delim(f, check.names=FALSE)
  d$suffix <- sub("^ko","", d$feature)
  d
}
res <- list()
run <- function(arm, kepi){
  cat(sprintf("\n================ %s ================\n", arm))
  if(is.null(kepi)){ cat("  no epigenomic KEGG result; test not available\n"); return(NULL) }
  m <- mb(arm)
  pool <- intersect(kepi$suffix, m$suffix)      # pathways testable in BOTH layers
  epi_sig <- intersect(kepi$suffix[kepi$qvalue < 0.05], pool)
  mic_sig <- intersect(m$suffix[m$qval_refiltered < 0.05], pool)
  ov <- intersect(epi_sig, mic_sig)
  cat(sprintf("  pathways testable in both layers : %d\n", length(pool)))
  cat(sprintf("  significant, epigenomic side     : %d\n", length(epi_sig)))
  cat(sprintf("  significant, microbial side      : %d\n", length(mic_sig)))
  cat(sprintf("  overlap                          : %d\n", length(ov)))
  if(!length(pool) || !length(epi_sig) || !length(mic_sig)){
    cat("  one side has no significant pathway in the shared pool; no test possible\n")
    return(data.frame(arm=arm, pool=length(pool), epi=length(epi_sig), mic=length(mic_sig),
                      overlap=length(ov), p=NA_real_, fold=NA_real_))
  }
  p <- phyper(length(ov)-1, length(epi_sig), length(pool)-length(epi_sig), length(mic_sig), lower.tail=FALSE)
  fold <- (length(ov)/length(mic_sig)) / (length(epi_sig)/length(pool))
  cat(sprintf("  hypergeometric P = %.4f | %.1f-fold\n", p, fold))
  if(length(ov)){
    nm <- kepi$Description[match(ov, kepi$suffix)]
    cat("  shared pathways:\n"); for(i in seq_along(ov)) cat(sprintf("    ko%s / hsa%s  %s\n", ov[i], ov[i], nm[i]))
  }
  write.csv(data.frame(suffix=pool,
                       epi_q = kepi$qvalue[match(pool, kepi$suffix)],
                       mic_q = m$qval_refiltered[match(pool, m$suffix)],
                       description = kepi$Description[match(pool, kepi$suffix)]),
            file.path(OUT, sprintf("microbiome_epigenomics_%s.csv", arm)), row.names=FALSE)
  data.frame(arm=arm, pool=length(pool), epi=length(epi_sig), mic=length(mic_sig),
             overlap=length(ov), p=p, fold=fold)
}
res$ASD <- run("ASD", k_asd)
res$SCZ <- run("SCZ", k_scz)
cat("\n================ MDD ================\n")
cat("  no microbial pathway survived class-restricted correction; no test available\n")

S <- do.call(rbind, res)
if(!is.null(S)){ cat("\n================ SUMMARY ================\n"); print(S, row.names=FALSE)
  write.csv(S, file.path(OUT,"microbiome_epigenomics_summary.csv"), row.names=FALSE) }
cat(sprintf("\nwritten to %s\n", OUT))
