.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)
                  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)})
set.seed(1)
D   <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final"
M   <- "/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-methylation"
INT <- file.path(D, "INTEGRATION_inputs"); AX <- file.path(INT, "06_kegg_axis")
OUT <- file.path(INT, "14_microbiome_epigenomics"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

p2g <- read.csv(file.path(INT,"07_chd8_axis","CHD8_peak_to_gene.csv"), stringsAsFactors=FALSE)
asd_ens <- unique(sub("\\..*$","",p2g$ens))
asd_ent <- na.omit(unique(AnnotationDbi::mapIds(org.Hs.eg.db, asd_ens, "ENTREZID","ENSEMBL", multiVals="first")))
cat(sprintf("ASD CHD8 targets: %d genes -> %d Entrez\n", length(asd_ens), length(asd_ent)))

dmp <- read.csv(file.path(M,"GSE74193/results/DMCpGs_SCZ_postnatal_autosomal_213.csv"), row.names=1, check.names=FALSE)
ann <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
sym <- unique(unlist(strsplit(ann[rownames(dmp),"UCSC_RefGene_Name"], ";")))
sym <- sym[sym != "" & !is.na(sym)]
scz_ent <- na.omit(unique(AnnotationDbi::mapIds(org.Hs.eg.db, sym, "ENTREZID","SYMBOL", multiVals="first")))
cat(sprintf("SCZ DMP genes: %d symbols -> %d Entrez\n\n", length(sym), length(scz_ent)))

kegg_of <- function(entrez, label){
  e <- enrichKEGG(gene=entrez, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
  if(is.null(e) || !nrow(as.data.frame(e))){ cat(sprintf("  %s: no KEGG result\n", label)); return(NULL) }
  r <- as.data.frame(e); r$suffix <- sub("^hsa","",r$ID)
  cat(sprintf("  %s: %d pathways tested, %d at q<0.05\n", label, nrow(r), sum(r$qvalue<0.05, na.rm=TRUE))); r }
k_asd <- kegg_of(asd_ent, "ASD CHD8 targets")
k_scz <- kegg_of(scz_ent, "SCZ DMP genes")

mb <- function(a){ d <- read.delim(file.path(AX, sprintf("KO_pathway_DA_%s_filtered.tsv", a)), check.names=FALSE)
  d$suffix <- sub("^ko","", d$feature); d }
res <- list()
run <- function(arm, kepi){
  cat(sprintf("\n=========== %s ===========\n", arm))
  if(is.null(kepi)){ cat("  no epigenomic KEGG result; test not available\n"); return(NULL) }
  m <- mb(arm); pool <- intersect(kepi$suffix, m$suffix)
  epi_sig <- intersect(kepi$suffix[kepi$qvalue < 0.05], pool)
  mic_sig <- intersect(m$suffix[m$qval_refiltered < 0.05], pool)
  ov <- intersect(epi_sig, mic_sig)
  cat(sprintf("  testable in both : %d\n  epigenomic sig   : %d\n  microbial sig    : %d\n  overlap          : %d\n",
              length(pool), length(epi_sig), length(mic_sig), length(ov)))
  if(!length(pool) || !length(epi_sig) || !length(mic_sig)){ cat("  no test possible\n")
    return(data.frame(arm=arm, pool=length(pool), epi=length(epi_sig), mic=length(mic_sig), overlap=length(ov), p=NA, fold=NA)) }
  p <- phyper(length(ov)-1, length(epi_sig), length(pool)-length(epi_sig), length(mic_sig), lower.tail=FALSE)
  fold <- (length(ov)/length(mic_sig)) / (length(epi_sig)/length(pool))
  cat(sprintf("  hypergeometric P = %.4f | %.1f-fold\n", p, fold))
  if(length(ov)){ nm <- kepi$Description[match(ov, kepi$suffix)]
    for(i in seq_along(ov)) cat(sprintf("    ko%s / hsa%s  %s\n", ov[i], ov[i], nm[i])) }
  write.csv(data.frame(suffix=pool, epi_q=kepi$qvalue[match(pool,kepi$suffix)],
                       mic_q=m$qval_refiltered[match(pool,m$suffix)],
                       description=kepi$Description[match(pool,kepi$suffix)]),
            file.path(OUT, sprintf("microbiome_epigenomics_%s.csv", arm)), row.names=FALSE)
  data.frame(arm=arm, pool=length(pool), epi=length(epi_sig), mic=length(mic_sig), overlap=length(ov), p=p, fold=fold) }
res$ASD <- run("ASD", k_asd); res$SCZ <- run("SCZ", k_scz)
cat("\n=========== MDD ===========\n  no microbial pathway survived correction; no test available\n")
S <- do.call(rbind, res)
if(!is.null(S)){ cat("\n=========== SUMMARY ===========\n"); print(S, row.names=FALSE)
  write.csv(S, file.path(OUT,"microbiome_epigenomics_summary.csv"), row.names=FALSE) }
cat(sprintf("\nwritten to %s\n", OUT))
