.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)})
logf <- function(...) cat(sprintf(...),"\n")
CS <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/counts_scz"
OUT <- file.path(CS,"META"); SYN <- "GO:0097060"

## --- 1. per-cohort GSEA: is the synaptic signature in all three, or one? ---
files <- list(CA1  =file.path(CS,"HIPP_donor","dec_ca1.csv"),
              DLPFC=file.path(CS,"PERCOHORT_prua","GSE80655_D2_full.csv"),
              BA11 =file.path(CS,"PERCOHORT_prua","GSE78936_BA11_F1_full.csv"))
logf("=========== PER-COHORT GSEA (GO-CC, synaptic membrane) ===========")
for (nm in names(files)) {
  d <- read.csv(files[[nm]], stringsAsFactors=FALSE)
  gcol <- if ("gene" %in% names(d)) "gene" else names(d)[1]
  d$g <- sub("\\..*","", d[[gcol]])
  d <- d[!is.na(d$pvalue) & !is.na(d$log2FoldChange),]
  d$e <- AnnotationDbi::mapIds(org.Hs.eg.db, d$g, "ENTREZID","ENSEMBL", multiVals="first")
  d <- d[!is.na(d$e) & !duplicated(d$e),]
  gl <- sort(setNames(sign(d$log2FoldChange)*-log10(pmax(d$pvalue,1e-300)), d$e), decreasing=TRUE)
  g <- tryCatch(gseGO(gl, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont="CC",
        pAdjustMethod="BH", pvalueCutoff=1, verbose=FALSE), error=function(e) NULL)
  if (is.null(g)) { logf("  %-6s : failed", nm); next }
  r <- g@result; hit <- r[r$ID==SYN,]
  logf("  %-6s : %4d terms at FDR<0.05 | synaptic membrane NES=%s p.adj=%s",
       nm, sum(r$p.adjust<0.05),
       if(nrow(hit)) sprintf("%+.2f",hit$NES) else "NA",
       if(nrow(hit)) format.pval(hit$p.adjust,3) else "NA")
}

## --- 2. does the meta signal survive neuronal-proportion adjustment? ---
logf("\n=========== NEURONAL PROPORTION BY DIAGNOSIS (unique donors) ===========")
RT <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/GSE80655/metadata/PRJNA319583_SraRunTable.csv"
pr <- read.csv(file.path(CS,"Deconvolution","cell_type_proportions_long.csv"), stringsAsFactors=FALSE)
pr <- pr[, setdiff(names(pr), c("diagnosis","study"))]
neu <- unique(pr$cell_type)[grepl("Neuro|Excit|Inhib", unique(pr$cell_type))]
logf("  neuronal cell types found: %s", paste(neu, collapse=", "))
sh <- do.call(rbind, lapply(c("GSE80655_DLPFC","GSE78936"), function(s){
  x <- read.delim(file.path(CS,"sheets",paste0(s,"_samplesheet.tsv")), check.names=FALSE)
  data.frame(srr=x$srr, diagnosis=x$diagnosis, region=x$brain_region, study=sub("_DLPFC","",s))}))
m <- merge(subset(pr, cell_type %in% neu), sh, by="srr")
rt <- read.csv(RT, check.names=FALSE); mp <- setNames(rt$BioSample, rt$Run)
m$bio <- ifelse(is.na(mp[m$srr]), m$srr, mp[m$srr])
m <- aggregate(proportion ~ bio+signature+cell_type+diagnosis+region+study, m, mean)
u <- subset(m, study=="GSE80655" | region=="Brodmann's Area 11")
for (s in unique(u$signature)) for (ct in neu) {
  x <- subset(u, signature==s & cell_type==ct); if (nrow(x)<10) next
  logf("  %-3s %-12s n=%d | Ctrl %.4f vs SCZ %.4f | Wilcoxon p = %.4f", s, ct, nrow(x),
       median(x$proportion[x$diagnosis=="Control"]), median(x$proportion[x$diagnosis=="SCZ"]),
       wilcox.test(proportion~diagnosis, data=x)$p.value)
}
