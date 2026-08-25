## Is the 146-gene hippocampal result driven by DG, or by method choice?
## Same DESeq2 model (R1) throughout; only the sample subset changes.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(DESeq2)})
logf <- function(...) cat(sprintf(...),"\n")
NP  <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
CS  <- file.path(NP,"counts_scz"); SHD <- file.path(CS,"sheets")
OR  <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
OUT <- file.path(CS,"HIPP_donor"); dir.create(OUT, showWarnings=FALSE)

## ---- series-matrix metadata ----
ln <- readLines(gzfile(file.path(SHD,"GSE138082_series_matrix.txt.gz")))
rows <- list()
for (l in ln[grepl("^!Sample_", ln)]) {
  p <- gsub('"','', strsplit(l,"\t")[[1]]); k <- sub("^!","",p[1])
  if (k=="Sample_title") rows$title <- p[-1]
  else if (k=="Sample_characteristics_ch1") {
    tag <- trimws(sub(":.*","",p[2]))
    rows[[tag]] <- trimws(sub("^[^:]*:","",p[-1]))
  }
}
md <- as.data.frame(rows, stringsAsFactors=FALSE)
md$gsm <- gsub('"','', strsplit(grep("^!Sample_geo_accession", ln, value=TRUE)[1],"\t")[[1]][-1])
md$group <- ifelse(grepl("_SZ_", md$title), "SCZ", "Control")
md$donor <- paste(md$gender, md$age, md$pmi, md$group, sep="_")

rt <- read.csv(file.path(NP,"metadata/PRJNA574470_SraRunTable.csv"), check.names=FALSE)
gcol <- grep("GEO_Accession|Sample.Name|Sample Name", names(rt), value=TRUE)[1]
ss <- read.delim(file.path(SHD,"GSE138082_samplesheet.tsv"), check.names=FALSE)
ss$gsm <- setNames(rt[[gcol]], rt$Run)[ss$srr]
ss <- merge(ss, md[,c("gsm","title","hippocampal_subfield","drug","donor","group")], by="gsm")
rownames(ss) <- ss$srr
logf("metadata joined: %d samples", nrow(ss)); print(table(ss$group, ss$hippocampal_subfield))

## ---- counts ----
bt <- read.delim(file.path(OR,"gene_id2biotype.tsv"), header=FALSE, col.names=c("gid","bt"))
cm <- read.delim(file.path(NP,"PRJNA574470/counts/PRJNA574470_gene_counts_matrix.tsv"),
                 comment.char="#", check.names=FALSE)
g <- cm[[1]]; cm[[1]] <- NULL
cm <- cm[, !colnames(cm) %in% c("Chr","Start","End","Strand","Length"), drop=FALSE]
colnames(cm) <- regmatches(colnames(cm), regexpr("SRR[0-9]+", colnames(cm)))
cnt <- as.matrix(cm); rownames(cnt) <- g; storage.mode(cnt) <- "integer"
k <- intersect(colnames(cnt), rownames(ss)); cnt <- cnt[,k]; ss <- ss[k,]
cnt <- cnt[rownames(cnt) %in% bt$gid[bt$bt %in% c("protein_coding","lncRNA")],]

run <- function(sel, lab, tag) {
  s <- ss[sel,,drop=FALSE]; c <- cnt[,rownames(s),drop=FALSE]
  s$group <- relevel(factor(s$group),"Control")
  for (v in c("age","pmi","rin")) s[[paste0(v,"_z")]] <- as.numeric(scale(as.numeric(s[[v]])))
  c <- c[rowSums(c>=10) >= min(table(s$group)),,drop=FALSE]
  dds <- suppressMessages(DESeq(DESeqDataSetFromMatrix(c, s, ~ age_z+sex+pmi_z+rin_z+group), quiet=TRUE))
  r <- as.data.frame(results(dds, name=tail(resultsNames(dds),1)))
  write.csv(r, file.path(OUT,paste0(tag,".csv")))
  f <- function(l) sum(!is.na(r$padj) & r$padj<0.05 & abs(r$log2FoldChange)>=l)
  logf("  %-34s n=%2d (%s) donors=%2d | padj05:%5d | >=0.3:%5d | >=1.0:%4d",
       lab, nrow(s), paste(table(s$group), collapse="/"), length(unique(s$donor)), f(0), f(0.3), f(1))
  rownames(r)[!is.na(r$padj) & r$padj<0.05]
}

logf("\n=========== SAME MODEL (R1), DIFFERENT SUBSETS ===========")
all78 <- run(rep(TRUE,nrow(ss)),                              "ALL 78 (original R1)",   "dec_all78")
ca1c3 <- run(ss$hippocampal_subfield %in% c("CA1","CA3"),     "CA1+CA3 (52)",           "dec_ca1ca3")
ca1   <- run(ss$hippocampal_subfield=="CA1",                  "CA1 only (26)",          "dec_ca1")
ca3   <- run(ss$hippocampal_subfield=="CA3",                  "CA3 only (26)",          "dec_ca3")
dg    <- run(ss$hippocampal_subfield=="DG",                   "DG only (26)",           "dec_dg")

logf("\n=========== OVERLAP WITH THE ORIGINAL 78-SAMPLE LIST ===========")
for (nm in c("CA1+CA3","CA1","CA3","DG")) {
  x <- switch(nm, "CA1+CA3"=ca1c3, "CA1"=ca1, "CA3"=ca3, "DG"=dg)
  logf("  %-8s : %4d genes | %4d shared with ALL78 (%.0f%% of ALL78)",
       nm, length(x), length(intersect(x,all78)),
       100*length(intersect(x,all78))/max(1,length(all78)))
}
logf("\nOutputs -> %s", OUT)
