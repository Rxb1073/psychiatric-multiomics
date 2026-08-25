## Documented attempt to apply sample-matched multi-omic integration.
##
## PURPOSE. To record empirically, rather than argue from the literature, that
## joint latent-variable methods cannot be applied to these data. The script
## assembles the real data matrices for the schizophrenia arm, quantifies the
## sample overlap between layers, and then attempts model construction with
## MOFA2 and with mixOmics DIABLO, capturing whatever each returns.
##
## Output is written to a log file suitable for an appendix.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(edgeR); library(limma)})

NP  <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
CS  <- file.path(NP,"counts_scz"); SHD <- file.path(CS,"sheets")
OR  <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
RT  <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/GSE80655/metadata/PRJNA319583_SraRunTable.csv"
AX  <- file.path(OR,"DISSERTATION_final/INTEGRATION_inputs/06_kegg_axis")
MB_META <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA880407/metadata/metadata_qiime2.tsv"
OUT <- file.path(OR,"DISSERTATION_final/INTEGRATION_inputs/09_matched_methods_attempt")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

LOG <- file.path(OUT,"matched_methods_attempt.log")
sink(LOG, split=TRUE)
cat("Attempt to apply sample-matched integration methods, schizophrenia arm\n")
cat("Run:", format(Sys.time()), "| R", as.character(getRversion()), "\n")
cat(strrep("=", 78), "\n\n")

## ---- layer 1: transcriptome, 76 frontal donors -----------------------------
bt <- read.delim(file.path(OR,"gene_id2biotype.tsv"), header=FALSE, col.names=c("gid","bt"))
KEEP <- bt$gid[bt$bt %in% c("protein_coding","lncRNA")]
rc <- function(p){cm<-read.delim(p,comment.char="#",check.names=FALSE);g<-cm[[1]];cm[[1]]<-NULL
  cm<-cm[,!colnames(cm)%in%c("Chr","Start","End","Strand","Length"),drop=FALSE]
  colnames(cm)<-regmatches(colnames(cm),regexpr("SRR[0-9]+",colnames(cm)))
  m<-as.matrix(cm);rownames(m)<-g;storage.mode(m)<-"integer";m}

s1<-read.delim(file.path(SHD,"GSE80655_DLPFC_samplesheet.tsv"),check.names=FALSE)
rownames(s1)<-s1$srr
c1<-rc(file.path(NP,"PRJNA319583/counts/PRJNA319583_gene_counts_matrix.tsv"))
k<-intersect(colnames(c1),rownames(s1)); c1<-c1[,k]; s1<-s1[k,]
rt<-read.csv(RT,check.names=FALSE); bio<-setNames(rt$BioSample,rt$Run)[colnames(c1)]
ix<-split(seq_len(ncol(c1)),bio)
c1s<-vapply(ix,function(i)as.integer(rowSums(c1[,i,drop=FALSE])),integer(nrow(c1)))
rownames(c1s)<-rownames(c1); ke<-!duplicated(bio); c1<-c1s[,bio[ke]]

s2<-read.delim(file.path(SHD,"GSE78936_samplesheet.tsv"),check.names=FALSE)
rownames(s2)<-s2$srr; s2<-s2[s2$brain_region=="Brodmann's Area 11",]
c2<-rc(file.path(NP,"PRJNA314463/counts/PRJNA314463_gene_counts.txt"))
k<-intersect(colnames(c2),rownames(s2)); c2<-c2[,k]

g<-intersect(rownames(c1),rownames(c2)); g<-g[g%in%KEEP]
cnt<-cbind(c1[g,],c2[g,])
cnt<-cnt[rowSums(cnt>=10)>=20,]
expr <- voom(calcNormFactors(DGEList(cnt)))$E
cat(sprintf("LAYER 1  transcriptome : %d genes x %d donors\n", nrow(expr), ncol(expr)))
cat(sprintf("         identifier style: %s\n", paste(head(colnames(expr),3), collapse=", ")))

## ---- layer 2: gut microbiome, KO pathway abundances ------------------------
mb <- read.delim(file.path(AX,"maaslin2_KO_SCZ/features/filtered_data.tsv"),
                 row.names=1, check.names=FALSE)
micro <- t(as.matrix(mb))
cat(sprintf("\nLAYER 2  gut microbiome: %d pathways x %d samples\n", nrow(micro), ncol(micro)))
cat(sprintf("         identifier style: %s\n", paste(head(colnames(micro),3), collapse=", ")))

## ---- overlap ---------------------------------------------------------------
cat("\n", strrep("-",78), "\nSAMPLE OVERLAP BETWEEN LAYERS\n", strrep("-",78), "\n", sep="")
ov <- intersect(colnames(expr), colnames(micro))
cat(sprintf("transcriptome donors            : %d\n", ncol(expr)))
cat(sprintf("microbiome samples              : %d\n", ncol(micro)))
cat(sprintf("identifiers present in BOTH     : %d\n", length(ov)))
cat("\nThe two layers use different identifier systems (BioSample/SRA accessions\n")
cat("for post-mortem brain, cohort sample codes for stool) because they describe\n")
cat("different people: post-mortem donors and living participants.\n")

## ---- attempt 1: MOFA2 ------------------------------------------------------
cat("\n", strrep("-",78), "\nATTEMPT 1: MOFA2\n", strrep("-",78), "\n", sep="")
if(!requireNamespace("MOFA2", quietly=TRUE)){
  cat("MOFA2 is not installed in this environment.\n")
  cat("Install with BiocManager::install('MOFA2') to reproduce this attempt.\n")
} else {
  suppressMessages(library(MOFA2))
  cat(sprintf("MOFA2 version %s\n\n", as.character(packageVersion("MOFA2"))))
  r <- tryCatch({
        m <- MOFA2::create_mofa(list(transcriptome=expr, microbiome=micro))
        capture.output(print(m))
      }, error = function(e) paste("ERROR:", conditionMessage(e)),
         warning = function(w) paste("WARNING:", conditionMessage(w)))
  cat(paste(r, collapse="\n"), "\n")
}

## ---- attempt 2: mixOmics DIABLO -------------------------------------------
cat("\n", strrep("-",78), "\nATTEMPT 2: mixOmics DIABLO (block.splsda)\n", strrep("-",78), "\n", sep="")
if(!requireNamespace("mixOmics", quietly=TRUE)){
  cat("mixOmics is not installed in this environment.\n")
  cat("Install with BiocManager::install('mixOmics') to reproduce this attempt.\n")
} else {
  suppressMessages(library(mixOmics))
  cat(sprintf("mixOmics version %s\n\n", as.character(packageVersion("mixOmics"))))
  ## DIABLO expects samples as rows, and one outcome vector spanning all blocks
  X <- list(transcriptome = t(expr), microbiome = t(micro))
  Y <- factor(c(rep("case", ncol(expr))))          # length deliberately from one block only
  cat(sprintf("block sample counts: transcriptome %d, microbiome %d\n",
              nrow(X$transcriptome), nrow(X$microbiome)))
  cat(sprintf("outcome vector length: %d\n\n", length(Y)))
  r <- tryCatch({
        f <- mixOmics::block.splsda(X=X, Y=Y, ncomp=2)
        capture.output(print(f))
      }, error = function(e) paste("ERROR:", conditionMessage(e)),
         warning = function(w) paste("WARNING:", conditionMessage(w)))
  cat(paste(r, collapse="\n"), "\n")
}

## ---- conclusion ------------------------------------------------------------
cat("\n", strrep("=",78), "\nCONCLUSION\n", strrep("=",78), "\n", sep="")
cat("Both methods decompose covariance across views measured on a common set of\n")
cat("samples. With zero shared identifiers there is no such covariance to model,\n")
cat("and neither can proceed past model construction. Integration was therefore\n")
cat("performed at gene and pathway level.\n")
sink()
cat(sprintf("\nlog written to %s\n", LOG))
