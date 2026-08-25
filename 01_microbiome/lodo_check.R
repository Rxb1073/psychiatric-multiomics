suppressMessages({library(ANCOMBC); library(TreeSummarizedExperiment); library(Maaslin2)})
setwd("/rds/projects/e/elhamsak-dissertation-autism/Microbiome/pooled")

otu <- read.delim("merged_export/feature-table.tsv", skip=1, row.names=1, check.names=FALSE)
otu <- otu[!grepl("^Unassigned", rownames(otu)), , drop=FALSE]
otu <- round(as.matrix(otu))
meta <- read.delim("metadata_pooled.tsv", row.names=1, check.names=FALSE)
meta <- meta[rownames(meta) != "#q2:types", , drop=FALSE]
meta <- meta[colnames(otu), , drop=FALSE]
meta$diagnosis <- factor(meta$diagnosis, levels=c("Control","ASD"))
meta$study <- factor(meta$study)

core <- read.csv("stage5_concordant_core.csv")
core_taxa <- core$taxon.x
core_key  <- make.names(core_taxa)
studies <- levels(meta$study)

run_concordance <- function(otu_s, meta_s, tag){
  meta_s$study <- droplevels(factor(meta_s$study))
  formula <- if(nlevels(meta_s$study) >= 2) "study + diagnosis" else "diagnosis"
  tse <- TreeSummarizedExperiment(assays=list(counts=otu_s), colData=meta_s)
  ab <- tryCatch(ancombc2(data=tse, assay_name="counts", fix_formula=formula,
                 group="diagnosis", struc_zero=TRUE, prv_cut=0.10,
                 p_adj_method="BH", verbose=FALSE),
                 error=function(e){cat("ANCOM fail",tag,conditionMessage(e),"\n");NULL})
  if(is.null(ab)) return(NULL)
  ar <- ab$res
  lc <- grep("^lfc_diagnosis",names(ar),value=TRUE)[1]
  qc <- grep("^diff_diagnosis",names(ar),value=TRUE)[1]
  a_sig <- ar[ar[[qc]]==TRUE & !is.na(ar[[qc]]), ]
  a <- data.frame(key=make.names(a_sig$taxon), alfc=a_sig[[lc]])
  refstudy <- levels(meta_s$study)[1]
  m <- tryCatch({
    fit <- Maaslin2(input_data=as.data.frame(t(otu_s)), input_metadata=meta_s,
       output=paste0("lodo_maaslin_",tag),
       fixed_effects=if(nlevels(meta_s$study)>=2) c("study","diagnosis") else c("diagnosis"),
       reference=paste0("study,",refstudy,";diagnosis,Control"),
       normalization="TSS", transform="LOG", min_prevalence=0.10,
       plot_heatmap=FALSE, plot_scatter=FALSE)
    fit$results
  }, error=function(e){cat("MaAsLin fail",tag,conditionMessage(e),"\n");NULL})
  if(is.null(m)) return(NULL)
  m_dx <- m[m$metadata=="diagnosis" & m$qval<0.05, ]
  mm <- data.frame(key=make.names(m_dx$feature), mcoef=m_dx$coef)
  both <- merge(a, mm, by="key")
  both <- both[sign(both$alfc)==sign(both$mcoef), ]
  both$key
}

survmat <- matrix(0, nrow=length(core_key), ncol=length(studies),
                  dimnames=list(core_key, studies))
for(s in studies){
  keep <- rownames(meta)[meta$study != s]
  otu_s <- otu[, keep]; meta_s <- meta[keep, ]
  otu_s <- otu_s[rowSums(otu_s)>0, ]
  cat("=== dropping", s, "| n =", ncol(otu_s), "===\n")
  conc <- run_concordance(otu_s, meta_s, gsub("PRJNA","",s))
  if(!is.null(conc)) survmat[core_key %in% conc, s] <- 1
}

surv <- rowSums(survmat)
out <- data.frame(genus=sub(".*;g__","g__",core_taxa), taxon=core_taxa,
                  survives_of_4=surv,
                  broken_by=apply(survmat,1,function(r) paste(gsub("PRJNA","",studies[r==0]),collapse=",")))
out <- out[order(-out$survives_of_4), ]
cat("\n===== LEAVE-ONE-DATASET-OUT survival =====\n")
for(i in seq_len(nrow(out)))
  cat(sprintf("  %-34s %d/4  %s\n", substr(out$genus[i],1,34), out$survives_of_4[i],
              ifelse(out$survives_of_4[i]<4, paste("(lost when dropping:",out$broken_by[i],")"), "")))
write.csv(out, "stage5_lodo_survival.csv", row.names=FALSE)
cat("\nSaved -> stage5_lodo_survival.csv\n")
