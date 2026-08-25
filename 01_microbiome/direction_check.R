setwd("/rds/projects/e/elhamsak-dissertation-autism/Microbiome/pooled")
otu <- read.delim("merged_export/feature-table.tsv", skip=1, row.names=1, check.names=FALSE)
otu <- otu[!grepl("^Unassigned", rownames(otu)), , drop=FALSE]
otu <- as.matrix(otu)                                  # matrix => single-row indexing gives a vector
meta <- read.delim("metadata_pooled.tsv", row.names=1, check.names=FALSE)
meta <- meta[rownames(meta) != "#q2:types", , drop=FALSE]
meta <- meta[colnames(otu), , drop=FALSE]

core <- read.csv("stage5_concordant_core.csv")
studies <- c("PRJNA453621","PRJNA529598","PRJNA769228","PRJNA1022510")
rel <- sweep(otu, 2, colSums(otu), "/")                # relative abundance

cat(sprintf("%-34s %-7s | %s\n","genus","pooled",paste(sprintf("%-9s",sub("PRJNA","",studies)),collapse=" ")))
cat(strrep("-",82),"\n")
for(i in seq_len(nrow(core))){
  tax <- core$taxon.x[i]
  if(!tax %in% rownames(rel)) next
  gname <- sub(".*;g__","g__",tax)
  pooled_sign <- ifelse(core$ancom_lfc[i]>0,"+","-")
  signs <- c()
  for(s in studies){
    samp <- rownames(meta)[meta$study==s]
    asd_s <- samp[meta[samp,"diagnosis"]=="ASD"]
    ctl_s <- samp[meta[samp,"diagnosis"]=="Control"]
    asd <- as.numeric(rel[tax, asd_s])
    ctl <- as.numeric(rel[tax, ctl_s])
    if((sum(asd,na.rm=T)==0 && sum(ctl,na.rm=T)==0)){
      signs <- c(signs,"·")
    } else {
      d <- median(asd,na.rm=T)-median(ctl,na.rm=T)
      signs <- c(signs, ifelse(d>0,"+",ifelse(d<0,"-","0")))
    }
  }
  agree <- sum(signs==pooled_sign); present <- sum(signs %in% c("+","-"))
  cat(sprintf("%-34s %-7s | %s  %d/%d\n",
      substr(gname,1,34), paste0("ASD",pooled_sign),
      paste(sprintf("%-9s",signs),collapse=" "), agree, present))
}
cat("\n(+ = higher in ASD in that cohort, - = lower, · = ~absent) | last col = cohorts agreeing with pooled direction\n")
