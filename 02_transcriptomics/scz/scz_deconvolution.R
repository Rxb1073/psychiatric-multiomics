#!/usr/bin/env Rscript
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({ library(dtangle); library(brainyR) })
logf <- function(...) cat(sprintf(...), "\n"); CASE <- "SCZ"
NEURO <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
CS <- file.path(NEURO,"counts_scz"); SH <- file.path(CS,"sheets")
OR <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
OUT <- file.path(CS,"Deconvolution"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
e <- new.env(); data("sigsBrain", package="brainyR", envir=e); sigs <- get("sigsBrain", envir=e)
USE_SIGS <- c("MB","CA")
for (s in USE_SIGS) logf("signature %-3s : %d genes x %d types (%s)", s, nrow(sigs[[s]]), ncol(sigs[[s]]), paste(colnames(sigs[[s]]),collapse=", "))
bt <- read.delim(file.path(OR,"gene_id2biotype.tsv"), header=FALSE, col.names=c("gid","bt")); btv <- setNames(bt$bt,bt$gid)
len <- read.delim(file.path(OR,"gene_lengths.tsv"), header=FALSE, col.names=c("gid","len")); lenv <- setNames(as.numeric(len$len),len$gid)
read_counts <- function(p){ cm<-read.delim(p,comment.char="#",check.names=FALSE); g<-cm[[1]];cm[[1]]<-NULL
  cm<-cm[,!(colnames(cm)%in%c("Chr","Start","End","Strand","Length")),drop=FALSE]
  h<-regmatches(colnames(cm),regexpr("SRR[0-9]+",colnames(cm))); cm<-cm[,nchar(h)>0,drop=FALSE]
  colnames(cm)<-regmatches(colnames(cm),regexpr("SRR[0-9]+",colnames(cm))); rownames(cm)<-g
  m<-as.matrix(cm);m[is.na(m)]<-0;m }
datasets <- list(
  GSE80655_DLPFC=list(counts=file.path(NEURO,"PRJNA319583/counts/PRJNA319583_gene_counts_matrix.tsv"),
    sheet=file.path(SH,"GSE80655_DLPFC_samplesheet.tsv"), covars=c("age","sex")),
  GSE78936=list(counts=file.path(NEURO,"PRJNA314463/counts/PRJNA314463_gene_counts.txt"),
    sheet=file.path(SH,"GSE78936_samplesheet.tsv"), covars=c("brain_region")),
  GSE138082=list(counts=file.path(NEURO,"PRJNA574470/counts/PRJNA574470_gene_counts_matrix.tsv"),
    sheet=file.path(SH,"GSE138082_samplesheet.tsv"), covars=c("age","sex")))
all_props <- list(); all_tests <- list()
for (dsn in names(datasets)) {
  d <- datasets[[dsn]]; logf("\n=========== %s ===========", dsn)
  ss <- read.delim(d$sheet, check.names=FALSE); rownames(ss) <- ss$srr
  cm <- read_counts(d$counts); k <- intersect(colnames(cm),rownames(ss))
  cm <- cm[,k,drop=FALSE]; ss <- ss[k,,drop=FALSE]; logf("samples: %d", ncol(cm))
  g <- rownames(cm)[!is.na(btv[rownames(cm)]) & btv[rownames(cm)]%in%c("protein_coding","lncRNA")]
  cm <- cm[g,,drop=FALSE]
  lnc <- g[btv[g]=="lncRNA"]; pc <- g[btv[g]=="protein_coding"]
  ss$q_lnc <- colSums(cm[lnc,,drop=FALSE])/colSums(cm)
  Lp <- lenv[pc]; ok <- !is.na(Lp)&Lp>0; W <- cm[pc,,drop=FALSE][ok,,drop=FALSE]; Lo <- Lp[ok]
  ss$q_len <- apply(W,2,function(x){s<-sum(x); if(is.na(s)||s==0) NA_real_ else sum(x*Lo)/s})
  cpm <- t(t(cm)/(colSums(cm)/1e6))
  DROP_CT <- list(MB=character(0), CA=c("Neurons"))
  for (sg in USE_SIGS) {
    ref <- as.matrix(sigs[[sg]]); drop <- intersect(DROP_CT[[sg]],colnames(ref))
    if (length(drop)){logf("  [%s] drop nested: %s",sg,paste(drop,collapse=", ")); ref<-ref[,setdiff(colnames(ref),drop),drop=FALSE]}
    common <- intersect(rownames(cpm),rownames(ref)); logf("  [%s] shared genes: %d/%d",sg,length(common),nrow(ref))
    if (length(common)<1000){logf("  [%s] too few - skip",sg); next}
    Y <- t(cpm[common,,drop=FALSE]); Ref <- t(ref[common,,drop=FALSE])
    Ylog <- log2(Y+1); Reflog <- log2(Ref+1)
    pure <- lapply(seq_len(nrow(Reflog)),function(i) i); names(pure) <- rownames(Reflog)
    combined <- rbind(Reflog,Ylog)
    mk <- dtangle::find_markers(combined, pure_samples=pure, data_type="rna-seq", marker_method="ratio")
    nmk <- as.numeric(pmax(3, round(0.1*lengths(mk$L))))
    cat(sprintf("  [%s] markers/type: %s\n",sg,paste(nmk,collapse=", ")))
    res <- dtangle::dtangle(Y=combined, pure_samples=pure, n_markers=nmk, markers=mk$L, data_type="rna-seq")
    P <- res$estimates[(nrow(Reflog)+1):nrow(combined),,drop=FALSE]; rownames(P)<-rownames(Ylog); colnames(P)<-rownames(Reflog)
    if (any(is.na(P))){logf("  [%s] NA props - skip",sg); next}
    logf("  [%s] median props: %s",sg,paste(sprintf("%s=%.3f",colnames(P),apply(P,2,median)),collapse="  "))
    df <- data.frame(srr=rownames(P), dataset=dsn, signature=sg, P, check.names=FALSE)
    keep_meta <- c("diagnosis","sex","age","q_lnc","q_len", intersect(c("brain_region","pmi","rin","ph"),names(ss)))
    df <- merge(df, ss[,c("srr",keep_meta)], by="srr")
    all_props[[paste(dsn,sg,sep="_")]] <- df
    if (length(unique(df$diagnosis))<2){logf("  [%s] single group - no test",sg); next}
    df$diagnosis <- relevel(factor(df$diagnosis), ref="Control")
    if ("age"%in%names(df)) df$age_z <- as.numeric(scale(suppressWarnings(as.numeric(df$age))))
    df$q_lnc_z <- as.numeric(scale(df$q_lnc)); df$q_len_z <- as.numeric(scale(df$q_len))
    extra <- intersect(c("brain_region"),names(df))
    covs <- c(if("age"%in%names(df))"age_z", if("sex"%in%names(df))"sex", extra)
    covs <- covs[sapply(covs,function(v) v%in%names(df) && length(unique(df[[v]]))>1)]
    for (ct in colnames(P)) {
      y <- df[[ct]]; if (all(is.na(y))||sd(y,na.rm=TRUE)==0) next
      f_adj <- as.formula(paste("y ~ diagnosis +", paste(c(covs,"q_len_z","q_lnc_z"),collapse=" + ")))
      m1 <- tryCatch(lm(y~diagnosis,data=df),error=function(e)NULL)
      m2 <- tryCatch(lm(f_adj,data=df),error=function(e)NULL)
      if (is.null(m1)||is.null(m2)){logf("  [%s] %s: lm failed",sg,ct); next}
      c1<-summary(m1)$coefficients; c2<-summary(m2)$coefficients
      r <- grep("^diagnosisSCZ",rownames(c2),value=TRUE)
      if (length(r)!=1) next
      all_tests[[length(all_tests)+1]] <- data.frame(dataset=dsn,signature=sg,cell_type=ct,
        median_SCZ=median(y[df$diagnosis=="SCZ"]), median_Ctrl=median(y[df$diagnosis=="Control"]),
        beta_naive=c1["diagnosisSCZ","Estimate"], p_naive=c1["diagnosisSCZ","Pr(>|t|)"],
        beta_adj=c2[r,"Estimate"], p_adj=c2[r,"Pr(>|t|)"]) }
  }
}
META <- c("srr","dataset","signature","diagnosis","sex","age","brain_region","pmi","rin","ph","q_lnc","q_len")
to_long <- function(df){ for(m in META) if(!m%in%names(df)) df[[m]]<-NA
  ct_cols <- setdiff(names(df),META)
  do.call(rbind, lapply(ct_cols, function(ct) data.frame(df[,META,drop=FALSE], cell_type=ct, proportion=df[[ct]], row.names=NULL, check.names=FALSE))) }
long <- do.call(rbind, lapply(all_props, to_long))
write.csv(long, file.path(OUT,"cell_type_proportions_long.csv"), row.names=FALSE)
logf("proportions written: %d rows", nrow(long))
tests <- do.call(rbind, all_tests)
tests$p_naive_BH <- p.adjust(tests$p_naive,"BH"); tests$p_adj_BH <- p.adjust(tests$p_adj,"BH")
tests <- tests[order(tests$p_adj),]
write.csv(tests, file.path(OUT,"proportion_tests.csv"), row.names=FALSE)
logf("\n=== PROPORTION TESTS (adj = PRIMARY) ===")
print(format(tests[,c("dataset","signature","cell_type","median_SCZ","median_Ctrl","p_naive","p_adj","p_adj_BH")],digits=3), row.names=FALSE)
logf("\n * significant in BOTH signatures + adjusted model = candidate real composition diff")
logf(" * naive-sig but not adjusted = explained by quality/covariates")
pdf(file.path(OUT,"proportions_by_diagnosis.pdf"), width=11, height=6)
for (n in names(all_props)) { df <- all_props[[n]]; if (length(unique(df$diagnosis))<2) next
  cts <- setdiff(colnames(df), META)
  par(mfrow=c(2,ceiling(length(cts)/2)), mar=c(4,4,3,1))
  for (ct in cts) boxplot(df[[ct]]~df$diagnosis, main=paste(n,ct), ylab="proportion", xlab="", col=c("grey85","indianred"))
  par(mfrow=c(1,1)) }
dev.off()
sink(file.path(OUT,"sessionInfo.txt")); print(sessionInfo()); sink()
logf("\noutputs -> %s", OUT)
