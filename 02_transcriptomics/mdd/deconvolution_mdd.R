#!/usr/bin/env Rscript
.libPaths(c("/rds/projects/e/elhamsak-dissertation-autism-2/Rlibs", .libPaths()))
suppressMessages({ library(dtangle); library(brainyR) })
logf <- function(...) cat(sprintf(...), "\n")

A2   <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics"
DE   <- file.path(A2, "DE")
ANN  <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
OUT  <- file.path(A2, "Deconvolution")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

e <- new.env(); data("sigsBrain", package="brainyR", envir=e)
sigs <- get("sigsBrain", envir=e)
USE_SIGS <- c("MB","CA")
for (s in USE_SIGS)
  logf("signature %-3s : %d genes x %d cell types (%s)", s,
       nrow(sigs[[s]]), ncol(sigs[[s]]), paste(colnames(sigs[[s]]), collapse=", "))

bt   <- read.delim(file.path(ANN,"gene_id2biotype.tsv"), header=FALSE, col.names=c("gid","bt"))
btv  <- setNames(bt$bt, bt$gid)
len  <- read.delim(file.path(ANN,"gene_lengths.tsv"), header=FALSE, col.names=c("gid","len"))
lenv <- setNames(as.numeric(len$len), len$gid)

col_all <- read.delim(file.path(DE,"rna_colData.tsv"), check.names=FALSE)
rownames(col_all) <- col_all$`sample-id`
cnt_all <- read.delim(file.path(DE,"merged_counts.tsv"), row.names=1, check.names=FALSE)
datasets <- c("GSE80655","GSE101521","GSE102556")
all_props <- list(); all_tests <- list()

for (dsn in datasets) {
  logf("\n=================== %s ===================", dsn)
  ss <- col_all[col_all$study==dsn, , drop=FALSE]
  cm <- as.matrix(cnt_all[, rownames(ss), drop=FALSE]); mode(cm) <- "numeric"
  cm[is.na(cm)] <- 0
  logf("samples: %d", ncol(cm))
  g  <- rownames(cm)[!is.na(btv[rownames(cm)]) & btv[rownames(cm)] %in% c("protein_coding","lncRNA")]
  cm <- cm[g, , drop=FALSE]
  lnc <- g[btv[g]=="lncRNA"]; pc <- g[btv[g]=="protein_coding"]
  ss$q_lnc <- colSums(cm[lnc, , drop=FALSE]) / colSums(cm)
  Lp <- lenv[pc]; ok <- !is.na(Lp) & Lp>0
  W  <- cm[pc, , drop=FALSE][ok, , drop=FALSE]; Lo <- Lp[ok]
  ss$q_len <- apply(W, 2, function(x){ s<-sum(x); if(is.na(s)||s==0) NA_real_ else sum(x*Lo)/s })
  cpm <- t(t(cm)/(colSums(cm)/1e6))
  DROP_CT <- list(MB=character(0), CA=c("Neurons"))

  for (sg in USE_SIGS) {
    ref <- as.matrix(sigs[[sg]])
    drop <- intersect(DROP_CT[[sg]], colnames(ref))
    if (length(drop)) { logf("  [%s] dropping nested: %s", sg, paste(drop,collapse=", "))
                        ref <- ref[, setdiff(colnames(ref),drop), drop=FALSE] }
    common <- intersect(rownames(cpm), rownames(ref))
    logf("  [%s] genes shared with signature: %d / %d", sg, length(common), nrow(ref))
    if (length(common) < 1000) { logf("  [%s] TOO FEW shared genes - skip", sg); next }
    Y   <- t(cpm[common, , drop=FALSE]); Ref <- t(ref[common, , drop=FALSE])
    Ylog <- log2(Y+1); Reflog <- log2(Ref+1)
    pure <- lapply(seq_len(nrow(Reflog)), function(i) i); names(pure) <- rownames(Reflog)
    combined <- rbind(Reflog, Ylog)
    mk  <- dtangle::find_markers(combined, pure_samples=pure, data_type="rna-seq", marker_method="ratio")
    nmk <- as.numeric(pmax(3, round(0.1*lengths(mk$L))))
    logf("  [%s] markers per cell type: %s", sg, paste(nmk, collapse=", "))
    res <- dtangle::dtangle(Y=combined, pure_samples=pure, n_markers=nmk, markers=mk$L, data_type="rna-seq")
    P <- res$estimates[(nrow(Reflog)+1):nrow(combined), , drop=FALSE]
    rownames(P) <- rownames(Ylog); colnames(P) <- rownames(Reflog)
    if (any(is.na(P))) { logf("  [%s] WARNING NA proportions - skip", sg); next }
    logf("  [%s] median proportions: %s", sg,
         paste(sprintf("%s=%.3f", colnames(P), apply(P,2,median)), collapse="  "))
    df <- data.frame(`sample-id`=rownames(P), dataset=dsn, signature=sg, P, check.names=FALSE)
    keep_meta <- intersect(c("disease","sex","age","pmi","rin","ph","q_lnc","q_len"), names(ss))
    ss2 <- ss; ss2$`sample-id` <- rownames(ss2)
    df <- merge(df, ss2[, c("sample-id", keep_meta)], by="sample-id")
    all_props[[paste(dsn,sg,sep="_")]] <- df
    if (length(unique(df$disease)) < 2) { logf("  [%s] single group - no test", sg); next }
    df$disease <- relevel(factor(df$disease), ref="Control")
    numify <- function(v){ n<-suppressWarnings(as.numeric(v)); if(all(!is.na(n))) n else v }
    for (cv in c("age","q_lnc","q_len")) if (cv %in% names(df)) df[[paste0(cv,"_z")]] <- as.numeric(scale(numify(df[[cv]])))
    covs <- c("age_z","sex"); covs <- covs[sapply(covs, function(v) v %in% names(df) && length(unique(df[[v]]))>1)]
    qterms <- intersect(c("q_len_z","q_lnc_z"), names(df))
    for (ct in colnames(P)) {
      y <- df[[ct]]
      if (all(is.na(y)) || sd(y,na.rm=TRUE)==0) next
      f_naive <- y ~ disease
      f_adj   <- as.formula(paste("y ~ disease +", paste(c(covs,qterms), collapse=" + ")))
      m1 <- tryCatch(lm(f_naive, data=df), error=function(e) NULL)
      m2 <- tryCatch(lm(f_adj,   data=df), error=function(e) NULL)
      if (is.null(m1)||is.null(m2)) next
      c1 <- summary(m1)$coefficients; c2 <- summary(m2)$coefficients
      r  <- grep("^diseaseMDD", rownames(c2), value=TRUE)
      if (!length(r) || !"diseaseMDD" %in% rownames(c1)) next
      all_tests[[length(all_tests)+1]] <- data.frame(
        dataset=dsn, signature=sg, cell_type=ct,
        median_MDD =median(y[df$disease=="MDD"]),
        median_Ctrl=median(y[df$disease=="Control"]),
        beta_naive=c1["diseaseMDD","Estimate"], p_naive=c1["diseaseMDD","Pr(>|t|)"],
        beta_adj  =c2[r,"Estimate"],           p_adj  =c2[r,"Pr(>|t|)"])
    }
  }
}

META <- c("sample-id","dataset","signature","disease","sex","age","pmi","rin","ph","q_lnc","q_len")
to_long <- function(df){
  for (m in META) if (!m %in% names(df)) df[[m]] <- NA
  ct_cols <- setdiff(names(df), META)
  do.call(rbind, lapply(ct_cols, function(ct)
    data.frame(df[,META,drop=FALSE], cell_type=ct, proportion=df[[ct]], row.names=NULL, check.names=FALSE)))
}
long <- do.call(rbind, lapply(all_props, to_long))
write.csv(long, file.path(OUT,"cell_type_proportions_long.csv"), row.names=FALSE)
for (sg in unique(long$signature)) {
  sub <- long[long$signature==sg, ]
  w <- reshape(sub[,c("sample-id","dataset","disease","sex","age","cell_type","proportion")],
               idvar=c("sample-id","dataset","disease","sex","age"), timevar="cell_type", direction="wide")
  names(w) <- sub("^proportion\\.","", names(w))
  write.csv(w, file.path(OUT, sprintf("cell_type_proportions_%s.csv", sg)), row.names=FALSE)
}
logf("proportions written: %d rows (long)", nrow(long))

tests <- do.call(rbind, all_tests)
tests$p_naive_BH <- p.adjust(tests$p_naive, "BH")
tests$p_adj_BH   <- p.adjust(tests$p_adj,   "BH")
tests$p_naive_BH_sig <- ave(tests$p_naive, tests$signature, FUN=function(p) p.adjust(p,"BH"))
tests$p_adj_BH_sig   <- ave(tests$p_adj,   tests$signature, FUN=function(p) p.adjust(p,"BH"))
tests <- tests[order(tests$p_adj), ]
write.csv(tests, file.path(OUT,"proportion_tests.csv"), row.names=FALSE)
logf("\n================ PROPORTION TESTS (MDD vs Control) ================")
logf("naive = proportion ~ disease   |   adj = + age + sex + q_len + q_lnc  (PRIMARY)")
print(format(tests[, c("dataset","signature","cell_type","median_MDD","median_Ctrl","p_naive","p_adj","p_adj_BH_sig")], digits=3), row.names=FALSE)
logf("\n================ NEURONAL PROPORTION (key for IEG interpretation) ================")
neuro <- tests[tests$cell_type %in% c("Neurons","Excitatory","Inhibitory"), ]
print(format(neuro[, c("dataset","signature","cell_type","median_MDD","median_Ctrl","p_naive","p_adj")], digits=3), row.names=FALSE)
logf("\nINTERPRETATION:")
logf("  * Neuronal proportion NOT lower in MDD  -> IEG downregulation is per-neuron (regulatory).")
logf("  * Neuronal proportion LOWER in MDD      -> part of IEG signal is compositional.")
logf("\n================ SUMMARY ================")
for (col in c("p_adj_BH","p_adj_BH_sig")) {
  s <- tests[tests[[col]] < 0.05, ]
  logf("  %-16s < 0.05: %d test(s)%s", col, nrow(s),
       if(nrow(s)) paste0("  [", paste(unique(paste(s$dataset,s$signature,s$cell_type)),collapse="; "), "]") else "")
}
pdf(file.path(OUT,"proportions_by_diagnosis.pdf"), width=11, height=6)
for (n in names(all_props)) {
  df <- all_props[[n]]; if (length(unique(df$disease))<2) next
  cts <- setdiff(colnames(df), META)
  par(mfrow=c(2, ceiling(length(cts)/2)), mar=c(4,4,3,1))
  for (ct in cts) boxplot(df[[ct]] ~ df$disease, main=paste(n,ct), ylab="proportion", xlab="", col=c("grey85","indianred"))
  par(mfrow=c(1,1))
}
dev.off()
sink(file.path(OUT,"sessionInfo.txt")); print(sessionInfo()); sink()
logf("\noutputs -> %s", OUT); logf("=== done ===")
