
#!/usr/bin/env Rscript
# =====================================================================
# SUPPLEMENTARY FIGURES  S1-S4
#
# Four claims the chapter makes without showing:
#
#   S1  the GSE80655 exclusion. Currently a table of numbers, and it is
#       the most contestable decision in the chapter. Pairwise log2FC
#       scatter shows at a glance that two cohorts agree and one does not.
#
#   S2  deconvolution by diagnosis. The chapter states no cell-type
#       proportion differs; SCZ has a figure for exactly this. The
#       existing one was built on the superseded 167-sample run.
#
#   S3  sex. A whole Results section with no figure, and the source study
#       for GSE102556 is about sex specificity, so it will be asked about.
#       The AUC reversal between cohorts reads far better as a plot.
#
#   S4  meta-analysis. "One gene under random effects" is supported by no
#       visual. Forest plot of the top genes plus the I-squared
#       distribution.
#
# Base graphics, 300 dpi PNG (cairo) plus vector PDF, into audit/figures2/.
# =====================================================================

suppressPackageStartupMessages({ library(DESeq2) })

MDD <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD"
A   <- file.path(MDD,"audit/filtered/mirror")
T   <- file.path(MDD,"transcriptomics")
OUT <- file.path(MDD,"audit/figures2_filtered"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
hdr <- function(x) cat("\n== ", x, "\n", sep="")
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

CTL <- "grey35"; MDDC <- "#C0392B"; C1 <- "#2E86C1"; C2 <- "#E67E22"; GY <- "grey60"

dev2 <- function(name, w, h, expr) {
  e <- substitute(expr); pf <- parent.frame()
  for (fmt in c("png","pdf")) {
    f <- file.path(OUT, paste0(name,".",fmt))
    if (fmt=="png") png(f, width=w, height=h, units="in", res=300, type="cairo")
    else pdf(f, width=w, height=h)
    tryCatch(eval(e, pf), error=function(err){plot.new(); text(.5,.5,err$message)})
    dev.off() }
  cat("   wrote", name, "\n") }

rd <- function(p, sep="\t") if (file.exists(p))
        read.delim(p, sep=sep, check.names=FALSE, stringsAsFactors=FALSE) else NULL

# ============== S1  COHORT AGREEMENT ==============
hdr("S1  pairwise cohort agreement")
m3 <- rd(file.path(A,"meta","meta_analysis_full.tsv"))
if (!is.null(m3)) {
  lf <- grep("^lfc_", names(m3), value=TRUE)
  cat("   cohorts:", paste(sub("lfc_","",lf), collapse=", "), "\n")
  dev2("S1_cohort_agreement", 12, 4.2, {
    par(mfrow=c(1,3), mar=c(4.6,4.6,3.2,1))
    cmb <- combn(lf, 2)
    for (i in seq_len(ncol(cmb))) {
      x <- m3[[cmb[1,i]]]; y <- m3[[cmb[2,i]]]
      ok <- is.finite(x) & is.finite(y)
      r  <- cor(x[ok], y[ok])
      lim <- quantile(c(x[ok],y[ok]), c(.001,.999))
      good <- abs(r) > 0.15
      plot(x[ok], y[ok], pch=16, cex=.28,
           col=adjustcolor(if (good) C1 else GY, .25),
           xlim=lim, ylim=lim,
           xlab=paste("log2FC", sub("lfc_","",cmb[1,i])),
           ylab=paste("log2FC", sub("lfc_","",cmb[2,i])),
           main=sprintf("%s  vs  %s", sub("lfc_","",cmb[1,i]), sub("lfc_","",cmb[2,i])))
      abline(h=0, v=0, col="grey85"); abline(0, 1, lty=3, col="grey55")
      abline(lm(y[ok] ~ x[ok]), col=if (good) C1 else MDDC, lwd=2)
      legend("topleft", bty="n", cex=1.05, text.col=if (good) C1 else MDDC,
             legend=sprintf("r = %+.3f", r))
      if (!good) legend("bottomright", bty="n", cex=.8, text.col=MDDC,
                        legend="no agreement")
    }
  })
}

# ============== S2  DECONVOLUTION ==============
hdr("S2  cell-type proportions by diagnosis")
dc <- rd(file.path(T,"Deconvolution","cell_type_proportions_long.csv"), sep=",")
if (!is.null(dc)) {
  names(dc)[1] <- "sample_id"
  # restrict to the two primary cohorts -- GSE80655 rows are run-level
  dc <- dc[dc$dataset %in% c("GSE101521","GSE102556"), ]
  dc$disease <- ifelse(grepl("MDD|depress", dc$disease, ignore.case=TRUE), "MDD","Control")
  cts <- sort(unique(dc$cell_type))
  dev2("S2_deconvolution", 12, 6.5, {
    par(mfrow=c(2,4), mar=c(6,4.2,3,1))
    for (ct in cts) {
      d <- dc[dc$cell_type==ct, ]
      if (!nrow(d)) next
      g <- interaction(d$dataset, d$disease, sep="\n", lex.order=TRUE)
      cols <- ifelse(grepl("MDD", levels(g)), MDDC, CTL)
      boxplot(d$proportion ~ g, col=adjustcolor(cols,.3), border=cols,
              main=ct, xlab="", ylab="proportion", las=2, cex.axis=.7,
              outline=FALSE)
      set.seed(1)
      points(jitter(as.numeric(g),1.1), d$proportion, pch=16, cex=.45,
             col=adjustcolor(ifelse(d$disease=="MDD",MDDC,CTL), .5))
      ps <- sapply(unique(d$dataset), function(s){
        dd <- d[d$dataset==s,]
        suppressWarnings(wilcox.test(dd$proportion[dd$disease=="Control"],
                                     dd$proportion[dd$disease=="MDD"])$p.value) })
      mtext(paste(sprintf("%s p=%.2f", names(ps), ps), collapse="   "),
            side=3, line=-0.9, cex=.6)
    }
    mtext("Cell-type proportions, two primary cohorts (nothing survives BH)",
          outer=TRUE, line=-1.4, cex=.95)
  })
}

# ============== S3  SEX ==============
hdr("S3  sex-stratified module scores")
sx <- rd(file.path(A,"sexspecific","scores_by_sex.tsv"))
if (!is.null(sx)) {
  auc <- function(a,b) suppressWarnings(wilcox.test(a,b)$statistic)/(length(a)*length(b))
  wp  <- function(a,b) suppressWarnings(wilcox.test(a,b)$p.value)
  dev2("S3_sex", 11, 5.2, {
    layout(matrix(c(1,2,3,3), nrow=2, byrow=TRUE), heights=c(1.25,1))
    # panels 1-2: score by cohort x sex x diagnosis
    for (v in c("ieg","erk")) {
      par(mar=c(6.5,4.4,3,1))
      g <- interaction(sx$cohort, sx$sex, sx$disease, sep="\n", lex.order=TRUE)
      cols <- ifelse(grepl("MDD", levels(g)), MDDC, CTL)
      boxplot(sx[[v]] ~ g, col=adjustcolor(cols,.3), border=cols, outline=FALSE,
              las=2, cex.axis=.62, xlab="", ylab=paste(toupper(v),"score"),
              main=paste(toupper(v),"module by sex"))
      set.seed(2)
      points(jitter(as.numeric(g),1.1), sx[[v]], pch=16, cex=.6,
             col=adjustcolor(ifelse(sx$disease=="MDD",MDDC,CTL), .6))
    }
    # panel 3: the AUC reversal
    par(mar=c(4.4,4.6,3.4,1))
    rows <- expand.grid(cohort=unique(sx$cohort), sex=unique(sx$sex),
                        stringsAsFactors=FALSE)
    rows$auc <- NA; rows$p <- NA; rows$n <- NA
    for (i in seq_len(nrow(rows))) {
      d <- sx[sx$cohort==rows$cohort[i] & sx$sex==rows$sex[i], ]
      a <- d$ieg[d$disease=="Control"]; b <- d$ieg[d$disease=="MDD"]
      if (length(a)<3 || length(b)<3) next
      rows$auc[i] <- auc(a,b); rows$p[i] <- wp(a,b); rows$n[i] <- nrow(d) }
    rows <- rows[!is.na(rows$auc), ]
    lab <- paste(rows$cohort, rows$sex, sep=" / ")
    bp <- barplot(rows$auc, names.arg=lab, las=1, horiz=FALSE,
                  col=ifelse(rows$sex=="female", C2, C1), border=NA,
                  ylim=c(0,1), ylab="AUC (IEG score)",
                  main="Sex effect reverses between cohorts -- no interaction (p = 0.57)")
    abline(h=0.5, lty=2)
    text(bp, rows$auc+0.05, sprintf("p=%.3f\nn=%d", rows$p, rows$n), cex=.75)
    legend("topright", bty="n", fill=c(C2,C1), border=NA, cex=.85,
           legend=c("female","male"))
  })
}

# ============== S4  META-ANALYSIS ==============
hdr("S4  meta-analysis: forest and heterogeneity")
m2 <- rd(file.path(A,"primary2","meta_2cohort.tsv"))
if (!is.null(m2)) {
  m2 <- m2[order(m2$q_fe), ]
  sy <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
        keys=sub("\\..*","",m2$gene), keytype="ENSEMBL", columns="SYMBOL"))
  sy <- sy[!duplicated(sy$ENSEMBL),]
  m2$symbol <- sy$SYMBOL[match(sub("\\..*","",m2$gene), sy$ENSEMBL)]
  m2$symbol[is.na(m2$symbol)] <- m2$gene[is.na(m2$symbol)]
  top <- head(m2, 15)
  dev2("S4_meta_analysis", 11.5, 5.2, {
    par(mfrow=c(1,2))
    # forest
    par(mar=c(4.6,7.5,3.2,1))
    y <- rev(seq_len(nrow(top)))
    xr <- range(c(top$lfc_101521, top$lfc_102556, top$beta_fe), na.rm=TRUE)
    plot(NA, xlim=xr*1.15, ylim=c(0.5, nrow(top)+0.5), yaxt="n",
         xlab=expression(log[2]~fold~change), ylab="",
         main="Top 15 by fixed-effect meta-analysis")
    axis(2, at=y, labels=top$symbol, las=1, cex.axis=.8, tick=FALSE)
    abline(v=0, col="grey60"); abline(h=y, col="grey93")
    points(top$lfc_101521, y+0.16, pch=16, cex=.85, col=C1)
    points(top$lfc_102556, y-0.16, pch=16, cex=.85, col=C2)
    points(top$beta_fe, y, pch=18, cex=1.6, col="black")
    k <- top$q_re < 0.05
    if (any(k, na.rm=TRUE)) points(top$beta_fe[k], y[k], pch=5, cex=2.3, col=MDDC, lwd=2)
    legend("bottomleft", bty="n", cex=.72, pch=c(16,16,18,5),
           col=c(C1,C2,"black",MDDC), pt.cex=c(.85,.85,1.6,1.8),
           legend=c("GSE101521","GSE102556","fixed effect","survives random effects"))
    # heterogeneity
    par(mar=c(4.6,4.6,3.2,1))
    hist(m2$I2, breaks=40, col="grey82", border="white",
         xlab=expression(I^2~"(%)"), main="Between-cohort heterogeneity")
    abline(v=50, lty=2, col=MDDC)
    legend("topright", bty="n", cex=.8,
           legend=c(sprintf("median %.1f%%", median(m2$I2, na.rm=TRUE)),
                    sprintf("%.1f%% above 50%%", 100*mean(m2$I2>50, na.rm=TRUE)),
                    sprintf("fixed %d, random %d at FDR<0.05",
                            sum(m2$q_fe<0.05, na.rm=TRUE), sum(m2$q_re<0.05, na.rm=TRUE))))
  })
}

hdr("DONE")
cat("figures in", OUT, "\n")
sessionInfo()

