
#!/usr/bin/env Rscript
# =====================================================================
# FIGURES for the two-cohort primary analysis.
#
# The eight existing figures are built on the 167-column three-cohort
# matrix, superseded twice. These replace them and add what was missing
# against the SCZ arm's figure list: MA plots, a heatmap, and panels for
# the new results (ERK module, robustness, mode of death).
#
# Base graphics throughout, to avoid depending on packages that may not
# be in the module stack. 300 dpi PNG plus vector PDF, matching the
# existing convention.
#
# Inputs, all produced by earlier scripts:
#   audit/exploratory/pca_*.tsv
#   audit/primary2/{pooled_2cohort,meta_2cohort,IEG_2cohort}.tsv
#   audit/primary2/camera_pooled2.csv
#   audit/robust/{leave_one_gene_out,leave_one_donor_out,empirical_null}.tsv
#   audit/suicide/suicide_decomposition.tsv
#   audit/sexspecific/scores_by_sex.tsv
#
# Writes to MDD/audit/figures2/. Existing figures are NOT touched.
# =====================================================================

suppressPackageStartupMessages({ library(DESeq2) })

MDD <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD"
A   <- file.path(MDD,"audit/filtered/mirror")
OUT <- file.path(MDD,"audit/figures2_filtered"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
hdr <- function(x) cat("\n== ", x, "\n", sep="")

IEG <- c("NPAS4","FOS","FOSB","EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3",
         "ATF3","BTG2","ARC","JUN","JUNB","DUSP1","DUSP6")
CTL <- "grey35"; MDDC <- "#C0392B"; C1 <- "#2E86C1"; C2 <- "#E67E22"

dev2 <- function(name, w=9, h=6, expr) {
  e <- substitute(expr); pf <- parent.frame()
  for (fmt in c("png","pdf")) {
    f <- file.path(OUT, paste0(name, ".", fmt))
    if (fmt=="png") png(f, width=w, height=h, units="in", res=300, type="cairo")
    else pdf(f, width=w, height=h)
    tryCatch(eval(e, pf), error=function(e) { plot.new(); text(.5,.5,paste("failed:",e$message)) })
    dev.off()
  }
  cat("   wrote", name, ".png/.pdf\n") }

rd <- function(p, sep="\t") if (file.exists(p))
        read.delim(p, sep=sep, check.names=FALSE, stringsAsFactors=FALSE) else NULL

# ===================== F1  PCA =====================
hdr("F1  PCA, per cohort and pooled")
dev2("F1_PCA", 12, 4.2, {
  par(mfrow=c(1,3), mar=c(4.2,4.2,3,1), cex.axis=.9)
  for (s in c("GSE101521","GSE102556","POOLED")) {
    d <- rd(file.path(A,"exploratory",paste0("pca_",s,".tsv")))
    if (is.null(d)) { plot.new(); title(paste(s,"missing")); next }
    cols <- ifelse(d$disease=="MDD", MDDC, CTL)
    pch  <- if ("study" %in% names(d)) ifelse(d$study=="GSE101521",16,17) else 16
    plot(d$PC1, d$PC2, col=cols, pch=pch, cex=1.3, xlab="PC1", ylab="PC2",
         main=paste0(s, " (n=", nrow(d), ")"))
    abline(h=0, v=0, col="grey85", lty=3)
    legend("topright", bty="n", cex=.85,
           legend=c("Control","MDD"), col=c(CTL,MDDC), pch=16)
    if (s=="POOLED") legend("bottomright", bty="n", cex=.8,
           legend=c("GSE101521","GSE102556"), pch=c(16,17), col="grey40")
  }
})

# ===================== F2  volcano + MA =====================
hdr("F2  volcano and MA, pooled two-cohort")
res <- rd(file.path(A,"primary2","pooled_2cohort.tsv"))
if (!is.null(res)) {
  res$sig <- !is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 0.3
  sy <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
        keys=sub("\\..*","",res$gene), keytype="ENSEMBL", columns="SYMBOL"))
  sy <- sy[!duplicated(sy$ENSEMBL),]
  res$symbol <- sy$SYMBOL[match(sub("\\..*","",res$gene), sy$ENSEMBL)]
  res$isieg <- res$symbol %in% IEG

  dev2("F2_volcano_MA", 11, 5.2, {
    par(mfrow=c(1,2), mar=c(4.4,4.4,3,1))
    # volcano
    y <- -log10(res$pvalue)
    plot(res$log2FoldChange, y, pch=16, cex=.35, col="grey78",
         xlab=expression(log[2]~fold~change), ylab=expression(-log[10]~p),
         main="Pooled two-cohort DE (107 donors)")
    points(res$log2FoldChange[res$sig], y[res$sig], pch=16, cex=.5, col="grey30")
    k <- res$isieg & !is.na(res$pvalue)
    points(res$log2FoldChange[k], y[k], pch=21, cex=1.5, bg=MDDC, col="black", lwd=.7)
    abline(v=c(-0.3,0.3), lty=3, col="grey60")
    lab <- res[k & y > 3.5, ]
    if (nrow(lab)) text(lab$log2FoldChange, -log10(lab$pvalue), lab$symbol,
                        pos=2, cex=.75, col=MDDC)
    legend("topleft", bty="n", cex=.85, pch=c(16,21), pt.bg=c(NA,MDDC),
           col=c("grey30","black"), legend=c("FDR<0.05 & |LFC|>0.3","IEG"))
    # MA
    plot(log10(res$baseMean+1), res$log2FoldChange, pch=16, cex=.35, col="grey78",
         xlab=expression(log[10]~mean~expression), ylab=expression(log[2]~fold~change),
         main="MA plot")
    points(log10(res$baseMean[res$sig]+1), res$log2FoldChange[res$sig],
           pch=16, cex=.5, col="grey30")
    points(log10(res$baseMean[k]+1), res$log2FoldChange[k],
           pch=21, cex=1.5, bg=MDDC, col="black", lwd=.7)
    abline(h=0, col="grey50"); abline(h=c(-0.3,0.3), lty=3, col="grey60")
  })
}

# ===================== F3  IEG heatmap =====================
hdr("F3  IEG expression heatmap")
CO <- list(GSE101521=list(cnt=file.path(MDD,"transcriptomics/GSE101521/counts/GSE101521_gene_counts_matrix.tsv"),
                          col=file.path(MDD,"transcriptomics/DE/GSE101521_colData.tsv")),
           GSE102556=list(cnt=file.path(MDD,"transcriptomics/GSE102556/counts/GSE102556_gene_counts_matrix.tsv"),
                          col=file.path(MDD,"transcriptomics/DE/GSE102556_colData.tsv")))
Z <- NULL; ann <- NULL
for (n in names(CO)) {
  cn <- read.delim(CO[[n]]$cnt, row.names=1, check.names=FALSE)
  cl <- read.delim(CO[[n]]$col, check.names=FALSE); names(cl)[1] <- "sample_id"
  rownames(cl) <- cl$sample_id
  k <- intersect(colnames(cn), rownames(cl)); cn <- round(as.matrix(cn[,k])); cl <- cl[k,]
  cl$disease <- ifelse(grepl("MDD|depress", cl$disease, ignore.case=TRUE),"MDD","Control")
  keep <- rowSums(cn>=10) >= min(table(cl$disease))
  dds <- estimateSizeFactors(DESeqDataSetFromMatrix(cn[keep,], cl, design=~1))
  lc  <- log2(counts(dds, normalized=TRUE)+1)
  m2  <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
          keys=sub("\\..*","",rownames(lc)), keytype="ENSEMBL", columns="SYMBOL"))
  m2 <- m2[!duplicated(m2$ENSEMBL),]
  s2 <- m2$SYMBOL[match(sub("\\..*","",rownames(lc)), m2$ENSEMBL)]
  idx <- which(s2 %in% IEG)
  sub <- lc[idx,,drop=FALSE]; rownames(sub) <- s2[idx]
  sub <- t(scale(t(sub)))                       # z-score within cohort
  o   <- order(cl$disease)
  Z   <- if (is.null(Z)) sub[,o,drop=FALSE] else
           cbind(Z, sub[match(rownames(Z), rownames(sub)), o, drop=FALSE])
  ann <- rbind(ann, data.frame(cohort=n, disease=cl$disease[o]))
}
if (!is.null(Z)) {
  dev2("F3_IEG_heatmap", 11, 5.5, {
    layout(matrix(c(1,2), nrow=2), heights=c(1,7))
    par(mar=c(0.3,7,2.5,6))
    image(seq_len(ncol(Z)), 1, matrix(ifelse(ann$disease=="MDD",1,0), ncol=1),
          col=c(CTL,MDDC), axes=FALSE, xlab="", ylab="")
    b <- which(diff(as.numeric(factor(ann$cohort)))!=0)
    if (length(b)) abline(v=b+0.5, lwd=2)
    title("Immediate-early genes, z-scored within cohort", cex.main=1.1)
    mtext("Control", side=2, las=1, line=.4, cex=.7, col=CTL)
    par(mar=c(4,7,0.4,6))
    zz <- pmax(pmin(Z, 2.5), -2.5)
    image(seq_len(ncol(zz)), seq_len(nrow(zz)), t(zz),
          col=colorRampPalette(c("#2166AC","white","#B2182B"))(64),
          axes=FALSE, xlab="donors (Control | MDD, within cohort)", ylab="")
    axis(2, at=seq_len(nrow(zz)), labels=rownames(zz), las=1, tick=FALSE, cex.axis=.8)
    if (length(b)) abline(v=b+0.5, lwd=2)
    box()
  })
}

# ===================== F4  module scores =====================
hdr("F4  IEG and ERK module scores")
sc <- rd(file.path(A,"sexspecific","scores_by_sex.tsv"))
if (!is.null(sc)) {
  dev2("F4_module_scores", 10, 5, {
    par(mfrow=c(1,2), mar=c(5,4.4,3,1))
    for (v in c("ieg","erk")) {
      grp <- interaction(sc$cohort, sc$disease, sep="\n", lex.order=TRUE)
      cols <- ifelse(grepl("MDD", levels(grp)), MDDC, CTL)
      boxplot(sc[[v]] ~ grp, col=adjustcolor(cols,.35), border=cols,
              ylab=paste(toupper(v),"score (median percentile rank)"),
              xlab="", main=paste(toupper(v),"module"), outline=FALSE, las=2, cex.axis=.8)
      set.seed(1)
      points(jitter(as.numeric(grp), 1.1), sc[[v]],
             pch=16, cex=.7, col=adjustcolor(ifelse(sc$disease=="MDD",MDDC,CTL), .65))
      for (ch in unique(sc$cohort)) {
        d <- sc[sc$cohort==ch,]
        p <- wilcox.test(d[[v]][d$disease=="Control"], d[[v]][d$disease=="MDD"])$p.value
        mtext(sprintf("%s p=%.3f", ch, p), side=3, line=-1.2,
              at=mean(which(grepl(ch, levels(grp)))), cex=.75)
      }
    }
  })
}

# ===================== F5  CAMERA =====================
hdr("F5  set-level testing, corrected vs uncorrected")
cam <- rd(file.path(A,"primary2","camera_pooled2.csv"), sep=",")
if (!is.null(cam)) {
  nm <- if ("set" %in% names(cam)) cam$set else rownames(cam)
  dev2("F5_CAMERA", 8.5, 5.5, {
    par(mar=c(5,15,3,2))
    p <- cam$PValue; o <- order(p, decreasing=TRUE)
    cols <- ifelse(grepl("NEG", nm[o]), "grey55",
             ifelse(nm[o]=="IEG", MDDC, C1))
    bp <- barplot(-log10(p[o]), horiz=TRUE, col=cols, border=NA,
                  names.arg=nm[o], las=1, cex.names=.85,
                  xlab=expression(-log[10]~p~(CAMERA)),
                  main="Set-level testing, pooled two-cohort")
    abline(v=-log10(0.05), lty=2, col="black")
    text(-log10(0.05), max(bp), " p=0.05", pos=4, cex=.75, xpd=NA)
    legend("bottomright", bty="n", cex=.8, fill=c(MDDC,C1,"grey55"), border=NA,
           legend=c("IEG","other sets","negative controls"))
  })
}

# ===================== F6  robustness =====================
hdr("F6  robustness panel")
logo <- rd(file.path(A,"robust","leave_one_gene_out.tsv"))
lodo <- rd(file.path(A,"robust","leave_one_donor_out.tsv"))
enul <- rd(file.path(A,"robust","empirical_null.tsv"))
dev2("F6_robustness", 12, 4.2, {
  par(mfrow=c(1,3), mar=c(4.5,4.5,3,1))
  base <- 3.577e-03
  if (!is.null(logo)) {
    o <- order(logo$p)
    barplot(-log10(logo$p[o]), names.arg=logo$dropped[o], las=2, cex.names=.7,
            col=C1, border=NA, ylab=expression(-log[10]~p),
            main="Leave-one-gene-out")
    abline(h=-log10(0.05), lty=2); abline(h=-log10(base), lty=3, col=MDDC)
  }
  if (!is.null(lodo)) {
    hist(-log10(lodo$p), breaks=15, col=C1, border="white",
         xlab=expression(-log[10]~p), main="Leave-one-donor-out (50 refits)")
    abline(v=-log10(0.05), lty=2); abline(v=-log10(base), lty=3, col=MDDC, lwd=2)
  }
  if (!is.null(enul)) {
    np <- enul$p[grepl("^rand", enul$set)]; ob <- enul$p[enul$set=="IEG"]
    hist(-log10(np), breaks=40, col="grey80", border="white",
         xlab=expression(-log[10]~p), main="Empirical null (1000 matched sets)")
    abline(v=-log10(ob), col=MDDC, lwd=2.5)
    text(-log10(ob), par("usr")[4]*.85, " observed IEG", pos=2, col=MDDC, cex=.85)
  }
})

# ===================== F7  mode of death =====================
hdr("F7  mode of death decomposition")
sd_ <- rd(file.path(A,"suicide","suicide_decomposition_bothcohorts.tsv"))
if (is.null(sd_)) sd_ <- rd(file.path(A,"suicide","suicide_decomposition.tsv"))
if (!is.null(sd_)) {
  i <- sd_[sd_$set=="IEG",]
  dev2("F7_mode_of_death", 10, 6, {
    par(mar=c(13,4.6,3,2))
    v <- i$ieg_median_logFC; o <- order(v)
    cols <- ifelse(v[o] > 0, C1, MDDC)
    lab <- if ("test" %in% names(i)) paste0(substr(i$test[o],4,9),"  ",i$contrast[o]) else i$contrast[o]
    bp <- barplot(v[o], names.arg=lab, las=2, cex.names=.68,
                  col=cols, border=NA, ylab=expression(median~log[2]~FC~(IEG~set)),
                  main="IEG module by contrast")
    abline(h=0)
    mtext("cohort prefix: 101521 = GSE101521, 102556 = GSE102556",
          side=1, line=11.5, cex=.62, col="grey35")
    text(bp, rep(0.015, length(v)),
         sprintf("p=%.3f\n%d/16", i$p_corr[o], i$ieg_negative[o]),
         cex=.6, pos=3, offset=0, xpd=NA)
         cex=.6, pos=ifelse(v[o]>0,3,1), offset=.35, xpd=NA)
  })
}

hdr("DONE")
cat("figures in", OUT, "\n")
cat("\nStill to add (needs new analysis, not just plotting):\n")
cat("  - three-way enrichment panel (all / up-only / down-only ORA),\n")
cat("    as suggested for the ASD arm\n")
sessionInfo()

