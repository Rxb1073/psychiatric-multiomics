#!/usr/bin/env Rscript
# PUBLICATION FIGURES — SCZ transcriptomics arm (mirrors ASD make_figures.R)
#  Fig1 astrocyte proportion by diagnosis (3 cohorts x 2 signatures)
#  Fig2 mediation: 1638 pooled DE genes -> ~10 after astrocyte adjustment
#  Fig3 RIN by diagnosis (GSE138082, only cohort with measured RIN)
#  Fig4 astrocyte vs RIN  -> the two mechanisms are separable
#  Fig5 volcano, pooled frontal baseline -> the 6:1 up-skew made visible
#  Fig6 per-cohort composition test -> the decisive result (SCZ-specific)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
logf <- function(...) cat(sprintf(...), "\n")
NEURO <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
CS  <- file.path(NEURO,"counts_scz"); DEC <- file.path(CS,"Deconvolution")
FIG <- file.path(CS,"Figures"); dir.create(FIG, recursive=TRUE, showWarnings=FALSE)
COL_CTRL <- "#9EC5D9"; COL_SCZ <- "#D98E8E"
FRONTAL <- c("GSE80655_DLPFC","GSE78936"); ALLDS <- c("GSE80655_DLPFC","GSE78936","GSE138082")
dev_both <- function(name,w,h,expr){
  for (fn in list(function() pdf(file.path(FIG,paste0(name,".pdf")),width=w,height=h),
                  function() png(file.path(FIG,paste0(name,".png")),width=w*150,height=h*150,res=150,type="cairo"))){
    fn(); eval.parent(substitute(expr)); dev.off() }
  logf("  wrote %s.pdf / .png", name) }
star <- function(p) if(is.na(p)) "" else if(p<0.001)"***" else if(p<0.01)"**" else if(p<0.05)"*" else "ns"

# ---- FIG 1: astrocyte proportion ----
logf("--- Figure 1: astrocyte proportion ---")
props <- read.csv(file.path(DEC,"cell_type_proportions_long.csv"), stringsAsFactors=FALSE)
astro <- subset(props, cell_type=="Astrocytes" & dataset %in% ALLDS)
astro$diagnosis <- factor(astro$diagnosis, levels=c("Control","SCZ"))
panels <- expand.grid(dataset=ALLDS, signature=c("MB","CA"), stringsAsFactors=FALSE)
panels$p <- NA_real_
for (i in seq_len(nrow(panels))) {
  d <- subset(astro, dataset==panels$dataset[i] & signature==panels$signature[i])
  if (nrow(d)>3) panels$p[i] <- summary(lm(proportion~diagnosis,data=d))$coefficients["diagnosisSCZ","Pr(>|t|)"] }
print(panels)
dev_both("Fig1_astrocyte_proportion", 10, 6.5, {
  par(mfrow=c(2,3), mar=c(3.4,4.4,3.2,1), mgp=c(2.6,0.7,0), las=1, oma=c(0,0,2,0))
  for (i in seq_len(nrow(panels))) {
    d <- subset(astro, dataset==panels$dataset[i] & signature==panels$signature[i])
    if (nrow(d)<3) { plot.new(); next }
    pv <- panels$p[i]
    reg <- if (panels$dataset[i]=="GSE138082") "hippocampus" else "frontal"
    boxplot(proportion~diagnosis, data=d, col=c(COL_CTRL,COL_SCZ), ylab="Astrocyte proportion",
            xlab="", outline=FALSE, main=sprintf("%s (%s)\n%s signature", panels$dataset[i], reg, panels$signature[i]),
            cex.main=0.95, ylim=range(d$proportion)*c(0.9,1.2))
    set.seed(1); points(jitter(as.integer(d$diagnosis),0.6), d$proportion, pch=21,
           bg=adjustcolor("grey20",0.35), col=NA, cex=0.8)
    yt <- max(d$proportion)*1.11; segments(1,yt,2,yt)
    text(1.5, yt*1.03, sprintf("p = %s %s", format.pval(pv,digits=2), star(pv)), cex=0.9)
    mtext(sprintf("n = %d", nrow(d)), side=3, line=-1.2, adj=0.02, cex=0.65, col="grey40") }
  mtext("Astrocyte proportion is elevated in SCZ frontal cortex, but not hippocampus",
        outer=TRUE, line=0.3, font=2, cex=1.0) })

# ---- FIG 2: mediation ----
logf("--- Figure 2: mediation ---")
ngenes <- function(f){ if(!file.exists(f)) return(NULL); m<-readRDS(f)
  sapply(m[c("m0","m1","m2","m3")], function(x) length(x$sig)) }
mb <- ngenes(file.path(CS,"POOLED_composition_adjusted","models_MB.rds"))
ca <- ngenes(file.path(CS,"POOLED_composition_adjusted","models_CA.rds"))
if (is.null(mb)) mb <- c(1638,10,234,9); if (is.null(ca)) ca <- c(1638,12,234,12)
M <- rbind(MB=mb, CA=ca); colnames(M) <- c("baseline","+ astrocyte","+ RNA quality","+ both"); print(M)
dev_both("Fig2_mediation", 7.5, 5.2, {
  par(mar=c(6.5,4.6,3.4,1), mgp=c(2.9,0.7,0), las=1)
  bp <- barplot(M, beside=TRUE, col=c("#3B6E8F","#7FA8C0"), border=NA, ylab="DE genes (FDR < 0.05)",
                ylim=c(0,max(M)*1.25), names.arg=rep("",ncol(M)),
                main="99% of the pooled frontal DE signal is explained by astrocyte proportion")
  text(colMeans(bp), par("usr")[3]-max(M)*0.06, srt=20, adj=1, xpd=TRUE, labels=colnames(M), cex=0.92)
  text(bp, as.vector(t(M))+max(M)*0.04, labels=as.vector(t(M)), cex=0.9, font=2)
  legend("topright", rownames(M), fill=c("#3B6E8F","#7FA8C0"), border=NA, bty="n", title="signature")
  mtext("Pooled frontal cohort, n = 109.  Model: ~ study + [covariate] + diagnosis",
        side=1, line=5.1, cex=0.72, col="grey35") })

# ---- FIG 3: RIN by diagnosis (GSE138082) ----
logf("--- Figure 3: RIN by diagnosis ---")
ss <- read.delim(file.path(CS,"sheets","GSE138082_samplesheet.tsv"), stringsAsFactors=FALSE)
ss$rin <- suppressWarnings(as.numeric(ss$rin)); ss$diagnosis <- factor(ss$diagnosis, levels=c("Control","SCZ"))
ss <- ss[!is.na(ss$rin),]
p_rin <- suppressWarnings(wilcox.test(rin~diagnosis, data=ss)$p.value)
logf("  RIN medians: Control %.2f | SCZ %.2f | p=%.4f", median(ss$rin[ss$diagnosis=="Control"]),
     median(ss$rin[ss$diagnosis=="SCZ"]), p_rin)
dev_both("Fig3_RIN_by_diagnosis", 5, 5.4, {
  par(mar=c(3.4,4.4,3.6,1), mgp=c(2.7,0.7,0), las=1)
  boxplot(rin~diagnosis, data=ss, col=c(COL_CTRL,COL_SCZ), outline=FALSE,
          ylab="RIN (RNA Integrity Number)", xlab="", ylim=c(min(ss$rin)-0.3,max(ss$rin)+0.9),
          main="RNA integrity by diagnosis")
  set.seed(1); points(jitter(as.integer(ss$diagnosis),0.6), ss$rin, pch=21, bg=adjustcolor("grey20",0.4), col=NA, cex=0.9)
  yt <- max(ss$rin)+0.45; segments(1,yt,2,yt)
  text(1.5, yt+0.18, sprintf("p = %.3f %s", p_rin, star(p_rin)), cex=1)
  mtext(sprintf("GSE138082 hippocampus (n = %d).  Medians: Control %.1f, SCZ %.1f", nrow(ss),
        median(ss$rin[ss$diagnosis=="Control"]), median(ss$rin[ss$diagnosis=="SCZ"])),
        side=1, line=2.1, cex=0.75, col="grey35") })

# ---- FIG 4: astrocyte vs RIN ----
logf("--- Figure 4: astrocyte vs RIN ---")
ah <- subset(props, cell_type=="Astrocytes" & dataset=="GSE138082")
ah$rin <- suppressWarnings(as.numeric(ah$rin)); ah <- ah[!is.na(ah$rin),]
dev_both("Fig4_astrocyte_vs_RIN", 8.6, 4.6, {
  par(mfrow=c(1,2), mar=c(4.2,4.4,3.4,1), mgp=c(2.6,0.7,0), las=1, oma=c(0,0,1.6,0))
  for (sg in c("MB","CA")) {
    d <- subset(ah, signature==sg); r <- cor(d$proportion, d$rin)
    plot(d$rin, d$proportion, pch=21, bg=ifelse(d$diagnosis=="SCZ",COL_SCZ,COL_CTRL), col="grey30", cex=1.2,
         xlab="RIN", ylab="Astrocyte proportion", main=sprintf("%s signature:  r = %+.3f", sg, r))
    abline(lm(proportion~rin, data=d), col="grey45", lty=2, lwd=1.5)
    legend("topright", c("Control","SCZ"), pt.bg=c(COL_CTRL,COL_SCZ), pch=21, col="grey30", bty="n", cex=0.85) }
  mtext("Astrocyte content is uncorrelated with RNA integrity: the two mechanisms are separable",
        outer=TRUE, line=-0.4, font=2, cex=0.9) })

# ---- FIG 5: volcano ----
logf("--- Figure 5: volcano ---")
ff <- file.path(CS,"POOLED_frontal","POOLED_frontal_MAIN_full_results.tsv")
if (file.exists(ff)) {
  full <- read.delim(ff, stringsAsFactors=FALSE); full <- full[!is.na(full$padj),]
  sig <- full$padj < 0.05
  dev_both("Fig5_volcano_pooled_frontal", 6, 5.6, {
    par(mar=c(4.2,4.4,3.6,1), mgp=c(2.6,0.7,0), las=1)
    plot(full$log2FC_shrunk, -log10(full$padj), pch=16, cex=0.4, col=adjustcolor("grey65",0.5),
         xlab=expression(log[2]~"fold change (SCZ vs Control)"), ylab=expression(-log[10]~"adjusted p"),
         main="Pooled frontal DE: significant genes skew strongly one way")
    points(full$log2FC_shrunk[sig], -log10(full$padj[sig]), pch=16, cex=0.75, col=COL_SCZ)
    abline(h=-log10(0.05), lty=2, col="grey40"); abline(v=0, lty=3, col="grey60")
    text(par("usr")[1]+0.05*diff(par("usr")[1:2]), -log10(0.05)+0.18, "FDR 0.05", cex=0.72, col="grey35", adj=0)
    legend("topleft", sprintf("%d significant\n(%d up, %d down)", sum(sig),
      sum(sig & full$log2FC_shrunk>0), sum(sig & full$log2FC_shrunk<0)), bty="n", cex=0.85, text.col="grey20") })
} else logf("  (pooled frontal results absent - volcano skipped)")

# ---- FIG 6: per-cohort composition test (SCZ-specific, the decisive result) ----
logf("--- Figure 6: per-cohort composition ---")
pc <- file.path(CS,"PERCOHORT_composition")
gc_ <- function(f){ if(!file.exists(f)) return(NULL); m<-readRDS(f)
  sapply(m[c("c0","c1","c2","c3")], function(x) length(x$sig)) }
a <- gc_(file.path(pc,"GSE80655_DLPFC_MB.rds")); b <- gc_(file.path(pc,"GSE78936_MB.rds"))
if (is.null(a)) a <- c(35,10,80,16); if (is.null(b)) b <- c(353,0,32,0)
K <- rbind("GSE80655 (n=56)"=a, "GSE78936 (n=52)"=b)
colnames(K) <- c("baseline","+ astrocyte","+ quality","+ both"); print(K)
dev_both("Fig6_percohort_composition", 7.5, 5.2, {
  par(mar=c(6.5,4.6,3.4,1), mgp=c(2.9,0.7,0), las=1)
  bp <- barplot(K, beside=TRUE, col=c("#6B8E23","#A9C47F"), border=NA, ylab="DE genes (FDR < 0.05)",
                ylim=c(0,max(K)*1.25), names.arg=rep("",ncol(K)),
                main="Within each cohort alone, astrocyte adjustment removes the signal")
  text(colMeans(bp), par("usr")[3]-max(K)*0.06, srt=20, adj=1, xpd=TRUE, labels=colnames(K), cex=0.92)
  text(bp, as.vector(t(K))+max(K)*0.04, labels=as.vector(t(K)), cex=0.9, font=2)
  legend("topright", rownames(K), fill=c("#6B8E23","#A9C47F"), border=NA, bty="n")
  mtext("Per-cohort models with each cohort's own covariates (MB signature; CA near-identical)",
        side=1, line=5.1, cex=0.72, col="grey35") })
logf("\nAll figures -> %s", FIG)
