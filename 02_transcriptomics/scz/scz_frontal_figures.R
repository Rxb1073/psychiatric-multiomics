## Frontal arm figures - composition, nulls, and the pseudoreplication demonstration
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(ggplot2)})
logf <- function(...) cat(sprintf(...),"\n")
NP  <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
CS  <- file.path(NP,"counts_scz"); SHD <- file.path(CS,"sheets")
PC  <- file.path(CS,"PERCOHORT_prua")
OUT <- file.path(CS,"FRONTAL_final"); dir.create(OUT, showWarnings=FALSE)
sv <- function(p,w=7,h=6){ ggsave(file.path(OUT,paste0(p,".pdf")),width=w,height=h)
                           ggsave(file.path(OUT,paste0(p,".png")),width=w,height=h,dpi=150) }

## ---- F1: astrocyte proportion by diagnosis, both signatures, corrected samples ----
pr <- read.csv(file.path(CS,"Deconvolution","cell_type_proportions_long.csv"), stringsAsFactors=FALSE)
pr <- subset(pr, cell_type=="Astrocytes")
sh <- do.call(rbind, lapply(c("GSE80655_DLPFC","GSE78936"), function(s){
  x <- read.delim(file.path(SHD,paste0(s,"_samplesheet.tsv")), check.names=FALSE)
  data.frame(srr=x$srr, diagnosis=x$diagnosis, study=sub("_DLPFC","",s)) }))
pr <- pr[, setdiff(names(pr), c("diagnosis","study")), drop=FALSE]
m <- merge(pr, sh, by="srr")
## collapse GSE80655 run-splits: average proportion per BioSample
RT <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/GSE80655/metadata/PRJNA319583_SraRunTable.csv"
rt <- read.csv(RT, check.names=FALSE); mp <- setNames(rt$BioSample, rt$Run)
m$bio <- ifelse(is.na(mp[m$srr]), m$srr, mp[m$srr])
m <- aggregate(proportion ~ bio + signature + diagnosis + study, data=m, FUN=mean)
logf("after collapse: %d rows (%d per signature)", nrow(m), nrow(m)/length(unique(m$signature)))
logf("merged: %d rows | cols: %s", nrow(m), paste(names(m), collapse=", "))
stopifnot("diagnosis" %in% names(m), "proportion" %in% names(m))
m$study <- factor(m$study); m$diagnosis <- factor(m$diagnosis, c("Control","SCZ"))
ggplot(m, aes(diagnosis, proportion, fill=diagnosis)) +
  geom_boxplot(outlier.shape=NA, alpha=.7) + geom_jitter(width=.16, size=1.1, alpha=.6) +
  facet_grid(signature ~ study, scales="free_y") + theme_bw() + guides(fill="none") +
  scale_fill_manual(values=c(Control="#5D8AA8", SCZ="#C0392B")) +
  labs(y="astrocyte mRNA fraction", x=NULL,
       title="Astrocyte fraction by diagnosis (frontal cohorts, both signatures)")
sv("Fig_F1_astrocyte", 8, 7)
for (s in unique(m$signature)) { mm <- subset(m, signature==s)
  logf("astro %s : Wilcoxon p = %.4f", s,
       suppressWarnings(wilcox.test(proportion~diagnosis, data=mm)$p.value)) }

## ---- F2: mediation ladder (counts verified from comp_*.out) ----
lad <- data.frame(
  model = factor(rep(c("M0 baseline","M1 +astrocyte","M2 +quality","M3 +both","M4 +composition"),2),
                 levels=c("M0 baseline","M1 +astrocyte","M2 +quality","M3 +both","M4 +composition")),
  sig   = rep(c("MultiBrain","CA"), each=5),
  genes = c(966,1,54,0,1, 966,1,54,5,1))
ggplot(lad, aes(model, genes, fill=sig)) +
  geom_col(position=position_dodge(.75), width=.7) +
  geom_text(aes(label=genes), position=position_dodge(.75), vjust=-.4, size=3.2) +
  scale_fill_manual(values=c(MultiBrain="#2471A3", CA="#7D6608")) + theme_bw() +
  theme(axis.text.x=element_text(angle=20, hjust=1)) +
  labs(y="DE genes (FDR<0.05)", x=NULL, fill="signature",
       title="Frontal pooled (n=100): astrocyte adjustment removes 99.9% of baseline signal")
sv("Fig_F2_mediation", 8, 5.5)

## ---- F3: the internal control - same covariate, opposite effect by region ----
ctr <- data.frame(
  region = factor(rep(c("Frontal (pooled, n=100)","Hippocampus (n=78)"), each=2),
                  levels=c("Frontal (pooled, n=100)","Hippocampus (n=78)")),
  model  = factor(rep(c("baseline","+ astrocyte"),2), levels=c("baseline","+ astrocyte")),
  genes  = c(966, 1, 146, 176))
ggplot(ctr, aes(model, genes, group=region, colour=region)) +
  geom_line(linewidth=1.1) + geom_point(size=3.4) +
  geom_text(aes(label=genes), vjust=-1, size=3.4, show.legend=FALSE) +
  scale_colour_manual(values=c("#C0392B","#1E8449")) + theme_bw() +
  labs(y="DE genes (FDR<0.05)", x=NULL, colour=NULL,
       title="Astrocyte adjustment: removes frontal signal, retains hippocampal")
sv("Fig_F3_region_contrast", 7.5, 5.5)

## ---- F4: null volcanoes for the two corrected frontal cohorts ----
vol <- function(f, lab, tag) {
  if (!file.exists(f)) { logf("missing %s", f); return(invisible(NULL)) }
  r <- read.csv(f)
  r$tier <- ifelse(is.na(r$padj)|r$padj>=0.05,"NS", ifelse(r$log2FoldChange>0,"UP","DOWN"))
  n <- sum(r$tier!="NS")
  ggplot(r, aes(log2FoldChange, -log10(pvalue), colour=tier)) +
    geom_point(alpha=.5, size=1) +
    scale_colour_manual(values=c(NS="grey80", UP="#C0392B", DOWN="#2471A3")) +
    geom_vline(xintercept=c(-1,1), linetype=2, colour="grey40") + theme_bw() +
    labs(title=sprintf("%s - %d genes FDR<0.05", lab, n))
  sv(tag)
  logf("%s : %d significant", lab, n)
}
vol(file.path(PC,"GSE80655_D2_full.csv"),      "GSE80655 DLPFC (n=48, full model)", "Fig_F4a_volcano_GSE80655")
vol(file.path(PC,"GSE78936_BA11_F1_full.csv"), "GSE78936 BA11 (n=28, full model)",  "Fig_F4b_volcano_BA11")

## ---- F5: pseudoreplication demonstration, 545 -> 0 ----
a <- file.path(PC,"GSE78936_E2_full.csv"); b <- file.path(PC,"GSE78936_BA11_F1_full.csv")
if (file.exists(a) && file.exists(b)) {
  x <- read.csv(a); y <- read.csv(b)
  z <- merge(x[,c("gene","log2FoldChange","padj")], y[,c("gene","log2FoldChange","padj")],
             by="gene", suffixes=c(".all52",".ba11"))
  z$grp <- ifelse(!is.na(z$padj.all52) & z$padj.all52<0.05, "sig with all 52 samples","not significant")
  r <- cor(z$log2FoldChange.all52, z$log2FoldChange.ba11, use="complete.obs")
  ggplot(z, aes(log2FoldChange.all52, log2FoldChange.ba11, colour=grp)) +
    geom_point(alpha=.45, size=1) + geom_abline(slope=1, intercept=0, linetype=2) +
    scale_colour_manual(values=c("not significant"="grey78","sig with all 52 samples"="#C0392B")) +
    theme_bw() + labs(colour=NULL, x="log2FC, all 52 samples (12 donors duplicated)",
      y="log2FC, BA11 only (28 unique donors)",
      title=sprintf("Pseudoreplication: 545 significant -> 0, yet fold-changes agree (r = %.3f)", r))
  sv("Fig_F5_pseudoreplication", 8, 6)
  logf("pseudoreplication figure: r = %.3f | sig all52 = %d | sig BA11 = %d", r,
       sum(z$padj.all52<0.05, na.rm=TRUE), sum(z$padj.ba11<0.05, na.rm=TRUE))
}
logf("\nOutputs -> %s", OUT)
