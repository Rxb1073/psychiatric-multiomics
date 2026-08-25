## Final SCZ figures: IEG finding, meta-analysis, composition, donor-structure demos
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(ggplot2)})
logf <- function(...) cat(sprintf(...),"\n")
CS  <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/counts_scz"
RT  <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/GSE80655/metadata/PRJNA319583_SraRunTable.csv"
OUT <- file.path(CS,"FIGURES_FINAL"); dir.create(OUT, showWarnings=FALSE)
sv <- function(p,w=7,h=6){ ggsave(file.path(OUT,paste0(p,".pdf")),width=w,height=h)
                           ggsave(file.path(OUT,paste0(p,".png")),width=w,height=h,dpi=150) }
IEG <- c("NPAS4","EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3","FOS","FOSB",
         "ARC","ATF3","DUSP1","DUSP6","JUN","JUNB")

## ---- 1. IEG forest: per-cohort effect sizes + meta ----
m <- read.csv(file.path(CS,"META","SCZ_meta_FE_full.csv"), stringsAsFactors=FALSE)
s <- m[m$symbol %in% IEG,]
fr <- do.call(rbind, lapply(c("CA1","DLPFC","BA11"), function(k)
  data.frame(symbol=s$symbol, cohort=k, beta=s[[paste0("beta.",k)]], se=s[[paste0("se.",k)]])))
fr <- rbind(fr, data.frame(symbol=s$symbol, cohort="META (102 donors)",
                           beta=s$beta_FE, se=s$se_FE))
fr$cohort <- factor(fr$cohort, levels=c("CA1","BA11","DLPFC","META (102 donors)"))
fr$symbol <- factor(fr$symbol, levels=s$symbol[order(s$beta_FE)])
ggplot(fr, aes(beta, symbol, colour=cohort)) +
  geom_vline(xintercept=0, linetype=2, colour="grey50") +
  geom_errorbarh(aes(xmin=beta-1.96*se, xmax=beta+1.96*se), height=0,
                 position=position_dodge(.7), alpha=.75) +
  geom_point(position=position_dodge(.7), size=1.9) +
  scale_colour_manual(values=c(CA1="#7FB3D5", BA11="#F5B041",
                               DLPFC="#5D6D7E", "META (102 donors)"="#922B21")) +
  theme_bw() + labs(x=expression(log[2]~fold~change~(SCZ~vs~control)), y=NULL, colour=NULL,
    title="Immediate-early genes: all 15 downregulated",
    subtitle="Per-cohort estimates with 95% CI, and fixed-effect meta-analysis")
sv("Fig_1_IEG_forest", 8.5, 7)
logf("IEG forest: %d genes, all beta_FE negative = %s", nrow(s), all(s$beta_FE<0))

## ---- 2. IEG robustness: effect size across nested models (from ieg log) ----
rb <- data.frame(
  symbol = rep(c("NR4A2","EGR1","NR4A1","ATF3","NPAS4","FOS","EGR2","DUSP6"), 3),
  model  = factor(rep(c("baseline","+ pH","+ pH + neuronal"), each=8),
                  levels=c("baseline","+ pH","+ pH + neuronal")),
  beta   = c(-0.635,-0.899,-0.835,-0.664,-1.051,-0.654,-0.834,-0.402,
             -0.669,-0.912,-0.904,-0.757,-1.187,-0.723,-0.814,-0.405,
             -0.601,-0.803,-0.823,-0.781,-1.065,-0.735,-0.622,-0.272))
ggplot(rb, aes(model, beta, group=symbol, colour=symbol)) +
  geom_hline(yintercept=0, linetype=2, colour="grey50") +
  geom_line(linewidth=.8) + geom_point(size=2.2) + theme_bw() +
  labs(y=expression(log[2]~fold~change), x=NULL, colour=NULL,
    title="IEG effects are robust to agonal state and cellular composition",
    subtitle="GSE80655, 48 unique donors; adjustment attenuates by ~10%, does not remove")
sv("Fig_2_IEG_robustness", 7.5, 5.5)

## ---- 3. cell composition, unique donors, both signatures ----
pr <- read.csv(file.path(CS,"Deconvolution","cell_type_proportions_long.csv"), stringsAsFactors=FALSE)
pr <- pr[, setdiff(names(pr), c("diagnosis","study"))]
pr <- subset(pr, cell_type %in% c("Astrocytes","Neurons","Excitatory","Inhibitory"))
sh <- do.call(rbind, lapply(c("GSE80655_DLPFC","GSE78936"), function(x){
  d <- read.delim(file.path(CS,"sheets",paste0(x,"_samplesheet.tsv")), check.names=FALSE)
  data.frame(srr=d$srr, diagnosis=d$diagnosis, region=d$brain_region, study=sub("_DLPFC","",x))}))
d <- merge(pr, sh, by="srr")
rt <- read.csv(RT, check.names=FALSE); mp <- setNames(rt$BioSample, rt$Run)
d$bio <- ifelse(is.na(mp[d$srr]), d$srr, mp[d$srr])
d <- aggregate(proportion ~ bio+signature+cell_type+diagnosis+region+study, d, mean)
d <- subset(d, study=="GSE80655" | region=="Brodmann's Area 11")
d$diagnosis <- factor(d$diagnosis, c("Control","SCZ"))
lab <- do.call(rbind, lapply(split(d, list(d$signature,d$cell_type), drop=TRUE), function(x)
  data.frame(signature=x$signature[1], cell_type=x$cell_type[1],
             p=wilcox.test(proportion~diagnosis,data=x)$p.value,
             y=max(x$proportion)*1.04)))
lab$txt <- sprintf("p = %.3f", lab$p)
ggplot(d, aes(diagnosis, proportion, fill=diagnosis)) +
  geom_boxplot(outlier.shape=NA, alpha=.75, width=.6) +
  geom_jitter(width=.15, size=.9, alpha=.5) +
  geom_text(data=lab, aes(x=1.5, y=y, label=txt), inherit.aes=FALSE, size=3) +
  facet_wrap(signature ~ cell_type, scales="free_y", ncol=5) +
  scale_fill_manual(values=c(Control="#5D8AA8", SCZ="#C0392B")) +
  theme_bw() + guides(fill="none") +
  labs(x=NULL, y="estimated mRNA fraction",
    title="Cellular composition, frontal cortex (76 unique donors)",
    subtitle="Neuronal fraction reduced in schizophrenia; astrocyte fraction raised at borderline significance")
sv("Fig_3_composition", 13, 3.8)
print(lab[,c("signature","cell_type","p")])

## ---- 4. donor-structure demonstration, hippocampal subfields ----
dl <- data.frame(
  subset = factor(c("all 78 samples","CA1+CA3 (52)","CA1 (26)","CA3 (26)","DG (26)"),
           levels=c("all 78 samples","CA1+CA3 (52)","CA1 (26)","CA3 (26)","DG (26)")),
  donors = c("13 donors x3","13 donors x2","13 donors x1","13 donors x1","13 donors x1"),
  genes  = c(146, 50, 6, 0, 0))
ggplot(dl, aes(subset, genes)) +
  geom_col(fill="#922B21", width=.65) +
  geom_text(aes(label=genes), vjust=-.4, size=3.6) +
  geom_text(aes(label=donors, y=-6), size=2.9, colour="grey35") +
  theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) +
  labs(x=NULL, y="DE genes (FDR<0.05)",
    title="GSE138082: apparent signal scales with samples per donor, not with donors",
    subtitle="Identical model throughout; 13 controls and 13 cases in every analysis")
sv("Fig_4_subfield_pseudoreplication", 8, 5.5)

## ---- 5. meta-analysis volcano ----
m$tier <- ifelse(m$padj_FE<0.05, "FDR<0.05", "ns")
top <- m[m$padj_FE<0.05,]
ggplot(m, aes(beta_FE, -log10(p_FE), colour=tier)) +
  geom_point(alpha=.45, size=1) +
  scale_colour_manual(values=c("FDR<0.05"="#922B21", ns="grey80")) +
  geom_text(data=top, aes(label=symbol), vjust=-.7, size=3, show.legend=FALSE) +
  geom_vline(xintercept=0, linetype=2, colour="grey50") + theme_bw() +
  labs(x=expression(meta~log[2]~fold~change), colour=NULL,
    title=sprintf("Fixed-effect meta-analysis, 102 donors across three regions (%d genes FDR<0.05)",
                  sum(m$padj_FE<0.05, na.rm=TRUE)))
sv("Fig_5_meta_volcano", 8, 6)

## ---- 6. summary of all clean analyses ----
sm <- data.frame(
  analysis = factor(c("Frontal pool (76)","GSE80655 DLPFC (48)","GSE78936 BA11 (28)",
                      "Hippocampus CA1 (26)","Hippocampus CA3 (26)","Hippocampus DG (26)",
                      "Meta-analysis (102)"),
             levels=rev(c("Meta-analysis (102)","Frontal pool (76)","GSE80655 DLPFC (48)",
                          "GSE78936 BA11 (28)","Hippocampus CA1 (26)","Hippocampus CA3 (26)",
                          "Hippocampus DG (26)"))),
  genes = c(2, 8, 0, 6, 0, 0, 7))
ggplot(sm, aes(genes, analysis)) +
  geom_col(fill="#34495E", width=.6) + geom_text(aes(label=genes), hjust=-.4, size=3.4) +
  xlim(0,10) + theme_bw() +
  labs(x="DE genes (FDR<0.05)", y=NULL,
    title="All analyses on verified unique donors",
    subtitle="No transcriptome-wide signature; the meta-analysed IEG module is the signal")
sv("Fig_6_summary", 8, 5)
logf("\nOutputs -> %s", OUT)
