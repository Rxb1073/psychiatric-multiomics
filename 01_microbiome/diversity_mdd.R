#!/usr/bin/env Rscript
suppressMessages({
  personal <- Sys.getenv("R_LIBS_USER"); if (nchar(personal)==0) personal <- path.expand("~/R/mmuphin_lib")
  .libPaths(c(personal, .libPaths()))
  need <- c("vegan","ggplot2")
  for (p in need) if (!requireNamespace(p, quietly=TRUE))
    install.packages(p, repos="https://cloud.r-project.org", lib=personal)
  library(vegan); library(ggplot2)
})
setwd("/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/merged_genus")

## ---- load RAW counts (for alpha) + metadata ----
raw <- read.delim("merged_genus_raw.tsv", skip=1, header=TRUE, row.names=1,
                  check.names=FALSE, comment.char="")
meta <- read.delim("merged_genus_metadata.tsv", header=TRUE, check.names=FALSE)
rownames(meta) <- meta$`sample-id`
raw <- raw[, rownames(meta)]                       # align order
meta$disease <- factor(meta$disease, levels=c("Control","MDD"))
meta$study   <- factor(meta$study)
counts <- t(as.matrix(raw))                        # samples x genera
cat("counts:", nrow(counts), "samples x", ncol(counts), "genera\n")

## ================= ALPHA DIVERSITY (raw counts) =================
alpha <- data.frame(
  `sample-id` = rownames(counts),
  observed = rowSums(counts > 0),
  shannon  = diversity(counts, index="shannon"),
  simpson  = diversity(counts, index="simpson"),
  disease  = meta$disease, study = meta$study,
  check.names = FALSE
)
write.table(alpha, "alpha_diversity.tsv", sep="\t", quote=FALSE, row.names=FALSE)

## test disease effect, controlling for study (two-way, study as covariate)
cat("\n===== ALPHA: disease effect controlling for study =====\n")
for (m in c("observed","shannon","simpson")) {
  fit <- aov(alpha[[m]] ~ study + disease, data=alpha)
  cat("\n---", m, "---\n"); print(summary(fit)[[1]][c("study","disease"),])
  # also per-study Wilcoxon (nonparametric, robust)
  for (s in levels(alpha$study)) {
    sub <- alpha[alpha$study==s,]
    if (length(unique(sub$disease))==2) {
      w <- wilcox.test(sub[[m]] ~ sub$disease)
      cat(sprintf("   %s within %s: W p = %.4f\n", m, s, w$p.value))
    }
  }
}

## alpha boxplots
p_alpha <- ggplot(alpha, aes(disease, shannon, fill=disease)) +
  geom_boxplot(outlier.size=0.6) + facet_wrap(~study) +
  labs(title="Shannon diversity by disease, per study", y="Shannon", x=NULL) +
  theme_bw() + theme(legend.position="none")
ggsave("alpha_shannon_boxplot.png", p_alpha, width=7, height=4, dpi=150)

## ================= BETA DIVERSITY =================
## (a) corrected table -> Aitchison/Bray ordination (the "after batch removal" view)
corr <- read.delim("merged_genus_corrected.tsv", header=TRUE, row.names=1, check.names=FALSE)
corr <- corr[, rownames(meta)]
corr_t <- t(as.matrix(corr))

## Bray-Curtis on corrected
bc <- vegdist(corr_t, method="bray")
permB <- adonis2(bc ~ study + disease, data=meta, permutations=999, by="margin")
cat("\n===== BETA (corrected, Bray-Curtis) PERMANOVA =====\n"); print(permB)

## PCoA plot
pcoa <- cmdscale(bc, k=2, eig=TRUE)
ev <- round(100*pcoa$eig[1:2]/sum(pcoa$eig[pcoa$eig>0]),1)
pdf_df <- data.frame(PC1=pcoa$points[,1], PC2=pcoa$points[,2],
                     disease=meta$disease, study=meta$study)
p_beta <- ggplot(pdf_df, aes(PC1, PC2, color=disease, shape=study)) +
  geom_point(size=2, alpha=0.8) +
  stat_ellipse(aes(group=disease), type="norm", linetype=2) +
  labs(title="PCoA (Bray-Curtis, batch-corrected)",
       x=paste0("PC1 (",ev[1],"%)"), y=paste0("PC2 (",ev[2],"%)")) +
  theme_bw()
ggsave("beta_pcoa_braycurtis_corrected.png", p_beta, width=7, height=5, dpi=150)

## (b) Aitchison (CLR + Euclidean) on RAW counts as compositional sensitivity
clr <- function(x){ x <- x + 0.5; log(x) - rowMeans(log(x)) }   # pseudocount 0.5
ait <- dist(clr(counts), method="euclidean")
permA <- adonis2(ait ~ study + disease, data=meta, permutations=999, by="margin")
cat("\n===== BETA (raw, Aitchison/CLR) PERMANOVA =====\n"); print(permA)

## dispersion check (PERMANOVA assumes homogeneous dispersion)
bd <- betadisper(bc, meta$disease)
cat("\n===== beta-dispersion (disease) =====\n"); print(anova(bd))
cat("(If dispersion p<0.05, PERMANOVA differences may be partly dispersion, not location — note it.)\n")

saveRDS(list(alpha=alpha, permB=permB, permA=permA, disp=anova(bd)), "diversity_result.rds")
cat("\nWROTE: alpha_diversity.tsv, *.png, diversity_result.rds\n")
