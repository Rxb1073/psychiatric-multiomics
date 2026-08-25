## Cross-omics KEGG axis: gut microbial pathway shift vs host cortical pathway shift.
##
## DESIGN. No donor is measured in both layers, so this is diagonal (results-level)
## integration: the two layers meet only in the KEGG pathway identifier space,
## where a microbial map koNNNNN and a human map hsaNNNNN share a five-digit suffix.
## The test is threshold-free by design -- Spearman correlation across ALL pathways
## present in both layers, with no significance cutoff on either side. This uses the
## full distribution rather than an intersection of two thresholded lists, matters
## because MDD has no KEGG pathway at q<0.05 after class filtering, and limits the
## influence of any single pathway driven by a handful of shared KOs.
##
## WHAT THE IDENTIFIERS DO AND DO NOT MEAN. ko00020 and hsa00020 are the same KEGG
## map annotated in two kingdoms. Correspondence of identifiers is not correspondence
## of biology: the bacterial and human enzyme sets differ, and a correlation here
## means the two layers shift in a consistent direction across shared metabolic
## categories, not that the same molecular process is being measured twice.
##
## Microbial side : MaAsLin2 coefficient (log abundance, case vs control)
## Host side      : GSEA normalised enrichment score (case vs control)
## Significance   : permutation on the pairing, which respects the observed
##                  marginal distributions and makes no normality assumption
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
set.seed(1)

INT <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final/INTEGRATION_inputs"
AX  <- file.path(INT, "06_kegg_axis")
OUT <- file.path(AX, "cross_omics"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
NPERM <- 10000

## ASD has no transcriptomic KEGG object: ranked-list enrichment returned zero
## significant terms in KEGG and GO, so there is no host pathway layer to pair.
## The autism arm enters at gene level (CHD8 targets) instead -- see constraint 4.
ARMS <- c("SCZ","MDD")

report <- list(); pooled <- list()

for(a in ARMS){
  cat(sprintf("\n================ %s ================\n", a))

  mb <- read.delim(file.path(AX, sprintf("KO_pathway_DA_%s_filtered.tsv", a)), check.names=FALSE)
  mb$suffix <- sub("^ko", "", mb$feature)
  cat(sprintf("  microbial: %d prokaryote-relevant pathways\n", nrow(mb)))

  hs <- read.csv(file.path(INT, "02_pathway_KEGG", sprintf("%s_KEGG.csv", a)), check.names=FALSE)
  hs$suffix <- sub("^hsa", "", hs$ID)
  cat(sprintf("  host:      %d KEGG pathways (GSEA)\n", nrow(hs)))

  ## length-adjustment survivors only, where that column is present: transcript
  ## length correlates with power and with functional category
  if("survives_length_adjustment" %in% names(hs)){
    keep <- hs$survives_length_adjustment %in% c(TRUE,"TRUE","True","yes")
    cat(sprintf("  host:      %d of %d survive length adjustment\n", sum(keep), nrow(hs)))
    hs <- hs[keep, ]
  }

  m <- merge(mb[,c("suffix","feature","pathway_name","sub","coef","pval","qval_refiltered")],
             hs[,c("suffix","ID","Description","NES","pvalue","p.adjust")],
             by="suffix")
  cat(sprintf("  SHARED:    %d pathways present in both layers\n", nrow(m)))
  if(nrow(m) < 10){ cat("  too few shared pathways to correlate\n"); next }

  rho  <- cor(m$coef, m$NES, method="spearman")
  ## permutation: shuffle the pairing, keep both marginals
  null <- replicate(NPERM, cor(m$coef, sample(m$NES), method="spearman"))
  pperm <- (sum(abs(null) >= abs(rho)) + 1) / (NPERM + 1)
  ct <- cor.test(m$coef, m$NES, method="spearman", exact=FALSE)

  cat(sprintf("\n  Spearman rho = %+.3f   permutation p = %.4f   (asymptotic p = %.4f)\n",
              rho, pperm, ct$p.value))
  cat(sprintf("  direction agreement: %d of %d pathways (%.1f%%) share sign\n",
              sum(sign(m$coef)==sign(m$NES)), nrow(m),
              100*mean(sign(m$coef)==sign(m$NES))))

  ## concordant pathways where BOTH layers are at least nominally moving,
  ## listed for interpretation only -- not a second test
  m$concordant <- sign(m$coef)==sign(m$NES) & m$pval < 0.05 & m$pvalue < 0.05
  cat(sprintf("\n  pathways nominally moving in BOTH layers, same direction: %d\n", sum(m$concordant)))
  if(any(m$concordant)){
    cc <- m[m$concordant, c("feature","ID","pathway_name","sub","coef","NES","qval_refiltered","p.adjust")]
    print(head(cc[order(-abs(cc$NES)), ], 15), row.names=FALSE)
  }

  m$arm <- a
  write.csv(m[order(m$qval_refiltered), ], file.path(OUT, sprintf("kegg_axis_%s.csv", a)), row.names=FALSE)
  pooled[[a]] <- m
  report[[a]] <- data.frame(arm=a, shared=nrow(m), rho=rho, p_perm=pperm,
                            p_asymptotic=ct$p.value,
                            pct_same_sign=100*mean(sign(m$coef)==sign(m$NES)),
                            both_nominal=sum(m$concordant))

  pdf(file.path(OUT, sprintf("kegg_axis_%s.pdf", a)), width=5.5, height=5.5)
  plot(m$coef, m$NES, pch=21, bg=ifelse(m$concordant,"grey30","white"), cex=0.9,
       xlab=sprintf("microbial pathway coefficient (%s vs control)", a),
       ylab="host cortical NES", main=sprintf("%s: KEGG axis, rho = %+.3f, p = %.3g", a, rho, pperm))
  abline(h=0, v=0, col="grey70"); abline(lm(NES ~ coef, m), col="firebrick", lwd=2)
  dev.off()
}

cat("\n================ summary ================\n")
if(length(report)){
  S <- do.call(rbind, report); print(S, row.names=FALSE)
  write.csv(S, file.path(OUT,"kegg_axis_summary.csv"), row.names=FALSE)
}
cat(sprintf("\nwritten to %s\n", OUT))
cat("\nNOTE: ASD is absent by design -- no transcriptomic KEGG object exists for that\n")
cat("arm (zero significant terms from the ranked list). Autism's cross-omics link is\n")
cat("at gene level via CHD8 binding targets.\n")
