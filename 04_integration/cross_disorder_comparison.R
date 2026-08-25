## Cross-disorder comparison of the three transcriptomic arms.
##
## DESIGN. Unlike the multi-omic integration, this comparison is between three
## analyses of the SAME data type, quantified against one shared annotation.
## Gene identifiers are byte-identical across arms (one Ensembl GRCh38.109 GTF,
## referenced by the same absolute path from every count job), so genes join
## exactly rather than by symbol matching.
##
## INDEPENDENCE. The schizophrenia and depression primary analyses share no
## samples: intersecting their run accessions returns zero (107 and 109
## accessions). GSE80655 contains both diagnoses and contributes to the
## schizophrenia frontal analysis, but was excluded from the depression primary,
## where it correlated with neither other cohort. Correlations between arms are
## therefore between independent sample sets.
##
## EFFECT SIZES. Unshrunken coefficients are used throughout. Shrinkage
## estimators differ between arms (ashr in autism, apeglm or normal in
## schizophrenia), so shrunken values are not comparable across arms; maximum
## likelihood coefficients are.
##
## THE CHECK THAT MATTERS. All three arms have effect size distributions skewed
## in one direction. A correlation between two such distributions can arise from
## their shapes rather than from gene-level correspondence -- the same artefact
## that produced spurious gut-brain correlations in the integration chapter.
## Every correlation is therefore recomputed within sign strata and with extreme
## values excluded. An association that reflects real correspondence persists;
## one that reflects distribution shape does not.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(limma)})
set.seed(1)

D   <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final"
INT <- file.path(D,"INTEGRATION_inputs")
OUT <- file.path(INT,"10_cross_disorder"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
NPERM <- 10000

## ---- load the three arms on a common gene key ------------------------------
asd <- read.csv(file.path(INT,"05_genelevel","ASD_meta.csv"), check.names=FALSE)
asd <- data.frame(ens=sub("\\..*$","",asd$ens), sym=asd$sym,
                  beta=asd$beta_FE, p=asd$p_FE, stringsAsFactors=FALSE)

scz <- read.csv(file.path(INT,"05_genelevel","SCZ_frontal.csv"), check.names=FALSE)
names(scz)[1] <- "ens"
scz <- data.frame(ens=sub("\\..*$","",scz$ens), sym=scz$symbol,
                  beta=scz$logFC, p=scz$P.Value, stringsAsFactors=FALSE)

## MDD: the pooled two-cohort DESeq2 result, which carries RAW p-values.
## The symlinked MDD_meta.tsv holds only adjusted values and cannot be ranked.
mdd <- read.delim(file.path(D,"Transcriptomics_MDD/tables/01_PRIMARY_pooled_DE.tsv"),
                  check.names=FALSE)
mdd <- data.frame(ens=sub("\\..*$","",mdd$gene), sym=mdd$symbol,
                  beta=mdd$log2FoldChange, p=mdd$pvalue, stringsAsFactors=FALSE)

arms <- list(ASD=asd, SCZ=scz, MDD=mdd)
for(a in names(arms)){
  x <- arms[[a]]
  x <- x[!is.na(x$beta) & !is.na(x$p) & !duplicated(x$ens), ]
  x$p[x$p <= 0] <- .Machine$double.xmin
  x$stat <- sign(x$beta) * -log10(x$p)
  arms[[a]] <- x
  cat(sprintf("%s: %d genes | median beta %+.4f | %.1f%% negative\n",
              a, nrow(x), median(x$beta), 100*mean(x$beta < 0)))
}

## ---- pairwise concordance, with the stratified check -----------------------
pairwise <- function(a, b){
  m <- merge(arms[[a]][,c("ens","sym","beta","stat")],
             arms[[b]][,c("ens","beta","stat")], by="ens", suffixes=c("_1","_2"))
  rho  <- cor(m$beta_1, m$beta_2, method="spearman")
  null <- replicate(NPERM, cor(m$beta_1, sample(m$beta_2), method="spearman"))
  pp   <- (sum(abs(null) >= abs(rho)) + 1)/(NPERM + 1)

  ## does the association hold WITHIN a sign stratum, or only across the boundary?
  d <- m[m$beta_1 < 0, ]; u <- m[m$beta_1 > 0, ]
  rho_d <- cor(d$beta_1, d$beta_2, method="spearman")
  rho_u <- cor(u$beta_1, u$beta_2, method="spearman")
  ## and without the tails, which carry disproportionate leverage
  q <- quantile(abs(m$beta_1), 0.95)
  k <- abs(m$beta_1) <= q & abs(m$beta_2) <= quantile(abs(m$beta_2), 0.95)
  rho_trim <- cor(m$beta_1[k], m$beta_2[k], method="spearman")

  cat(sprintf("\n--- %s vs %s: %d shared genes ---\n", a, b, nrow(m)))
  cat(sprintf("  rho (all genes)          = %+.3f   permutation p = %.4f\n", rho, pp))
  cat(sprintf("  rho within %s-down       = %+.3f  (n = %d)\n", a, rho_d, nrow(d)))
  cat(sprintf("  rho within %s-up         = %+.3f  (n = %d)\n", a, rho_u, nrow(u)))
  cat(sprintf("  rho excluding top 5%% tails= %+.3f\n", rho_trim))
  cat(sprintf("  same direction: %.1f%% of genes\n", 100*mean(sign(m$beta_1)==sign(m$beta_2))))
  verdict <- if(sign(rho_d)==sign(rho) && sign(rho_u)==sign(rho) &&
                min(abs(c(rho_d,rho_u))) > 0.5*abs(rho))
             "HOLDS within strata" else "DOES NOT HOLD within strata -- likely distributional"
  cat(sprintf("  => %s\n", verdict))
  write.csv(m[order(-abs(m$stat_1 + m$stat_2)), ],
            file.path(OUT, sprintf("genes_%s_vs_%s.csv", a, b)), row.names=FALSE)
  data.frame(pair=sprintf("%s vs %s", a, b), shared=nrow(m), rho=rho, p_perm=pp,
             rho_down=rho_d, rho_up=rho_u, rho_trimmed=rho_trim,
             pct_same_sign=100*mean(sign(m$beta_1)==sign(m$beta_2)), verdict=verdict)
}

cat("\n", strrep("=",70), "\nPAIRWISE EFFECT-SIZE CONCORDANCE\n", strrep("=",70), "\n", sep="")
P <- rbind(pairwise("ASD","SCZ"), pairwise("ASD","MDD"), pairwise("SCZ","MDD"))

## ---- set-level testing, identical sets across all three arms ---------------
## The transcriptomics chapter tested these sets in SCZ and MDD but not ASD.
## Running the same sets in all three makes the set-level comparison like-for-like.
SETS <- list(
  GABAergic      = c("SST","PVALB","GAD1","GAD2","VIP","NPY","CCK","LHX6","RELN","CNR1"),
  IEG            = c("EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3","FOS","FOSB","JUN","JUNB",
                     "ARC","ATF3","DUSP1","DUSP6","NPAS4"),
  excitatory     = c("SLC17A7","SLC17A6","CAMK2A","SATB2","TBR1","NEUROD6"),
  astrocyte      = c("GFAP","ALDH1L1","AQP4","S100B","SLC1A2","SLC1A3","VIM"),
  oligodendrocyte= c("MBP","PLP1","MOG","MAG","CNP","SOX10"),
  CalbindinCR    = c("CALB1","CALB2"))

cat("\n", strrep("=",70), "\nSET-LEVEL TESTING, SAME SETS IN ALL THREE ARMS\n", strrep("=",70), "\n", sep="")
cat("cameraPR on the ranked list; inter-gene correlation must be ASSUMED because\n")
cat("the input is a summary table. Reported at 0.01 and 0.05; the transcriptomics\n")
cat("chapter estimated 0.34-0.80 from data for these sets, so 0.01 is not defensible.\n")
S <- list()
for(a in names(arms)){
  x <- arms[[a]]
  st <- setNames(x$stat, x$sym)
  st <- st[!is.na(names(st)) & !duplicated(names(st))]
  idx <- lapply(SETS, function(s) which(names(st) %in% s))
  idx <- idx[sapply(idx, length) >= 2]
  for(r in c(0.01, 0.05)){
    cam <- cameraPR(st, idx, inter.gene.cor=r)
    S[[paste(a,r)]] <- data.frame(arm=a, inter_gene_cor=r, set=rownames(cam),
                                  NGenes=cam$NGenes, Direction=cam$Direction,
                                  PValue=cam$PValue, FDR=cam$FDR)
  }
}
SS <- do.call(rbind, S); rownames(SS) <- NULL
for(r in c(0.01,0.05)){
  cat(sprintf("\n--- assumed inter-gene correlation %.2f ---\n", r))
  w <- reshape(SS[SS$inter_gene_cor==r, c("arm","set","Direction","FDR")],
               idvar="set", timevar="arm", direction="wide")
  print(w, row.names=FALSE)
}

write.csv(P,  file.path(OUT,"cross_disorder_concordance.csv"), row.names=FALSE)
write.csv(SS, file.path(OUT,"cross_disorder_setlevel.csv"), row.names=FALSE)
cat("\n", strrep("=",70), "\nSUMMARY\n", strrep("=",70), "\n", sep="")
print(P[,c("pair","shared","rho","p_perm","rho_down","rho_up","pct_same_sign","verdict")], row.names=FALSE)
cat(sprintf("\nwritten to %s\n", OUT))
