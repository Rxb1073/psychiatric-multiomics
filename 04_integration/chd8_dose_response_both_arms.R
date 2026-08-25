## Is the CHD8 dose-response autism-specific?
##
## CHD8-bound promoters are displaced upward in the autism transcriptome
## (p < 0.0001 against expression-matched nulls) AND, more weakly, in the
## schizophrenia transcriptome (p = 0.032). A set that moves in both disorders
## is not by itself an autism finding: it may reflect the general transcriptional
## state of CHD8-bound promoters rather than anything disorder-specific.
##
## The autism result has one feature the cross-disorder test did not examine: the
## displacement scales monotonically with binding strength, rising from +0.019 to
## +0.028 log2 units as the set is restricted from 4,351 sites to the strongest
## 1,000. A gradient aligned with immunoprecipitation signal is difficult to
## attribute to a general property of active promoters.
##
## This script asks whether that gradient also appears in schizophrenia. If it
## does, the autism finding is a stronger instance of a general phenomenon. If it
## does not, the dose-response remains autism-specific and the finding survives
## in a narrower but defensible form.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
set.seed(1)

D   <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final"
INT <- file.path(D,"INTEGRATION_inputs")
OUT <- file.path(INT,"13_epigenomic_cross_disorder"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
NPERM <- 10000

p2g <- read.csv(file.path(INT,"07_chd8_axis","CHD8_peak_to_gene.csv"), stringsAsFactors=FALSE)
p2g <- p2g[order(-p2g$score), ]
cat(sprintf("CHD8 peak-to-gene: %d rows | score range %.2f-%.2f\n",
            nrow(p2g), min(p2g$score, na.rm=TRUE), max(p2g$score, na.rm=TRUE)))

load_arm <- function(){
  asd <- read.csv(file.path(INT,"05_genelevel","ASD_meta.csv"), check.names=FALSE)
  asd <- data.frame(key=sub("\\..*$","",asd$ens), beta=asd$beta_FE,
                    abund=rowMeans(asd[,c("baseMean_1","baseMean_2")], na.rm=TRUE),
                    stringsAsFactors=FALSE)
  scz <- read.csv(file.path(INT,"05_genelevel","SCZ_frontal.csv"), check.names=FALSE)
  names(scz)[1] <- "ens"
  scz <- data.frame(key=sub("\\..*$","",scz$ens), beta=scz$logFC,
                    abund=scz$AveExpr, stringsAsFactors=FALSE)
  lapply(list(ASD=asd, SCZ=scz), function(x){
    x <- x[!is.na(x$beta) & !is.na(x$abund) & !duplicated(x$key), ]
    x$stratum <- cut(rank(x$abund, ties.method="first"),
                     breaks=quantile(rank(x$abund, ties.method="first"), probs=seq(0,1,0.05)),
                     include.lowest=TRUE, labels=FALSE)
    x })
}
arms <- load_arm()
for(a in names(arms)) cat(sprintf("%s transcriptome: %d genes\n", a, nrow(arms[[a]])))

matched <- function(x, genes){
  idx <- which(x$key %in% genes)
  if(length(idx) < 20) return(NULL)
  mt   <- median(x$beta[idx])
  tab  <- table(x$stratum[idx]); rest <- setdiff(seq_len(nrow(x)), idx)
  pool <- split(rest, x$stratum[rest])
  null <- replicate(NPERM, {
    s <- unlist(lapply(names(tab), function(k)
      if(!is.null(pool[[k]])) sample(pool[[k]], min(tab[[k]], length(pool[[k]]))) else integer(0)))
    median(x$beta[s]) })
  data.frame(n=length(idx), median=mt, null=median(null), sd_null=sd(null),
             delta=mt - median(null),
             p=(sum(abs(null-median(null)) >= abs(mt-median(null))) + 1)/(NPERM+1))
}

SIZES <- c(1000, 2000, 3000, nrow(p2g))
res <- list()
for(a in names(arms)){
  cat(sprintf("\n================ %s ================\n", a))
  for(N in SIZES){
    g <- head(p2g$ens, N)
    r <- matched(arms[[a]], g)
    if(is.null(r)){ cat(sprintf("  top %d: too few genes\n", N)); next }
    r <- cbind(arm=a, sites=N, r)
    cat(sprintf("  top %5d sites | %4d genes | median %+.4f | null %+.4f | delta %+.4f | p = %.4f\n",
                N, r$n, r$median, r$null, r$delta, r$p))
    res[[paste(a,N)]] <- r
  }
}
S <- do.call(rbind, res); rownames(S) <- NULL

## ---- is the gradient present in each arm? ----------------------------------
cat("\n", strrep("=",72), "\nDOSE-RESPONSE WITH BINDING STRENGTH\n", strrep("=",72), "\n", sep="")
for(a in names(arms)){
  s <- S[S$arm==a, ]
  s <- s[order(s$sites), ]
  rho <- cor(s$sites, s$delta, method="spearman")
  cat(sprintf("\n%s: delta from %+.4f (all sites) to %+.4f (top 1,000)\n",
              a, s$delta[nrow(s)], s$delta[1]))
  cat(sprintf("  rank correlation between set size and effect: %+.2f  (%s)\n", rho,
              if(rho <= -0.8) "MONOTONIC -- effect grows as binding evidence strengthens"
              else "no consistent gradient"))
  cat(sprintf("  fold change in effect, all sites -> top 1,000: %.2fx\n",
              s$delta[1]/s$delta[nrow(s)]))
}

cat("\n", strrep("=",72), "\nSUMMARY\n", strrep("=",72), "\n", sep="")
print(S[,c("arm","sites","n","median","null","delta","p")], row.names=FALSE)
write.csv(S, file.path(OUT,"CHD8_dose_response_both_arms.csv"), row.names=FALSE)
cat(sprintf("\nwritten to %s\n", OUT))
