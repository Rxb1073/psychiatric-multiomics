#!/usr/bin/env Rscript

# =====================================================================

# IS THE IEG EFFECT SEX-SPECIFIC?

#

# The evidence has pointed both ways and needs settling, because the

# source study for GSE102556 (Labonte et al. 2017) is entirely about sex

# specificity in MDD, so this will be asked.

#

# FOR sex-specificity

#   GSE102556 per-gene IEG medians: female -0.827 (16/16 down, 9 nominal)

#                                   male   -0.130 (14/16 down, 0 nominal)

#   CAMERA:  GSE101521 female p=4.95e-03 (n=17) vs male p=0.372 (n=42)

#            -- the female stratum is 40% of the male size with a p-value

#               two orders of magnitude smaller, and the pattern repeats

#               in GSE102556 (female 0.070 vs male 0.167)

#

# AGAINST

#   sex x diagnosis interaction LRT: 0 genes at FDR<0.10

#   p<0.05 excess: female 1.11x, male 0.70x -- neither above chance

#   male/female log2FC correlation +0.107, sign concordance 49.8% (chance)

#   the ERK feedback module flips: female stronger in GSE101521 but MALE

#     stronger in GSE102556 (AUC 0.701 vs 0.577)

#   female strata are tiny: 6 and 9 controls

#

# WHY THE EARLIER TEST WAS INCONCLUSIVE

#   The interaction was tested PER GENE, where power is minimal at these

#   sample sizes. For a set-level effect the right test is the interaction

#   on the SCORE -- one number per donor. That is far better powered and

#   is what this script runs.

#

# ALSO RUN

#   the same test on the ERK-output score, as an independent module

#   a permutation test of the sex difference in effect size

#   the female-only effect at matched sample size, to check whether the

#     apparent female advantage is just small-group variance

#

# READ-ONLY. Writes to MDD/audit/sexspecific/.

# =====================================================================



suppressPackageStartupMessages({ library(DESeq2); library(org.Hs.eg.db); library(AnnotationDbi) })



MDD <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD"

DE  <- file.path(MDD, "transcriptomics/DE")

OUT <- file.path(MDD, "audit/sexspecific"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

hdr <- function(x) cat("\n", strrep("=",76), "\n", x, "\n", strrep("=",76), "\n", sep="")

num <- function(x) suppressWarnings(as.numeric(as.character(x)))

set.seed(42)



IEG <- c("NPAS4","FOS","FOSB","EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3",

         "ATF3","BTG2","ARC","JUN","JUNB","DUSP1","DUSP6")

ERK <- c("DUSP4","DUSP5","SPRY1","SPRY2","SPRY4","SPRED1","SPRED2",

         "ETV4","ETV5","PHLDA1","EPHA2")   # DUSP6 excluded: overlaps IEG



CO <- list(

  GSE101521 = list(cnt=file.path(MDD,"transcriptomics/GSE101521/counts/GSE101521_gene_counts_matrix.tsv"),

                   col=file.path(DE,"GSE101521_colData.tsv")),

  GSE102556 = list(cnt=file.path(MDD,"transcriptomics/GSE102556/counts/GSE102556_gene_counts_matrix.tsv"),

                   col=file.path(DE,"GSE102556_colData.tsv")))



load1 <- function(x) {

  cnt <- read.delim(x$cnt, row.names=1, check.names=FALSE)

  col <- read.delim(x$col, check.names=FALSE); names(col)[1] <- "sample_id"

  rownames(col) <- col$sample_id

  k <- intersect(colnames(cnt), rownames(col))

  cnt <- round(as.matrix(cnt[,k,drop=FALSE])); col <- col[k,,drop=FALSE]

  col$disease <- factor(ifelse(grepl("MDD|depress", col$disease, ignore.case=TRUE),

                               "MDD","Control"), levels=c("Control","MDD"))

  col$sex <- factor(tolower(col$sex))

  keep <- rowSums(cnt >= 10) >= min(table(col$disease))

  dds <- estimateSizeFactors(DESeqDataSetFromMatrix(cnt[keep,], col, design=~disease))

  m <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,

        keys=sub("\\..*","",rownames(cnt)[keep]), keytype="ENSEMBL", columns="SYMBOL"))

  m <- m[!duplicated(m$ENSEMBL),]

  list(nc=counts(dds, normalized=TRUE), col=col,

       sy=setNames(m$SYMBOL, m$ENSEMBL)[sub("\\..*","",rownames(cnt)[keep])]) }

D <- lapply(CO, load1)



score <- function(d, syms) {

  ids <- rownames(d$nc)[which(d$sy %in% syms)]

  if (length(ids) < 3) return(NULL)

  R <- t(apply(d$nc[ids,,drop=FALSE], 1, function(v) rank(v)/length(v)))

  apply(R, 2, median) }

auc <- function(a,b) suppressWarnings(wilcox.test(a,b)$statistic)/(length(a)*length(b))

wp  <- function(a,b) suppressWarnings(wilcox.test(a,b)$p.value)



hdr("0. DESIGN")

S <- list()

for (n in names(D)) {

  cl <- D[[n]]$col

  print(table(cl$sex, cl$disease, dnn=c(n,"")))

  S[[n]] <- data.frame(cohort=n, sample=colnames(D[[n]]$nc),

                       disease=as.character(cl$disease), sex=as.character(cl$sex),

                       ieg=score(D[[n]], IEG), erk=score(D[[n]], ERK),

                       stringsAsFactors=FALSE)

}

A <- do.call(rbind, S)

A$disease <- factor(A$disease, levels=c("Control","MDD"))

A$sex     <- factor(A$sex)

A$cohort  <- factor(A$cohort)

cat(sprintf("\n  total %d donors; IEG score from %d genes, ERK score from %d\n",

            nrow(A), sum(D[[1]]$sy %in% IEG), sum(D[[1]]$sy %in% ERK)))



# ============== 1. THE PROPERLY POWERED INTERACTION TEST ==============

hdr("1. SEX x DIAGNOSIS INTERACTION ON THE SCORE")

cat("  The earlier interaction test was per gene, where power is minimal.\n")

cat("  On the score -- one number per donor -- it is far better powered.\n")

cat("  Scores are computed within cohort, so cohort is included as a term.\n")



for (v in c("ieg","erk")) {

  cat(sprintf("\n  --- %s score ---\n", toupper(v)))

  f0 <- as.formula(paste(v, "~ cohort + sex + disease"))

  f1 <- as.formula(paste(v, "~ cohort + sex * disease"))

  m0 <- lm(f0, data=A); m1 <- lm(f1, data=A)

  an <- anova(m0, m1)

  cat(sprintf("    main effect of diagnosis : p = %.4f\n",

              summary(m0)$coefficients["diseaseMDD","Pr(>|t|)"]))

  cat(sprintf("    sex x diagnosis          : F = %.2f, p = %.4f%s\n",

              an$F[2], an$`Pr(>F)`[2],

              ifelse(an$`Pr(>F)`[2] < 0.05, "   <-- SEX-SPECIFIC", "   (no interaction)")))

  ct <- summary(m1)$coefficients

  ix <- grep(":", rownames(ct))

  if (length(ix)) cat(sprintf("    interaction coefficient  : %+.4f (SE %.4f)\n",

                              ct[ix[1],1], ct[ix[1],2]))

}



# per cohort, in case the interaction is cohort-specific

cat("\n  same test within each cohort separately:\n")

for (n in levels(A$cohort)) {

  d <- A[A$cohort==n,]

  for (v in c("ieg","erk")) {

    m0 <- lm(as.formula(paste(v,"~ sex + disease")), data=d)

    m1 <- lm(as.formula(paste(v,"~ sex * disease")), data=d)

    an <- anova(m0,m1)

    cat(sprintf("    %-10s %s  interaction p = %.4f%s\n", n, toupper(v),

                an$`Pr(>F)`[2], ifelse(an$`Pr(>F)`[2]<0.05,"  *","")))

  }

}



# ============== 2. EFFECT SIZE BY SEX ==============

hdr("2. EFFECT SIZE BY SEX, WITH SAMPLE SIZES")

cat(sprintf("  %-10s %-7s %-7s %5s %5s %8s %9s\n",

            "cohort","score","sex","nCtl","nMDD","AUC","p"))

for (n in levels(A$cohort)) for (v in c("ieg","erk")) for (s in levels(A$sex)) {

  d <- A[A$cohort==n & A$sex==s,]

  a <- d[[v]][d$disease=="Control"]; b <- d[[v]][d$disease=="MDD"]

  if (length(a)<4 || length(b)<4) next

  cat(sprintf("  %-10s %-7s %-7s %5d %5d %8.3f %9.4f%s\n",

              n, toupper(v), s, length(a), length(b), auc(a,b), wp(a,b),

              ifelse(wp(a,b)<0.05,"  *","")))

}

cat("\n  Note the female group sizes. Small groups give more variable AUCs,\n")

cat("  so a higher female AUC is not by itself evidence of sex specificity.\n")



# ============== 3. IS THE FEMALE ADVANTAGE JUST SMALL n? ==============

hdr("3. SUBSAMPLING: DOES A MALE GROUP OF FEMALE SIZE LOOK AS STRONG?")

cat("  Repeatedly draw a male subsample matched to the female group size and\n")

cat("  recompute the AUC. If the observed female AUC sits inside that null,\n")

cat("  the apparent female advantage is explained by sample size alone.\n")

NB <- 2000

for (n in levels(A$cohort)) {

  d <- A[A$cohort==n,]

  fem <- d[d$sex=="female",]; mal <- d[d$sex=="male",]

  nc_ <- sum(fem$disease=="Control"); nm <- sum(fem$disease=="MDD")

  mc <- mal[mal$disease=="Control",]; mm <- mal[mal$disease=="MDD",]

  if (nrow(mc)<nc_ || nrow(mm)<nm) { cat(sprintf("  %-10s male group too small\n", n)); next }

  obs <- auc(fem$ieg[fem$disease=="Control"], fem$ieg[fem$disease=="MDD"])

  null <- replicate(NB, {

    a <- sample(mc$ieg, nc_); b <- sample(mm$ieg, nm); auc(a,b) })

  p <- (sum(null >= obs)+1)/(NB+1)

  cat(sprintf("\n  %-10s female AUC %.3f (n=%d/%d)\n", n, obs, nc_, nm))

  cat(sprintf("             male subsamples of the same size: median %.3f, 95th pct %.3f\n",

              median(null), quantile(null,0.95)))

  cat(sprintf("             p = %.4f%s\n", p,

      ifelse(p<0.05, "   female genuinely stronger",

                     "   female advantage explained by sample size")))

}



# ============== 4. POOLED, SEX AS COVARIATE ==============

hdr("4. THE PRACTICAL QUESTION: DOES SEX CHANGE THE HEADLINE?")

for (v in c("ieg","erk")) {

  m <- lm(as.formula(paste(v,"~ cohort + sex + disease")), data=A)

  ct <- summary(m)$coefficients

  cat(sprintf("\n  %s score, ~ cohort + sex + disease\n", toupper(v)))

  cat(sprintf("    diagnosis %+.4f  p = %.4f\n", ct["diseaseMDD",1], ct["diseaseMDD",4]))

  if ("sexmale" %in% rownames(ct))

    cat(sprintf("    sex       %+.4f  p = %.4f\n", ct["sexmale",1], ct["sexmale",4]))

  a <- A[[v]][A$disease=="Control"]; b <- A[[v]][A$disease=="MDD"]

  cat(sprintf("    unadjusted AUC %.3f, p %.4f\n", auc(a,b), wp(a,b)))

}

write.table(A, file.path(OUT,"scores_by_sex.tsv"), sep="\t", quote=FALSE, row.names=FALSE)



hdr("VERDICT")

cat("

  interaction significant

      the effect is genuinely sex-specific and must be reported that way,

      in line with Labonte. Present male and female results separately.



  interaction NOT significant, and section 3 shows the female advantage is

  explained by sample size

      the effect is present in both sexes. Report it pooled with sex as a

      covariate, and state explicitly that a sex-specific effect was tested

      for and not found -- which is a meaningful negative given the source

      study's emphasis, not an omission.



  Either way, section 4 shows whether the headline changes once sex is in

  the model. If it does not, the sex question affects interpretation but

  not the result.

")

cat("\noutputs in", OUT, "\n")

sessionInfo()

