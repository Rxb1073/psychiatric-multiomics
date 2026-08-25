#!/usr/bin/env Rscript
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(DESeq2)})
logf <- function(...) cat(sprintf(...), "\n")

## ---- locate project roots (they differ between arms) ----
roots <- c("/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq",
           "/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/RNA_seq")
findf <- function(rel) { for (r in roots) if (file.exists(file.path(r,rel))) return(file.path(r,rel))
                         stop("not found: ", rel) }
CS  <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/counts_scz"
SH  <- file.path(CS,"sheets")
OUT <- file.path(CS,"PERCOHORT_prua"); dir.create(OUT, showWarnings=FALSE)

C80 <- findf("PRJNA319583/counts/PRJNA319583_gene_counts_matrix.tsv")
C78 <- findf("PRJNA314463/counts/PRJNA314463_gene_counts.txt")
BTF <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/gene_id2biotype.tsv"
logf("counts GSE80655: %s", C80); logf("counts GSE78936: %s", C78)

read_counts <- function(p){
  cm <- read.delim(p, comment.char="#", check.names=FALSE)
  g <- cm[[1]]; cm[[1]] <- NULL
  drop <- intersect(c("Chr","Start","End","Strand","Length"), colnames(cm))
  if (length(drop)) cm <- cm[, setdiff(colnames(cm), drop), drop=FALSE]
  colnames(cm) <- regmatches(colnames(cm), regexpr("[SED]RR[0-9]+", colnames(cm)))
  m <- as.matrix(cm); rownames(m) <- g; storage.mode(m) <- "integer"; m
}
bt  <- read.delim(BTF, header=FALSE, col.names=c("gid","bt"))
keep_bt <- bt$gid[bt$bt %in% c("protein_coding","lncRNA")]

## ---- assemble one cohort: counts + sheet + PRUA ----
build <- function(cfile, sname){
  ss <- read.delim(file.path(SH, paste0(sname,"_samplesheet.tsv")), check.names=FALSE)
  pr <- read.delim(file.path(SH, paste0(sname,"_prua.tsv")),        check.names=FALSE)
  ss <- merge(ss, pr[,c("srr","prua")], by="srr")
  rownames(ss) <- ss$srr
  cm <- read_counts(cfile)
  k  <- intersect(colnames(cm), rownames(ss))
  cm <- cm[, k, drop=FALSE]; ss <- ss[k, , drop=FALSE]

  ## --- collapse technical run-splits: sum counts per biological sample ---
  RT <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/GSE80655/metadata/PRJNA319583_SraRunTable.csv"
  if (file.exists(RT)) {
    rt  <- read.csv(RT, check.names=FALSE)
    map <- setNames(rt$BioSample, rt$Run)
    bio <- map[colnames(cm)]
    if (!any(is.na(bio)) && length(unique(bio)) < ncol(cm)) {
      logf("    run-splits: %d runs -> %d biological samples (summing counts)",
           ncol(cm), length(unique(bio)))
      cm2 <- vapply(split(seq_len(ncol(cm)), bio),
                    function(i) as.integer(rowSums(cm[, i, drop=FALSE])),
                    integer(nrow(cm)))
      rownames(cm2) <- rownames(cm)
      keep <- !duplicated(bio); ss <- ss[keep, , drop=FALSE]
      rownames(ss) <- bio[keep]
      cm <- cm2[, rownames(ss), drop=FALSE]
      storage.mode(cm) <- "integer"
    }
  }
  cm <- cm[rownames(cm) %in% keep_bt, , drop=FALSE]
  ss$diagnosis <- relevel(factor(ss$diagnosis), ref="Control")
  ss$prua_z <- as.numeric(scale(ss$prua))
  avail <- character(0)
  for (v in c("age","pmi","ph")) {
    if (!(v %in% names(ss))) next
    x <- suppressWarnings(as.numeric(as.character(ss[[v]])))
    frac <- sum(!is.na(x))/length(x)
    if (frac >= 0.9) {
      x[is.na(x)] <- median(x, na.rm=TRUE)
      ss[[paste0(v,"_z")]] <- as.numeric(scale(x))
      avail <- c(avail, v)
      logf("    covariate %-4s: %.0f%% complete -> USED", v, 100*frac)
    } else logf("    covariate %-4s: %.0f%% complete -> DROPPED", v, 100*frac)
  }
  attr(ss,"avail") <- avail
  ng <- min(table(ss$diagnosis))
  cm <- cm[rowSums(cm >= 10) >= ng, , drop=FALSE]
  logf("%s: n=%d (%s) | %d genes after biotype+lowcount",
       sname, nrow(ss), paste(names(table(ss$diagnosis)), table(ss$diagnosis),
       sep="=", collapse=" "), nrow(cm))
  list(cm=cm, ss=ss)
}

fit <- function(d, form, lab, tag){
  dds <- suppressMessages(DESeq(DESeqDataSetFromMatrix(d$cm, d$ss, form), quiet=TRUE))
  r   <- as.data.frame(results(dds, name=tail(resultsNames(dds),1)))   # MLE LFC, unshrunken
  r$gene <- rownames(r)
  s <- function(l) sum(!is.na(r$padj) & r$padj < 0.05 & abs(r$log2FoldChange) >= l)
  logf("  %-38s padj05: %5d | >=0.3: %5d | >=0.5: %5d | >=1.0: %5d",
       lab, s(0), s(0.3), s(0.5), s(1.0))
  write.csv(r, file.path(OUT, paste0(tag,"_full.csv")), row.names=FALSE)
  r
}

## ================= GSE80655 (full covariates) =================
logf("\n=========== GSE80655 DLPFC ===========")
d80 <- build(C80, "GSE80655_DLPFC")
fit(d80, ~ diagnosis,                                    "D0 diagnosis only",        "GSE80655_D0")
av <- attr(d80$ss,"avail"); cz <- paste0(av,"_z")
logf("  covariates used: %s", paste(av, collapse=", "))
fit(d80, as.formula(paste("~", paste(c(cz,"diagnosis"), collapse=" + "))),
    paste0("D1 +", paste(av, collapse="+")), "GSE80655_D1")
r80 <- fit(d80, as.formula(paste("~", paste(c(cz,"prua_z","diagnosis"), collapse=" + "))),
    "D2 +PRUA  [Ramaker design]", "GSE80655_D2")
fit(d80, as.formula(paste("~", paste(c("sex",cz,"prua_z","diagnosis"), collapse=" + "))),
    "D3 +sex", "GSE80655_D3")

## ================= GSE78936 (covariate-poor) =================
logf("\n=========== GSE78936 frontal ===========")
d78 <- build(C78, "GSE78936")
d78$ss$brain_region <- factor(d78$ss$brain_region)
fit(d78, ~ diagnosis,                                    "E0 diagnosis only",        "GSE78936_E0")
fit(d78, ~ brain_region + diagnosis,                     "E1 +region",               "GSE78936_E1")
r78 <- fit(d78, ~ brain_region + prua_z + diagnosis,     "E2 +region+PRUA",          "GSE78936_E2")

## ================= reproducible core =================
logf("\n=========== REPRODUCIBLE CORE (D2 vs E2) ===========")
for (L in c(0, 0.3, 0.5, 1.0)) {
  a <- r80$gene[!is.na(r80$padj) & r80$padj<0.05 & abs(r80$log2FoldChange)>=L]
  b <- r78$gene[!is.na(r78$padj) & r78$padj<0.05 & abs(r78$log2FoldChange)>=L]
  i <- intersect(a,b)
  cc <- if (length(i)) sum(sign(r80$log2FoldChange[match(i,r80$gene)]) ==
                           sign(r78$log2FoldChange[match(i,r78$gene)])) else 0
  logf("  |LFC|>=%.1f : GSE80655=%d  GSE78936=%d  INTERSECT=%d  same-direction=%d/%d",
       L, length(a), length(b), length(i), cc, length(i))
  if (L==0.3 && length(i)) write.csv(r80[match(i,r80$gene),],
      file.path(OUT,"CORE_intersect_LFC0.3.csv"), row.names=FALSE)
}
m <- merge(r80[,c("gene","log2FoldChange")], r78[,c("gene","log2FoldChange")], by="gene")
logf("  genome-wide LFC correlation between cohorts: r = %.3f",
     cor(m[,2], m[,3], use="complete.obs"))
logf("\nOutputs -> %s", OUT)

## ================= GSE78936, BA11 ONLY (unique donors, one region) =================
logf("\n=========== GSE78936 BA11-only ===========")
ba <- d78$ss$brain_region == "Brodmann's Area 11"
ssB <- d78$ss[ba, , drop=FALSE]
cmB <- d78$cm[, rownames(ssB), drop=FALSE]
ssB$diagnosis <- droplevels(factor(ssB$diagnosis))
ssB$prua_z <- as.numeric(scale(ssB$prua))          # rescale within this subset
ngB <- min(table(ssB$diagnosis))
cmB <- cmB[rowSums(cmB >= 10) >= ngB, , drop=FALSE]
logf("BA11: n=%d (%s) | %d genes after refilter", nrow(ssB),
     paste(names(table(ssB$diagnosis)), table(ssB$diagnosis), sep="=", collapse=" "),
     nrow(cmB))
dB <- list(cm=cmB, ss=ssB)
fit(dB, ~ diagnosis,           "F0 diagnosis only",  "GSE78936_BA11_F0")
rB <- fit(dB, ~ prua_z + diagnosis, "F1 +PRUA",      "GSE78936_BA11_F1")

## ---- reproducible core, corrected: GSE80655 D2 vs GSE78936 BA11 F1 ----
logf("\n=========== REPRODUCIBLE CORE (D2 vs BA11 F1) ===========")
for (L in c(0, 0.3, 0.5, 1.0)) {
  a <- r80$gene[!is.na(r80$padj) & r80$padj<0.05 & abs(r80$log2FoldChange)>=L]
  b <- rB$gene [!is.na(rB$padj)  & rB$padj <0.05 & abs(rB$log2FoldChange) >=L]
  i <- intersect(a,b)
  cc <- if (length(i)) sum(sign(r80$log2FoldChange[match(i,r80$gene)]) ==
                           sign(rB$log2FoldChange [match(i,rB$gene)])) else 0
  logf("  |LFC|>=%.1f : GSE80655=%d  BA11=%d  INTERSECT=%d  same-direction=%d/%d",
       L, length(a), length(b), length(i), cc, length(i))
  if (L==0.3 && length(i))
    write.csv(r80[match(i,r80$gene),],
              file.path(OUT,"CORE_intersect_BA11_LFC0.3.csv"), row.names=FALSE)
}
mB <- merge(r80[,c("gene","log2FoldChange")], rB[,c("gene","log2FoldChange")], by="gene")
logf("  genome-wide LFC correlation (GSE80655 vs BA11): r = %.3f",
     cor(mB[,2], mB[,3], use="complete.obs"))

## ---- did dropping the pseudoreplicated donors change GSE78936? ----
logf("\n--- GSE78936: all-52 (E2) vs BA11-only (F1) ---")
mE <- merge(r78[,c("gene","log2FoldChange","padj")], rB[,c("gene","log2FoldChange","padj")],
            by="gene", suffixes=c(".all",".ba11"))
sA <- mE$gene[!is.na(mE$padj.all)  & mE$padj.all <0.05]
sB <- mE$gene[!is.na(mE$padj.ba11) & mE$padj.ba11<0.05]
logf("  sig all-52=%d | sig BA11=%d | overlap=%d (%.0f%% of all-52)",
     length(sA), length(sB), length(intersect(sA,sB)),
     100*length(intersect(sA,sB))/max(1,length(sA)))
logf("  LFC correlation all-52 vs BA11: r = %.3f",
     cor(mE$log2FoldChange.all, mE$log2FoldChange.ba11, use="complete.obs"))
