#!/usr/bin/env Rscript
suppressMessages({ library(DESeq2); library(sva) })
setwd("/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/DE")
ANN <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis"
GTF <- "/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf"

col <- read.delim("rna_colData.tsv", check.names=FALSE)
rownames(col) <- col$`sample-id`
cnt <- read.delim("merged_counts.tsv", row.names=1, check.names=FALSE)
cnt <- round(as.matrix(cnt[, rownames(col)]))
stopifnot(all(colnames(cnt)==rownames(col)))

bt  <- read.delim(file.path(ANN,"gene_id2biotype.tsv"), header=FALSE, col.names=c("gid","bt"))
btv <- setNames(bt$bt, bt$gid)
len <- read.delim(file.path(ANN,"gene_lengths.tsv"), header=FALSE, col.names=c("gid","len"))
lenv<- setNames(as.numeric(len$len), len$gid)
g_bt <- rownames(cnt)[!is.na(btv[rownames(cnt)]) & btv[rownames(cnt)] %in% c("protein_coding","lncRNA")]
cmq  <- cnt[g_bt, , drop=FALSE]
lnc  <- g_bt[btv[g_bt]=="lncRNA"]; pc <- g_bt[btv[g_bt]=="protein_coding"]
col$q_lnc <- colSums(cmq[lnc,,drop=FALSE]) / colSums(cmq)
Lp <- lenv[pc]; ok <- !is.na(Lp) & Lp>0
W  <- cmq[pc,,drop=FALSE][ok,,drop=FALSE]; Lo <- Lp[ok]
col$q_len <- apply(W, 2, function(x){ s<-sum(x); if(is.na(s)||s==0) NA_real_ else sum(x*Lo)/s })

col$disease <- factor(col$disease, levels=c("Control","MDD"))
col$study   <- factor(col$study)
col$sex     <- factor(tolower(substr(col$sex,1,1)))
col$pmi     <- suppressWarnings(as.numeric(col$pmi))
col$ph      <- suppressWarnings(as.numeric(col$ph))
for (s in levels(col$study)) {
  idx <- col$study==s & is.na(col$pmi)
  if (any(idx)) col$pmi[idx] <- median(col$pmi[col$study==s], na.rm=TRUE)
}
for (cv in c("pmi","ph","q_len","q_lnc")) col[[paste0(cv,"_z")]] <- as.numeric(scale(col[[cv]]))

cat("samples:", nrow(col), "\n"); print(table(col$study, col$disease))
keep <- rowSums(cnt >= 10) >= 0.5*ncol(cnt)
cnt <- cnt[keep, ]
cat("genes after filter:", nrow(cnt), "\n")

sym <- {
  s <- c(); con <- file(GTF,"r")
  repeat { L <- readLines(con, n=50000); if(length(L)==0) break
    L <- L[grepl("\tgene\t", L)]
    gid <- sub('.*gene_id "([^"]+)".*','\\1',L); nm <- ifelse(grepl("gene_name",L), sub('.*gene_name "([^"]+)".*','\\1',L), gid)
    s[gid] <- nm }
  close(con); s
}
IEG <- c("NPAS4","FOS","FOSB","EGR1","EGR2","EGR3","EGR4","NR4A1","NR4A2","NR4A3",
         "ARC","ATF3","BTG2","DUSP1","JUNB","CYR61","CCN1","HOMER1")

runDE <- function(counts, meta, label, quality=TRUE) {
  form <- if (quality) as.formula("~ sex + pmi_z + ph_z + q_len_z + q_lnc_z + disease") else as.formula("~ sex + pmi_z + ph_z + disease")
  adj <- ComBat_seq(counts, batch=meta$study, group=meta$disease)
  dds <- DESeqDataSetFromMatrix(adj, meta, design=form)
  dds <- DESeq(dds)
  r <- results(dds, contrast=c("disease","MDD","Control"), alpha=0.05)
  r <- lfcShrink(dds, coef="disease_MDD_vs_Control", res=r, type="apeglm")
  r <- as.data.frame(r); r$gene <- rownames(r); r$symbol <- sym[r$gene]
  cat(sprintf("[%s] sig (padj<0.05 & |lfc|>0.3): %d\n", label, sum(r$padj<0.05 & abs(r$log2FoldChange)>0.3, na.rm=TRUE)))
  r
}

cat("\n===== QUALITY-ADJUSTED POOLED =====\n")
r_qual <- runDE(cnt, col, "quality-adjusted", TRUE)
r_noq  <- runDE(cnt, col, "NO quality adj",   FALSE)
write.table(r_qual[order(r_qual$padj), c("gene","symbol","log2FoldChange","lfcSE","pvalue","padj")],
            "deseq_pooled_qualityAdj_MDDvsControl.tsv", sep="\t", quote=FALSE, row.names=FALSE)

cat("\n===== IEG SURVIVAL (padj before vs after quality adj) =====\n")
ieg_ids <- names(sym)[sym %in% IEG]
tab <- data.frame(symbol=sym[ieg_ids],
  lfc_noq=r_noq$log2FoldChange[match(ieg_ids,r_noq$gene)], padj_noq=r_noq$padj[match(ieg_ids,r_noq$gene)],
  lfc_qual=r_qual$log2FoldChange[match(ieg_ids,r_qual$gene)], padj_qual=r_qual$padj[match(ieg_ids,r_qual$gene)])
tab <- tab[!is.na(tab$padj_noq) | !is.na(tab$padj_qual), ]
tab <- tab[order(tab$padj_qual), ]
print(format(tab, digits=3), row.names=FALSE)
cat(sprintf("\nIEG sig BEFORE: %d  |  AFTER quality adj: %d\n",
    sum(tab$padj_noq<0.05,na.rm=TRUE), sum(tab$padj_qual<0.05,na.rm=TRUE)))

cat("\n===== LEAVE-ONE-DATASET-OUT (quality-adjusted) =====\n")
core_qual <- r_qual$gene[which(r_qual$padj<0.05 & abs(r_qual$log2FoldChange)>0.3)]
cat("full core:", length(core_qual), "genes\n")
for (drop in levels(col$study)) {
  sub_m <- col[col$study!=drop, ]; sub_c <- cnt[, rownames(sub_m)]; sub_m$study <- droplevels(sub_m$study)
  r <- tryCatch(runDE(sub_c, sub_m, paste("drop",drop), TRUE), error=function(e){cat(" failed\n");NULL})
  if (!is.null(r)) {
    core_sub <- r$gene[which(r$padj<0.05 & abs(r$log2FoldChange)>0.3)]
    kf <- if(length(core_qual)) round(100*length(intersect(core_qual,core_sub))/length(core_qual)) else NA
    cat(sprintf("  drop %s: %d sig; retains %d%% of core; IEGs sig: %d\n",
        drop, length(core_sub), kf, sum(sym[core_sub] %in% IEG)))
  }
}
saveRDS(list(qual=r_qual, noq=r_noq, ieg=tab), "deseq_qualityAdj_results.rds")
cat("\nWROTE outputs. Trust IEG result IF survives quality adj AND leave-one-out.\n")
