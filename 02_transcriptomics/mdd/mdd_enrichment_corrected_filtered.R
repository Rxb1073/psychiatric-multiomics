
#!/usr/bin/env Rscript
# =====================================================================
# MDD enrichment, corrected -- and the integration-ready KEGG output.
#
# FOUR GAPS CLOSED HERE
#
#  1  WRONG INPUT. The earlier enrichment ran on the three-cohort
#     donor-level table, before GSE80655 was excluded. This uses the
#     two-cohort primary analysis (GSE101521 + GSE102556, 107 donors).
#
#  2  WRONG RANKING STATISTIC. The earlier run ranked on the Wald
#     statistic. The SCZ chapter states GSEA used
#         sign(log2FC) x -log10(p)
#     "the statistic used in the ASD arm". This uses that, so all three
#     arms are comparable. The Wald ranking is computed alongside for
#     reference only.
#
#  3  LENGTH DIAGNOSTIC NEVER RUN. GSEA previously returned 1,486 GO BP
#     terms -- about a fifth of the ontology -- with ribosomal genes UP
#     and extracellular matrix / collagen / angiogenesis DOWN. Ribosomal
#     transcripts are short and abundant; ECM and collagen transcripts are
#     among the longest in the genome. That is the shape of a transcript-
#     length axis, not a pathway result. The ASD arm returned the SAME
#     ECM/collagen/angiogenesis terms, which makes a shared technical
#     origin more likely than shared disease biology.
#     This correlates the ranking against length and re-runs GSEA on a
#     length-residualised ranking.
#
#  4  INTEGRATION OUTPUT. Writes a clean KEGG table for the pathway-level
#     convergence against the microbiome's MetaCyc/KEGG pathways, flagging
#     which terms survive length residualisation. Terms that do NOT survive
#     must not be used for convergence -- a technical axis converging with
#     a microbiome pathway would be a false positive.
#
# READ-ONLY. Writes to MDD/audit/enrichment2/.
# =====================================================================

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi) })

MDD <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD"
P2  <- file.path(MDD, "audit/filtered/mirror/primary2")
OUT <- file.path(MDD, "audit/enrichment2"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
GTF <- "/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf"
hdr <- function(x) cat("\n", strrep("=",76), "\n", x, "\n", strrep("=",76), "\n", sep="")
set.seed(42)

# ============================ 0. INPUT ============================
hdr("0. INPUT -- two-cohort primary analysis")
f <- file.path(P2, "pooled_2cohort.tsv")
if (!file.exists(f)) stop("missing ", f, " -- run mdd_primary_2cohort.R first")
res <- read.delim(f, check.names=FALSE)
gc0 <- if ("gene" %in% names(res)) "gene" else names(res)[1]
res$ensembl <- sub("\\..*", "", res[[gc0]])
res <- res[!is.na(res$pvalue) & !is.na(res$log2FoldChange), ]
cat(sprintf("  %s\n  %d genes with a p-value\n", f, nrow(res)))

# ============================ 1. RANKINGS ============================
hdr("1. RANKING STATISTIC")
cat("  primary : sign(log2FC) x -log10(p)   [SCZ and ASD convention]\n")
cat("  compare : Wald statistic             [what the earlier run used]\n")
eps <- .Machine$double.xmin
res$rank_scz  <- sign(res$log2FoldChange) * -log10(pmax(res$pvalue, eps))
res$rank_wald <- res$stat
r12 <- cor(res$rank_scz, res$rank_wald, method="spearman", use="complete.obs")
cat(sprintf("\n  the two rankings correlate rho = %+.3f\n", r12))
cat("  (high correlation means the earlier result was not badly wrong, but the\n")
cat("   SCZ-convention ranking is the one to report)\n")

map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,
        keys=unique(res$ensembl), keytype="ENSEMBL", columns=c("ENTREZID","SYMBOL")))
map <- map[!is.na(map$ENTREZID) & !duplicated(map$ENSEMBL), ]
res <- merge(res, map, by.x="ensembl", by.y="ENSEMBL")
cat(sprintf("  mapped to Entrez: %d genes\n", nrow(res)))

# ============================ 2. GENE LENGTH ============================
hdr("2. LENGTH DIAGNOSTIC")
cat("  summing exon widths per gene from the GTF (awk)\n")
len <- local({
  cmd <- sprintf(
   "awk -F'\\t' '$3==\"exon\"{match($9,/gene_id \"[^\"]+\"/); g=substr($9,RSTART+9,RLENGTH-10); L[g]+=$5-$4+1} END{for(k in L) print k\"\\t\"L[k]}' %s",
   GTF)
  tab <- read.delim(pipe(cmd), header=FALSE, col.names=c("gene","len"), stringsAsFactors=FALSE)
  setNames(tab$len, sub("\\..*", "", tab$gene)) })
res$len <- len[res$ensembl]
res <- res[is.finite(res$len) & res$len > 0, ]
res$log_len <- log10(res$len)
cat(sprintf("  lengths for %d genes; median %.0f bp\n", nrow(res), median(res$len)))

rl <- cor(res$rank_scz, res$log_len, method="spearman")
rb <- cor(res$rank_scz, log10(res$baseMean+1), method="spearman")
cat(sprintf("\n  rho(ranking, log10 length)    = %+.3f%s\n", rl,
            ifelse(abs(rl)>0.20, "   <-- STRONG length axis", "")))
cat(sprintf("  rho(ranking, log10 abundance) = %+.3f\n", rb))

fam <- function(pat, lab) {
  k <- grepl(pat, res$SYMBOL)
  if (sum(k) < 5) return(invisible())
  cat(sprintf("  %-24s n=%4d  median length %7.0f bp  median rank stat %+6.2f\n",
              lab, sum(k), median(res$len[k]), median(res$rank_scz[k]))) }
cat("\n  the families driving the earlier GSEA result:\n")
fam("^RP[LS][0-9]",       "ribosomal proteins")
fam("^COL[0-9]",          "collagens")
fam("^(FOS|JUN|EGR|NR4A|NPAS4|ATF3|BTG2|ARC|DUSP)[0-9]?$", "IEGs")

fit <- lm(rank_scz ~ log_len, data=res)
res$rank_resid <- residuals(fit)
cat(sprintf("\n  length explains %.1f%% of the variance in the ranking statistic\n",
            100*summary(fit)$r.squared))

mk <- function(v) { r <- v; names(r) <- res$ENTREZID
                    sort(r[!duplicated(names(r))], decreasing=TRUE) }
RK  <- mk(res$rank_scz)
RKR <- mk(res$rank_resid)

# ============================ 3. GSEA ============================
run_gsea <- function(rk, tag) {
  out <- list()
  for (ont in c("BP","MF","CC")) {
    g <- try(gseGO(rk, OrgDb=org.Hs.eg.db, ont=ont, keyType="ENTREZID",
                   minGSSize=15, maxGSSize=500, pvalueCutoff=0.05,
                   eps=0, verbose=FALSE), silent=TRUE)
    if (!inherits(g,"try-error")) out[[paste0("GO_",ont)]] <- as.data.frame(g) }
  k <- try(gseKEGG(rk, organism="hsa", minGSSize=15, maxGSSize=500,
                   pvalueCutoff=0.05, eps=0, verbose=FALSE), silent=TRUE)
  if (!inherits(k,"try-error")) out$KEGG <- as.data.frame(k)
  for (n in names(out)) if (nrow(out[[n]]))
    write.csv(out[[n]], file.path(OUT, sprintf("GSEA_%s_%s.csv", n, tag)), row.names=FALSE)
  out }

show <- function(d, n=10) {
  if (is.null(d) || !nrow(d)) { cat("      (nothing significant)\n"); return(invisible()) }
  d <- d[order(d$p.adjust), ]
  for (i in seq_len(min(n, nrow(d))))
    cat(sprintf("      %-54s NES %+5.2f  %.1e\n",
                substr(d$Description[i],1,54), d$NES[i], d$p.adjust[i])) }

hdr("3. GSEA -- SCZ-convention ranking, minGSSize=15 (matching SCZ/ASD)")
G1 <- run_gsea(RK, "primary")
for (n in names(G1)) {
  cat(sprintf("\n  --- %s: %d terms at FDR<0.05 ---\n", n, sum(G1[[n]]$p.adjust<0.05)))
  show(G1[[n]]) }

hdr("4. GSEA -- LENGTH-RESIDUALISED ranking")
cat("  Terms that vanish here were driven by transcript length, not biology.\n")
G2 <- run_gsea(RKR, "length_residualised")
for (n in names(G2)) {
  cat(sprintf("\n  --- %s: %d terms at FDR<0.05 ---\n", n, sum(G2[[n]]$p.adjust<0.05)))
  show(G2[[n]]) }

# ============================ 5. WHAT SURVIVES ============================
hdr("5. WHICH TERMS SURVIVE LENGTH RESIDUALISATION?")
for (n in intersect(names(G1), names(G2))) {
  a <- G1[[n]]$Description[G1[[n]]$p.adjust<0.05]
  b <- G2[[n]]$Description[G2[[n]]$p.adjust<0.05]
  cat(sprintf("\n  %-8s raw %4d  residualised %4d  surviving %4d (%.0f%%)\n",
              n, length(a), length(b), length(intersect(a,b)),
              ifelse(length(a), 100*length(intersect(a,b))/length(a), 0)))
  lost <- setdiff(a,b)
  if (length(lost)) {
    cat("    examples LOST (length-driven):\n")
    for (x in head(lost, 5)) cat(sprintf("      %s\n", substr(x,1,60))) }
  kept <- intersect(a,b)
  if (length(kept)) {
    cat("    examples SURVIVING:\n")
    for (x in head(kept, 5)) cat(sprintf("      %s\n", substr(x,1,60))) }
}

# ============================ 6. KEGG FOR INTEGRATION ============================
hdr("6. INTEGRATION-READY KEGG TABLE")
if (!is.null(G1$KEGG) && nrow(G1$KEGG)) {
  k1 <- G1$KEGG[G1$KEGG$p.adjust<0.05, c("ID","Description","NES","pvalue","p.adjust","setSize")]
  surv <- if (!is.null(G2$KEGG)) G2$KEGG$Description[G2$KEGG$p.adjust<0.05] else character(0)
  k1$survives_length_adjustment <- k1$Description %in% surv
  k1$direction <- ifelse(k1$NES<0, "down_in_MDD", "up_in_MDD")
  k1 <- k1[order(k1$p.adjust), ]
  write.csv(k1, file.path(OUT,"KEGG_for_integration.csv"), row.names=FALSE)
  cat(sprintf("  %d KEGG pathways at FDR<0.05; %d survive length adjustment\n",
              nrow(k1), sum(k1$survives_length_adjustment)))
  cat("\n  *** USE ONLY THE SURVIVING PATHWAYS FOR CONVERGENCE ***\n")
  cat("  A length-driven pathway converging with a microbiome pathway is a\n")
  cat("  false positive, not a finding.\n\n")
  cat(sprintf("  %-46s %6s %10s %s\n","pathway","NES","FDR","survives"))
  for (i in seq_len(min(25, nrow(k1))))
    cat(sprintf("  %-46s %+6.2f %10.2e %s\n", substr(k1$Description[i],1,46),
                k1$NES[i], k1$p.adjust[i],
                ifelse(k1$survives_length_adjustment[i], "yes", "-")))
} else cat("  no KEGG result (needs network access to rest.kegg.jp)\n")

# ============================ 7. IEG POSITION ============================
hdr("7. WHERE DO THE IEGs SIT?")
ieg <- c("NPAS4","FOS","FOSB","EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3",
         "ATF3","BTG2","ARC","JUN","JUNB","DUSP1","DUSP6")
im <- res[res$SYMBOL %in% ieg, ]
im$rank_raw   <- match(as.character(im$ENTREZID), names(RK))
im$rank_resid <- match(as.character(im$ENTREZID), names(RKR))
im <- im[order(im$rank_raw, decreasing=TRUE), ]
cat(sprintf("  %-8s %10s %10s %10s %9s\n","symbol","rank raw","rank resid","length","log2FC"))
for (i in seq_len(nrow(im)))
  cat(sprintf("  %-8s %10d %10d %10.0f %+9.3f\n", im$SYMBOL[i],
              im$rank_raw[i], im$rank_resid[i], im$len[i], im$log2FoldChange[i]))
cat(sprintf("\n  median IEG rank: raw %d, residualised %d, of %d\n",
            as.integer(median(im$rank_raw, na.rm=TRUE)),
            as.integer(median(im$rank_resid, na.rm=TRUE)), length(RK)))
cat("  A median rank near the bottom under BOTH rankings means the IEG result\n")
cat("  is independent of the length axis.\n")

hdr("DONE")
cat(sprintf("
  length-vs-ranking rho        %+.3f  (%.1f%% of variance)
  GO BP terms raw              %d
  GO BP terms residualised     %d
  IEG median rank raw          %d / %d
  IEG median rank residualised %d / %d
", rl, 100*summary(fit)$r.squared,
   ifelse(is.null(G1$GO_BP), 0, sum(G1$GO_BP$p.adjust<0.05)),
   ifelse(is.null(G2$GO_BP), 0, sum(G2$GO_BP$p.adjust<0.05)),
   as.integer(median(im$rank_raw, na.rm=TRUE)), length(RK),
   as.integer(median(im$rank_resid, na.rm=TRUE)), length(RKR)))
cat("\noutputs in", OUT, "\n")
sessionInfo()

