
#!/usr/bin/env Rscript
# =====================================================================
# STEP 10, final piece -- THREE-WAY ENRICHMENT PANEL
#
# Supervisor's suggestion on the ASD arm:
#   "We can try all dysregulated combined, then those up only and
#    down regulated only as one panel of the figure"
#
# and, on the ASD result:
#   "Combined gives fewer but better-supported terms"
#
# That is the standard reason for splitting: mixing up- and down-regulated
# genes in one ORA can either dilute a directional signal or manufacture
# terms that contain both directions. Running all three shows which.
#
# Applied here to the two-cohort primary. Gene lists at the arm's usual
# thresholds (FDR<0.05, |log2FC|>0.3), with a nominal-p list as a
# sensitivity analysis because the FDR list is small.
#
# CAUTION ALREADY ESTABLISHED. The GSEA on this data returns ~1,143 GO BP
# terms with trivial effect sizes (ribosomal proteins shift by +0.027,
# collagens by -0.093, against -0.591 for the IEGs). Those broad terms
# are a real but minuscule global axis and must NOT be carried into the
# integration. ORA on a thresholded list is less prone to this, but the
# same caution applies: check the effect sizes behind any term before
# reporting it.
#
# READ-ONLY. Writes to MDD/audit/enrichment3/.
# =====================================================================

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi) })

MDD <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD"
P2  <- file.path(MDD,"audit/primary2")
OUT <- file.path(MDD,"audit/enrichment3"); dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
hdr <- function(x) cat("\n", strrep("=",76), "\n", x, "\n", strrep("=",76), "\n", sep="")
set.seed(42)
FDR <- 0.05; LFC <- 0.3

IEG <- c("NPAS4","FOS","FOSB","EGR1","EGR2","EGR3","NR4A1","NR4A2","NR4A3",
         "ATF3","BTG2","ARC","JUN","JUNB","DUSP1","DUSP6")

# ============================ input ============================
hdr("0. GENE LISTS")
res <- read.delim(file.path(P2,"pooled_2cohort.tsv"), check.names=FALSE)
gc0 <- if ("gene" %in% names(res)) "gene" else names(res)[1]
res$ensembl <- sub("\\..*","", res[[gc0]])
map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=unique(res$ensembl),
        keytype="ENSEMBL", columns=c("ENTREZID","SYMBOL")))
map <- map[!is.na(map$ENTREZID) & !duplicated(map$ENSEMBL), ]
res <- merge(res, map, by.x="ensembl", by.y="ENSEMBL")
universe <- as.character(unique(res$ENTREZID))
cat(sprintf("  %d genes mapped to Entrez (universe)\n", length(universe)))

sig  <- !is.na(res$padj) & res$padj < FDR & abs(res$log2FoldChange) > LFC
nom  <- !is.na(res$pvalue) & res$pvalue < 0.05 & abs(res$log2FoldChange) > LFC
LISTS <- list(
  FDR_all  = as.character(res$ENTREZID[sig]),
  FDR_down = as.character(res$ENTREZID[sig & res$log2FoldChange < 0]),
  FDR_up   = as.character(res$ENTREZID[sig & res$log2FoldChange > 0]),
  nom_all  = as.character(res$ENTREZID[nom]),
  nom_down = as.character(res$ENTREZID[nom & res$log2FoldChange < 0]),
  nom_up   = as.character(res$ENTREZID[nom & res$log2FoldChange > 0]))
for (n in names(LISTS)) cat(sprintf("  %-9s %5d genes\n", n, length(LISTS[[n]])))

nomN <- sum(nom); exp_chance <- 0.05 * sum(!is.na(res$pvalue))
cat(sprintf("\n  nominal list: %d genes against ~%.0f expected by chance (%.2fx)\n",
            nomN, exp_chance, nomN/exp_chance))
cat("  Roughly half of a nominal list this size is noise -- the ASD arm made\n")
cat("  the same point. Treat nominal results as exploratory.\n")

# ============================ ORA ============================
run_ora <- function(genes, tag) {
  if (length(genes) < 5) { cat(sprintf("  %-9s too few genes\n", tag)); return(NULL) }
  out <- list()
  for (ont in c("BP","MF","CC")) {
    e <- try(enrichGO(genes, OrgDb=org.Hs.eg.db, ont=ont, keyType="ENTREZID",
                      universe=universe, pvalueCutoff=0.05, qvalueCutoff=0.10,
                      readable=TRUE), silent=TRUE)
    if (!inherits(e,"try-error") && !is.null(e)) out[[paste0("GO_",ont)]] <- as.data.frame(e) }
  k <- try(enrichKEGG(genes, organism="hsa", universe=universe,
                      pvalueCutoff=0.05, qvalueCutoff=0.10), silent=TRUE)
  if (!inherits(k,"try-error") && !is.null(k)) out$KEGG <- as.data.frame(k)
  for (n in names(out)) if (nrow(out[[n]]))
    write.csv(out[[n]], file.path(OUT, sprintf("ORA_%s_%s.csv", n, tag)), row.names=FALSE)
  out }

hdr("1. ORA -- FDR-SIGNIFICANT LISTS")
R1 <- lapply(c("FDR_all","FDR_down","FDR_up"), function(n) run_ora(LISTS[[n]], n))
names(R1) <- c("FDR_all","FDR_down","FDR_up")

hdr("2. ORA -- NOMINAL LISTS (exploratory)")
R2 <- lapply(c("nom_all","nom_down","nom_up"), function(n) run_ora(LISTS[[n]], n))
names(R2) <- c("nom_all","nom_down","nom_up")

# ============================ the panel ============================
summ <- function(R, tag) {
  cat(sprintf("\n  %-10s %8s %8s %8s %8s\n", tag, "GO BP","GO MF","GO CC","KEGG"))
  for (n in names(R)) {
    r <- R[[n]]
    cat(sprintf("  %-10s", n))
    for (o in c("GO_BP","GO_MF","GO_CC","KEGG")) {
      d <- if (is.null(r)) NULL else r[[o]]
      cat(sprintf(" %8s", if (is.null(d)) "-" else sum(d$p.adjust<0.05))) }
    cat("\n") } }

hdr("3. THE PANEL -- TERM COUNTS BY DIRECTION")
summ(R1, "FDR lists"); summ(R2, "nominal")
cat("\n  Combined giving FEWER but stronger terms than the split lists is the\n")
cat("  expected pattern when the signal is directional -- it is the result the\n")
cat("  ASD arm reported. Combined giving MORE suggests terms containing genes\n")
cat("  moving in both directions, which are weaker evidence.\n")

top <- function(R, n, o="GO_BP", k=8) {
  d <- if (is.null(R[[n]])) NULL else R[[n]][[o]]
  cat(sprintf("\n  --- %s / %s ---\n", n, o))
  if (is.null(d) || !nrow(d)) { cat("      (nothing significant)\n"); return(invisible()) }
  d <- d[order(d$p.adjust), ]
  for (i in seq_len(min(k, nrow(d))))
    cat(sprintf("      %-52s %3s  %.1e\n", substr(d$Description[i],1,52),
                d$Count[i], d$p.adjust[i])) }

hdr("4. TOP TERMS")
for (n in names(R1)) top(R1, n, "GO_BP")
for (n in names(R1)) top(R1, n, "KEGG", 5)
cat("\n  nominal lists, GO BP:\n")
for (n in names(R2)) top(R2, n, "GO_BP", 6)

# ============================ IEG contribution ============================
hdr("5. HOW MUCH OF THIS IS THE IEG SET?")
cat("  The down list is dominated by immediate-early genes, so terms recovered\n")
cat("  from it will largely restate the set-level result. Quantifying that.\n")
ieg_ent <- as.character(res$ENTREZID[res$SYMBOL %in% IEG])
for (n in c("FDR_down","nom_down")) {
  g <- LISTS[[n]]
  cat(sprintf("\n  %-9s %d genes, %d are IEGs (%.0f%%)\n", n, length(g),
              length(intersect(g, ieg_ent)),
              100*length(intersect(g, ieg_ent))/max(length(g),1)))
  r <- if (n=="FDR_down") R1$FDR_down else R2$nom_down
  d <- if (is.null(r)) NULL else r$GO_BP
  if (!is.null(d) && nrow(d)) {
    d <- d[order(d$p.adjust), ]
    hit <- sapply(head(d$geneID, 5), function(x)
             sum(strsplit(x,"/")[[1]] %in% IEG))
    tot <- sapply(head(d$Count, 5), identity)
    cat("    IEG share of the top 5 terms:\n")
    for (i in seq_along(hit))
      cat(sprintf("      %-46s %d/%d genes\n", substr(d$Description[i],1,46), hit[i], tot[i]))
  }
  # ORA with the IEGs removed
  g2 <- setdiff(g, ieg_ent)
  if (length(g2) >= 5) {
    e <- try(enrichGO(g2, OrgDb=org.Hs.eg.db, ont="BP", keyType="ENTREZID",
                      universe=universe, pvalueCutoff=0.05, qvalueCutoff=0.10,
                      readable=TRUE), silent=TRUE)
    nn <- if (inherits(e,"try-error") || is.null(e)) 0 else sum(as.data.frame(e)$p.adjust<0.05)
    cat(sprintf("    with IEGs removed (%d genes): %d GO BP terms\n", length(g2), nn))
    if (nn) write.csv(as.data.frame(e), file.path(OUT, sprintf("ORA_GO_BP_%s_noIEG.csv", n)),
                      row.names=FALSE) }
}

# ============================ figure ============================
hdr("6. FIGURE")
FIG <- file.path(MDD,"audit/figures2"); dir.create(FIG, recursive=TRUE, showWarnings=FALSE)
mkpanel <- function(R, tag, fname) {
  dat <- list()
  for (n in names(R)) {
    d <- if (is.null(R[[n]])) NULL else R[[n]]$GO_BP
    if (!is.null(d) && nrow(d)) {
      d <- d[order(d$p.adjust), ]; d <- head(d, 8)
      dat[[n]] <- data.frame(term=d$Description, p=d$p.adjust, cnt=d$Count) } }
  if (!length(dat)) { cat("   nothing to plot for", tag, "\n"); return(invisible()) }
  for (fmt in c("png","pdf")) {
    f <- file.path(FIG, paste0(fname, ".", fmt))
    if (fmt=="png") png(f, width=13, height=4.6, units="in", res=300, type="cairo")
    else pdf(f, width=13, height=4.6)
    par(mfrow=c(1,length(dat)), mar=c(4.5,19,3,1))
    for (n in names(dat)) {
      d <- dat[[n]]; o <- order(d$p, decreasing=TRUE)
      barplot(-log10(d$p[o]), horiz=TRUE, names.arg=substr(d$term[o],1,44),
              las=1, cex.names=.72, border=NA,
              col=ifelse(grepl("down", n), "#C0392B", ifelse(grepl("up", n), "#2E86C1", "grey45")),
              xlab=expression(-log[10]~adjusted~p), main=n)
      abline(v=-log10(0.05), lty=2) }
    dev.off() }
  cat("   wrote", fname, ".png/.pdf\n") }
mkpanel(R1, "FDR", "F8_enrichment_threeway_FDR")
mkpanel(R2, "nominal", "F8b_enrichment_threeway_nominal")

hdr("HOW TO READ THIS")
cat("
  This is the panel requested for the ASD arm: all dysregulated genes,
  then down-only and up-only, side by side.

  Section 5 matters most. If the down list is mostly IEGs and the terms
  vanish once they are removed, this analysis restates the set-level result
  in ontology language rather than adding independent evidence. Say so --
  it is still worth showing, but it is not a second finding.

  The same caution as before: check effect sizes behind any term before
  carrying it into the integration. The broad GSEA terms in this data rest
  on 2-9 percent shifts and should be excluded from convergence.
")
cat("\noutputs in", OUT, "\n  figures in", FIG, "\n")
sessionInfo()

