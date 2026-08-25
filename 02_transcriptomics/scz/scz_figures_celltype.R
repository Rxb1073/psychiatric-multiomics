## Figures for the GABAergic finding + CAMERA panel + IEG medication
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(ggplot2); library(DESeq2); library(AnnotationDbi); library(org.Hs.eg.db)})
logf <- function(...) cat(sprintf(...),"\n")
NP  <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
CS  <- file.path(NP,"counts_scz"); SHD <- file.path(CS,"sheets")
CT  <- file.path(CS,"CELLTYPE_SETS")
OUT <- file.path(CS,"FIGURES_CELLTYPE"); dir.create(OUT, showWarnings=FALSE)
sv <- function(p,w=7,h=6){ ggsave(file.path(OUT,paste0(p,".pdf")),width=w,height=h)
                           ggsave(file.path(OUT,paste0(p,".png")),width=w,height=h,dpi=150) }

## ---- Fig 1: marker forest, all sets, controls marked ----
m <- read.csv(file.path(CT,"marker_genes_per_set.csv"), stringsAsFactors=FALSE)
v <- read.csv(file.path(CT,"frontal_pool_voom_full.csv"), row.names=1, stringsAsFactors=FALSE)
m$se <- abs(v[match(m$symbol, v$symbol),"logFC"] / v[match(m$symbol, v$symbol),"t"])
m$set <- factor(m$set, levels=c("GABAergic","CalbindinCR","excitatory","astrocyte",
                                "oligodendrocyte","IEG"))
m$lab <- ifelse(m$set %in% c("CalbindinCR","oligodendrocyte"),
                "negative control", as.character(m$set))
m <- m[order(m$set, m$logFC),]
m$symbol <- factor(m$symbol, levels=unique(m$symbol))
ggplot(m, aes(logFC, symbol, colour=lab)) +
  geom_vline(xintercept=0, linetype=2, colour="grey55") +
  geom_errorbarh(aes(xmin=logFC-1.96*se, xmax=logFC+1.96*se), height=0, alpha=.7) +
  geom_point(size=2) +
  facet_grid(set ~ ., scales="free_y", space="free_y", switch="y") +
  scale_colour_manual(values=c(GABAergic="#922B21", IEG="#7D6608",
    excitatory="#1F618D", astrocyte="#1E8449", "negative control"="grey55")) +
  theme_bw() + theme(strip.text.y.left=element_text(angle=0), legend.position="bottom") +
  labs(x=expression(log[2]~fold~change~(SCZ~vs~control)), y=NULL, colour=NULL,
    title="Cell-type marker expression, frontal cortex (76 unique donors)",
    subtitle="All ten GABAergic markers reduced; both control sets flat")
sv("Fig_C1_marker_forest", 8, 11)
logf("Fig 1: %d markers across %d sets", nrow(m), nlevels(m$set))

## ---- Fig 2: CAMERA panel ----
cam <- read.csv(file.path(CT,"CAMERA_celltype_sets.csv"), row.names=1)
cam$set <- rownames(cam); cam$logp <- -log10(cam$FDR)
cam$ctrl <- cam$set %in% c("CalbindinCR","oligodendrocyte")
cam <- cam[order(cam$logp),]; cam$set <- factor(cam$set, levels=cam$set)
ggplot(cam, aes(logp, set, fill=Direction)) +
  geom_col(aes(alpha=!ctrl), width=.65) +
  geom_vline(xintercept=-log10(0.05), linetype=2, colour="grey40") +
  geom_text(aes(label=sprintf("%s  (n=%d)  FDR=%.1e", Direction, NGenes, FDR)),
            hjust=-0.03, size=3) +
  scale_fill_manual(values=c(Down="#922B21", Up="#1E8449")) +
  scale_alpha_manual(values=c(`TRUE`=1, `FALSE`=0.35), guide="none") +
  xlim(0, max(cam$logp)*1.75) + theme_bw() +
  labs(x=expression(-log[10]~FDR), y=NULL, fill=NULL,
    title="Cell-type set testing (CAMERA, corrects inter-gene correlation)",
    subtitle="Faded bars are negative controls; dashed line FDR = 0.05")
sv("Fig_C2_camera", 9, 5)

## ---- Fig 3: IEG by medication status, GSE138082 CA1 ----
ln <- readLines(gzfile(file.path(SHD,"GSE138082_series_matrix.txt.gz"))); rows <- list()
for (l in ln[grepl("^!Sample_", ln)]) {
  p <- gsub('"','', strsplit(l,"\t")[[1]]); k <- sub("^!","",p[1])
  if (k=="Sample_title") rows$title <- p[-1]
  else if (k=="Sample_characteristics_ch1") {
    tg <- trimws(sub(":.*","",p[2])); rows[[tg]] <- trimws(sub("^[^:]*:","",p[-1])) }
}
md <- as.data.frame(rows, stringsAsFactors=FALSE)
md$gsm <- gsub('"','', strsplit(grep("^!Sample_geo_accession", ln, value=TRUE)[1],"\t")[[1]][-1])
rt <- read.csv(file.path(NP,"metadata/PRJNA574470_SraRunTable.csv"), check.names=FALSE)
gc <- grep("GEO_Accession|Sample.Name|Sample Name", names(rt), value=TRUE)[1]
ss <- read.delim(file.path(SHD,"GSE138082_samplesheet.tsv"), check.names=FALSE)
ss$gsm <- setNames(rt[[gc]], rt$Run)[ss$srr]
ss <- merge(ss, md[,c("gsm","title","hippocampal_subfield","drug")], by="gsm")
ss$grp <- ifelse(grepl("_SZ_",ss$title), paste0("SCZ ",ss$drug), "Control")
rownames(ss) <- ss$srr; ss <- ss[ss$hippocampal_subfield=="CA1",]
cm <- read.delim(file.path(NP,"PRJNA574470/counts/PRJNA574470_gene_counts_matrix.tsv"),
                 comment.char="#", check.names=FALSE)
g <- cm[[1]]; cm[[1]] <- NULL
cm <- cm[, !colnames(cm) %in% c("Chr","Start","End","Strand","Length"), drop=FALSE]
colnames(cm) <- regmatches(colnames(cm), regexpr("SRR[0-9]+", colnames(cm)))
cnt <- as.matrix(cm); rownames(cnt) <- g; storage.mode(cnt) <- "integer"
k <- intersect(colnames(cnt), rownames(ss)); cnt <- cnt[,k]; ss <- ss[k,]
nz <- counts(estimateSizeFactors(DESeqDataSetFromMatrix(cnt, ss, ~1)), normalized=TRUE)
sym <- AnnotationDbi::mapIds(org.Hs.eg.db, sub("\\..*","",rownames(nz)),
                             "SYMBOL","ENSEMBL", multiVals="first")
show <- c("NPAS4","FOSB","EGR1","EGR3","JUN","NR4A1")
d <- do.call(rbind, lapply(show, function(gg){
  i <- which(!is.na(sym) & sym==gg)[1]; if (is.na(i)) return(NULL)
  data.frame(gene=gg, expr=nz[i,], grp=ss$grp) }))
d$grp <- factor(d$grp, levels=c("Control","SCZ OFF","SCZ ON"))
d$gene <- factor(d$gene, levels=show)
ggplot(d, aes(grp, expr, fill=grp)) +
  geom_boxplot(outlier.shape=NA, alpha=.7, width=.6) +
  geom_jitter(width=.15, size=1.2, alpha=.7) +
  facet_wrap(~gene, scales="free_y", nrow=2) +
  scale_y_log10() +
  scale_fill_manual(values=c(Control="#5D8AA8","SCZ OFF"="#7D6608","SCZ ON"="#922B21")) +
  theme_bw() + guides(fill="none") +
  labs(x=NULL, y="normalised counts (log scale)",
    title="Immediate-early genes by medication status (GSE138082 CA1)",
    subtitle="13 controls, 7 unmedicated and 6 medicated donors; effects appear only in the medicated group")
sv("Fig_C3_IEG_medication", 9, 6)
logf("Outputs -> %s", OUT)
