## Fixed-effect meta-analysis across clean, unique-donor SCZ cohorts.
## CA1 (26) + GSE80655 DLPFC (48) + GSE78936 BA11 (28) = 102 donors.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages({library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)})
logf <- function(...) cat(sprintf(...),"\n")
CS  <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/counts_scz"
OUT <- file.path(CS,"META"); dir.create(OUT, showWarnings=FALSE)

grab <- function(path, lab) {
  if (!file.exists(path)) { logf("MISSING: %s", path); return(NULL) }
  d <- read.csv(path, stringsAsFactors=FALSE)
  gcol <- if ("gene" %in% names(d)) "gene" else names(d)[1]
  d <- data.frame(gene = sub("\\..*","", d[[gcol]]),
                  beta = d$log2FoldChange, se = d$lfcSE, padj = d$padj,
                  stringsAsFactors=FALSE)
  d <- d[!is.na(d$beta) & !is.na(d$se) & d$se > 0 & is.finite(d$se), ]
  d <- d[!duplicated(d$gene), ]
  logf("  %-14s %6d genes | %d at FDR<0.05", lab, nrow(d), sum(d$padj<0.05, na.rm=TRUE))
  setNames(d, c("gene", paste0(c("beta.","se.","padj."), lab)))
}

logf("=========== INPUT COHORTS ===========")
a <- grab(file.path(CS,"HIPP_donor","dec_ca1.csv"),                    "CA1")
b <- grab(file.path(CS,"PERCOHORT_prua","GSE80655_D2_full.csv"),       "DLPFC")
c <- grab(file.path(CS,"PERCOHORT_prua","GSE78936_BA11_F1_full.csv"),  "BA11")
stopifnot(!is.null(a), !is.null(b), !is.null(c))

m <- Reduce(function(x,y) merge(x,y,by="gene"), list(a,b,c))
logf("\ngenes tested in all three cohorts: %d", nrow(m))

B  <- as.matrix(m[,grep("^beta\\.",colnames(m))])
SE <- as.matrix(m[,grep("^se\\.",  colnames(m))])
W  <- 1/SE^2
m$beta_FE <- rowSums(W*B)/rowSums(W)
m$se_FE   <- sqrt(1/rowSums(W))
m$z_FE    <- m$beta_FE/m$se_FE
m$p_FE    <- 2*pnorm(-abs(m$z_FE))
m$padj_FE <- p.adjust(m$p_FE, "BH")

## heterogeneity
Q <- rowSums(W*(B - m$beta_FE)^2)
m$Q <- Q; m$I2 <- pmax(0, 100*(Q-2)/Q)
m$p_Q <- pchisq(Q, df=2, lower.tail=FALSE)
m$concordant <- rowSums(sign(B) == sign(m$beta_FE))

sym <- AnnotationDbi::mapIds(org.Hs.eg.db, m$gene, "SYMBOL","ENSEMBL", multiVals="first")
m$symbol <- ifelse(is.na(sym), m$gene, sym)
m <- m[order(m$p_FE), ]
write.csv(m, file.path(OUT,"SCZ_meta_FE_full.csv"), row.names=FALSE)

logf("\n=========== META-ANALYSIS RESULT (102 donors) ===========")
f <- function(l) sum(m$padj_FE < 0.05 & abs(m$beta_FE) >= l, na.rm=TRUE)
logf("  FDR<0.05: %d | >=0.3: %d | >=0.5: %d | >=1.0: %d", f(0), f(0.3), f(0.5), f(1))
sig <- subset(m, padj_FE < 0.05)
if (nrow(sig)) {
  logf("  of these, all-3-cohort direction concordant: %d/%d", sum(sig$concordant==3), nrow(sig))
  logf("  significant heterogeneity (p_Q<0.05): %d", sum(sig$p_Q < 0.05, na.rm=TRUE))
  write.csv(sig, file.path(OUT,"SCZ_meta_FE_sig.csv"), row.names=FALSE)
  print(head(sig[,c("symbol","beta_FE","se_FE","p_FE","padj_FE","concordant","I2")], 25))
} else {
  logf("  no genes reach FDR<0.05")
  logf("\n  top 15 by p_FE (nominal only):")
  print(head(m[,c("symbol","beta_FE","p_FE","padj_FE","concordant")], 15))
  logf("\n  nominal p_FE<0.05: %d | of which 3/3 concordant: %d",
       sum(m$p_FE<0.05), sum(m$p_FE<0.05 & m$concordant==3))
}

logf("\n=========== CROSS-COHORT AGREEMENT ===========")
nm <- c("CA1","DLPFC","BA11")
for (i in 1:2) for (j in (i+1):3)
  logf("  beta correlation %-6s vs %-6s : r = %+.3f", nm[i], nm[j],
       cor(B[,i], B[,j], use="complete.obs"))
logf("  genes with all three same sign: %d / %d (%.1f%%; 25%% expected by chance)",
     sum(m$concordant==3), nrow(m), 100*mean(m$concordant==3))

logf("\n=========== GSEA ON THE META-ANALYSED RANKING ===========")
m$entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, m$gene, "ENTREZID","ENSEMBL", multiVals="first")
r <- m[!is.na(m$entrez) & !duplicated(m$entrez), ]
gl <- sort(setNames(sign(r$beta_FE) * -log10(pmax(r$p_FE, .Machine$double.xmin)), r$entrez),
           decreasing=TRUE)
logf("  ranked list: %d genes (ASD-harmonised statistic)", length(gl))
for (ont in c("BP","MF","CC")) {
  g <- tryCatch(gseGO(gl, OrgDb=org.Hs.eg.db, keyType="ENTREZID", ont=ont,
        pAdjustMethod="BH", pvalueCutoff=0.05, verbose=FALSE), error=function(e) NULL)
  n <- if (is.null(g)) 0 else nrow(g@result)
  logf("  GSEA GO-%s : %d terms", ont, n)
  if (n) { write.csv(g@result, file.path(OUT,sprintf("META_GSEA_GO%s.csv",ont)), row.names=FALSE)
    print(head(g@result[order(g@result$p.adjust), c("Description","NES","p.adjust")], 10)) }
}
k <- tryCatch(gseKEGG(gl, organism="hsa", pAdjustMethod="BH", pvalueCutoff=0.05, verbose=FALSE),
              error=function(e) NULL)
n <- if (is.null(k)) 0 else nrow(k@result)
logf("  GSEA KEGG : %d terms", n)
if (n) { write.csv(k@result, file.path(OUT,"META_GSEA_KEGG.csv"), row.names=FALSE)
  print(head(k@result[order(k@result$p.adjust), c("Description","NES","p.adjust")], 10)) }
logf("\nOutputs -> %s", OUT)
