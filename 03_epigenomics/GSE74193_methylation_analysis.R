suppressPackageStartupMessages({
  library(GEOquery); library(minfi); library(limma)
  library(DMRcate); library(missMethyl)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(IlluminaHumanMethylation450kmanifest)
  library(ggplot2); library(pheatmap)
})
outdir <- "./results"; plotdir <- "./plots"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
dir.create(plotdir, showWarnings=FALSE, recursive=TRUE)

message("Step 1: Downloading GSE74193 (526 samples)...")
if (length(list.files("./idat", pattern="\\.idat\\.gz$")) == 0) {
  getGEOSuppFiles("GSE74193", makeDirectory=FALSE, baseDir="./idat")
  for (t in list.files("./idat", pattern="\\.tar$", full.names=TRUE)) {
    message("  Extracting: ", basename(t))
    untar(t, exdir="./idat")
  }
} else {
  message("  idat files already present, skipping download")
}

message("Step 2: Reading data...")
gse     <- getGEO("GSE74193", GSEMatrix=TRUE, getGPL=FALSE)
pdata   <- pData(gse[[1]])
targets <- read.metharray.exp(base="./idat", recursive=TRUE)
sampleNames(targets) <- gsub(".*/(GSM[0-9]+).*","\\1",sampleNames(targets))
idx <- match(sampleNames(targets), rownames(pdata))
targets$diagnosis <- pdata$characteristics_ch1[idx]
targets$age       <- suppressWarnings(as.numeric(gsub("[^0-9.]","",pdata$characteristics_ch1.1[idx])))
targets$sex       <- gsub(".*: ","",pdata$characteristics_ch1.2[idx])
targets$period    <- ifelse(targets$age < 0, "Prenatal",
                    ifelse(targets$age < 18, "Postnatal_child","Postnatal_adult"))
message("  Total samples: ", ncol(targets))
print(table(targets$diagnosis))

message("Step 3: QC...")
detP   <- detectionP(targets)
failed <- colMeans(detP > 0.01) > 0.01
message("  Failing QC: ", sum(failed))
if(any(failed)){ targets <- targets[,!failed]; detP <- detP[,!failed] }

message("Step 4: Normalising...")
mSet.norm <- preprocessQuantile(targets)
mSet.norm <- dropLociWithSnps(mSet.norm)
keep      <- rowSums(detP[rownames(mSet.norm),] < 0.01) == ncol(mSet.norm)
mSet.norm <- mSet.norm[keep,]
mVals <- getM(mSet.norm); betaVals <- getBeta(mSet.norm)
message("  Probes retained: ", nrow(mSet.norm))
save(mVals, betaVals, file=file.path(outdir,"mVals_betaVals.RData"))

message("Step 5: MDS plot...")
is_scz <- grepl("schizophrenia",targets$diagnosis,ignore.case=TRUE)
pdf(file.path(plotdir,"mds_plot.pdf"),width=9,height=8)
plotMDS(mVals, col=ifelse(is_scz,"tomato","steelblue"),
  pch=c(Prenatal=17,Postnatal_child=15,Postnatal_adult=16)[targets$period],
  cex=0.8, main="MDS — GSE74193 (SCZ vs Control, PFC)")
legend("topright", legend=c("SCZ","Control","Prenatal","Child","Adult"),
  col=c("tomato","steelblue","black","black","black"),
  pch=c(16,16,17,15,16),bty="n",cex=0.8)
dev.off()

message("Step 6: Differential methylation...")
diagnosis <- factor(is_scz,levels=c(FALSE,TRUE),labels=c("Control","SCZ"))
age <- targets$age; sex <- factor(targets$sex); period <- factor(targets$period)
design   <- model.matrix(~0 + diagnosis + age + sex + period)
colnames(design)[1:2] <- c("Control","SCZ")
cont_mat <- makeContrasts(SCZ-Control, levels=design)
fit2     <- eBayes(contrasts.fit(lmFit(mVals,design), cont_mat))
top_cpgs <- topTable(fit2, num=Inf, coef=1, sort.by="p")
write.csv(top_cpgs, file=file.path(outdir,"DMCpGs_SCZ_vs_Control.csv"))
message("  Sig CpGs (FDR<0.05): ", sum(top_cpgs$adj.P.Val<0.05,na.rm=TRUE))

message("Step 7: DMRcate...")
cpg_anno <- cpg.annotate(object=mVals,datatype="array",what="M",
  analysis.type="differential",design=design,contrasts=TRUE,
  cont.matrix=cont_mat,coef="SCZ - Control",arraytype="450K")
dmr_ranges <- extractRanges(dmrcate(cpg_anno,lambda=1000,C=2))
write.csv(as.data.frame(dmr_ranges),
  file=file.path(outdir,"DMRs_SCZ_vs_Control.csv"),row.names=FALSE)
message("  DMRs: ", length(dmr_ranges))

message("Step 8: GO enrichment...")
sig_cpgs <- rownames(top_cpgs)[top_cpgs$adj.P.Val<0.05]
if(length(sig_cpgs)>10){
  go_res <- gometh(sig.cpg=sig_cpgs,all.cpg=rownames(top_cpgs),
                   plot.bias=FALSE,array.type="450K")
  write.csv(go_res[order(go_res$FDR),],
    file=file.path(outdir,"GO_enrichment.csv"),row.names=FALSE)
}

message("Step 9: Heatmap...")
top500 <- head(rownames(top_cpgs[order(top_cpgs$P.Value),]),500)
anno_col <- data.frame(Diagnosis=diagnosis,Sex=sex,Period=period,
                       row.names=colnames(betaVals))
anno_colors <- list(
  Diagnosis=c(Control="steelblue",SCZ="tomato"),
  Sex=c(M="gold",F="plum"),
  Period=c(Prenatal="darkgreen",Postnatal_child="orange",Postnatal_adult="purple"))
pdf(file.path(plotdir,"heatmap_top500.pdf"),width=20,height=10)
pheatmap(betaVals[top500,],annotation_col=anno_col,
  annotation_colors=anno_colors,
  show_rownames=FALSE,show_colnames=FALSE,
  clustering_distance_rows="correlation",
  clustering_distance_cols="correlation",
  color=colorRampPalette(c("navy","white","firebrick3"))(100),
  main="Top 500 Variable CpGs — GSE74193 (SCZ vs Control, PFC)")
dev.off()

message("Step 10: Saving...")
save(mSet.norm,mVals,betaVals,top_cpgs,dmr_ranges,
     file=file.path(outdir,"GSE74193_full_analysis.RData"))
message("=== GSE74193 COMPLETE ===")
