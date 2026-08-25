suppressPackageStartupMessages({
  library(GEOquery); library(minfi); library(limma)
  library(DMRcate); library(missMethyl)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(IlluminaHumanMethylation450kmanifest)
  library(ggplot2); library(pheatmap)
})
IDAT_DIR <- "/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-methylation/GSE80017/idat"
OUT_DIR  <- "/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-methylation/GSE80017/results"
PLOT_DIR <- "/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-methylation/GSE80017/plots"
dir.create(OUT_DIR,  showWarnings=FALSE, recursive=TRUE)
dir.create(PLOT_DIR, showWarnings=FALSE, recursive=TRUE)
message("Step 1: Checking idat files...")
grn_base <- gsub("_Grn[.]idat$","",list.files(IDAT_DIR, pattern="_Grn[.]idat$"))
red_base <- gsub("_Red[.]idat$","",list.files(IDAT_DIR, pattern="_Red[.]idat$"))
incomplete <- setdiff(grn_base, red_base)
if(length(incomplete)>0){
  message("  Removing incomplete: ", paste(incomplete,collapse=", "))
  file.remove(file.path(IDAT_DIR,paste0(incomplete,"_Grn.idat")))
}
message("  Complete pairs: ", length(intersect(grn_base,red_base)))
message("Step 2: Reading data...")
gse     <- getGEO("GSE80017", GSEMatrix=TRUE, getGPL=FALSE)
pdata   <- pData(gse[[1]])
targets <- read.metharray.exp(base=IDAT_DIR, recursive=TRUE)
sampleNames(targets) <- sub("_.*", "", sampleNames(targets))
idx <- match(sampleNames(targets), rownames(pdata))
targets$diagnosis <- pdata[["disease status:ch1"]][idx]
targets$tissue    <- gsub("tissue: ","",pdata$characteristics_ch1[idx])
targets$sex       <- gsub("gender: ","",pdata$characteristics_ch1.1[idx])
message("  Samples: ", ncol(targets))
print(table(targets$diagnosis, targets$tissue))
message("Step 3: QC...")
detP   <- detectionP(targets)
failed <- colMeans(detP > 0.01) > 0.01
if(any(failed)){ targets <- targets[,!failed]; detP <- detP[,!failed] }
message("  Passing QC: ", ncol(targets))
message("Step 4: Normalising...")
mSet.norm <- preprocessQuantile(targets)
mSet.norm <- dropLociWithSnps(mSet.norm)
keep      <- rowSums(detP[rownames(mSet.norm),] < 0.01) == ncol(mSet.norm)
mSet.norm <- mSet.norm[keep,]
mVals <- getM(mSet.norm); betaVals <- getBeta(mSet.norm)
message("  Probes retained: ", nrow(mSet.norm))
save(mVals, betaVals, file=file.path(OUT_DIR,"mVals_betaVals.RData"))
message("Step 5: MDS...")
is_asd  <- grepl("ASD|Autism", targets$diagnosis, ignore.case=TRUE)
regions <- targets$tissue
pdf(file.path(PLOT_DIR,"mds_by_region.pdf"),width=10,height=8)
plotMDS(mVals, pch=ifelse(is_asd,16,1), cex=1.2, main="MDS - GSE80017")
dev.off()
message("Step 6: Differential methylation...")
diagnosis <- factor(is_asd,levels=c(FALSE,TRUE),labels=c("Control","ASD"))
sex    <- factor(targets$sex)
tissue <- factor(regions)
design <- model.matrix(~0 + diagnosis)
if(length(levels(sex))>1)    design <- cbind(design, sex=as.numeric(sex))
if(length(levels(tissue))>1) design <- cbind(design, tissue=as.numeric(tissue))
colnames(design)[1:2] <- c("Control","ASD")
cont_mat <- makeContrasts(ASD-Control, levels=design)
fit2     <- eBayes(contrasts.fit(lmFit(mVals,design), cont_mat))
top_cpgs <- topTable(fit2, num=Inf, coef=1, sort.by="p")
write.csv(top_cpgs, file=file.path(OUT_DIR,"DMCpGs_ASD_vs_Control_allRegions.csv"))
message("  Sig CpGs (FDR<0.05): ", sum(top_cpgs$adj.P.Val<0.05,na.rm=TRUE))
message("Step 7: DMRcate...")
cpg_anno <- cpg.annotate(object=mVals,datatype="array",what="M",
  analysis.type="differential",design=design,contrasts=TRUE,
  cont.matrix=cont_mat,coef="ASD - Control",arraytype="450K")
dmr_ranges <- extractRanges(dmrcate(cpg_anno,lambda=1000,C=2))
write.csv(as.data.frame(dmr_ranges),file=file.path(OUT_DIR,"DMRs_ASD_vs_Control.csv"),row.names=FALSE)
message("  DMRs: ", length(dmr_ranges))
message("Step 8: GO...")
sig_cpgs <- rownames(top_cpgs)[top_cpgs$adj.P.Val<0.05]
if(length(sig_cpgs)>10){
  go_res <- gometh(sig.cpg=sig_cpgs,all.cpg=rownames(top_cpgs),plot.bias=FALSE,array.type="450K")
  write.csv(go_res[order(go_res$FDR),],file=file.path(OUT_DIR,"GO_enrichment.csv"),row.names=FALSE)
}
message("Step 9: Heatmap...")
top500 <- head(rownames(top_cpgs[order(top_cpgs$P.Value),]),500)
anno_col <- data.frame(Diagnosis=diagnosis, row.names=colnames(betaVals))
pdf(file.path(PLOT_DIR,"heatmap_top500.pdf"),width=18,height=10)
pheatmap(betaVals[top500,],annotation_col=anno_col,show_rownames=FALSE,show_colnames=FALSE,
  clustering_distance_rows="correlation",clustering_distance_cols="correlation",
  color=colorRampPalette(c("navy","white","firebrick3"))(100),main="Top 500 CpGs - GSE80017")
dev.off()
message("Step 10: Saving...")
save(mSet.norm,mVals,betaVals,top_cpgs,dmr_ranges,file=file.path(OUT_DIR,"GSE80017_full_analysis.RData"))
message("=== GSE80017 COMPLETE ===")
