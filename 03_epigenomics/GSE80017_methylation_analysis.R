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

message("Step 1: idat files already extracted — skipping download")
message("  idat files found: ", length(list.files("./idat", pattern="\\.idat\\.gz$")))

message("Step 2: Reading data...")
gse     <- getGEO("GSE80017", GSEMatrix=TRUE, getGPL=FALSE)
pdata   <- pData(gse[[1]])
targets <- read.metharray.exp(base="./idat", recursive=TRUE)
sampleNames(targets) <- gsub(".*/(GSM[0-9]+).*","\\1",sampleNames(targets))
idx <- match(sampleNames(targets), rownames(pdata))
targets$diagnosis <- pdata$characteristics_ch1[idx]
targets$tissue    <- pdata$characteristics_ch1.1[idx]
targets$age       <- suppressWarnings(as.numeric(gsub("[^0-9.]","",pdata$characteristics_ch1.2[idx])))
targets$sex       <- gsub(".*: ","",pdata$characteristics_ch1.3[idx])
message("  Total samples: ", ncol(targets))
print(table(targets$diagnosis, targets$tissue))

message("Step 3: QC...")
detP   <- detectionP(targets)
failed <- colMeans(detP > 0.01) > 0.01
if(any(failed)){ targets <- targets[,!failed]; detP <- detP[,!failed] }
message("  Samples passing QC: ", ncol(targets))

message("Step 4: Normalising...")
mSet.norm <- preprocessQuantile(targets)
mSet.norm <- dropLociWithSnps(mSet.norm)
keep      <- rowSums(detP[rownames(mSet.norm),] < 0.01) == ncol(mSet.norm)
mSet.norm <- mSet.norm[keep,]
mVals <- getM(mSet.norm); betaVals <- getBeta(mSet.norm)
message("  Probes retained: ", nrow(mSet.norm))
save(mVals, betaVals, file=file.path(outdir,"mVals_betaVals.RData"))

message("Step 5: MDS plot...")
is_asd  <- grepl("autism|asd",targets$diagnosis,ignore.case=TRUE)
regions <- gsub(".*: ","",targets$tissue)
reg_col <- c("PFC"="purple","TC"="darkorange","CB"="darkgreen")
pdf(file.path(plotdir,"mds_by_region.pdf"),width=10,height=8)
plotMDS(mVals, col=reg_col[regions], pch=ifelse(is_asd,16,1), cex=1.2,
        main="MDS — GSE80017 (shape=diagnosis, colour=region)")
legend("topright",legend=c(names(reg_col),"ASD","Control"),
       col=c(reg_col,"black","black"),pch=c(15,15,15,16,1),bty="n")
dev.off()

message("Step 6: Differential methylation (all regions)...")
diagnosis <- factor(is_asd,levels=c(FALSE,TRUE),labels=c("Control","ASD"))
age <- targets$age; sex <- factor(targets$sex); tissue <- factor(regions)
design   <- model.matrix(~0 + diagnosis + age + sex + tissue)
colnames(design)[1:2] <- c("Control","ASD")
cont_mat <- makeContrasts(ASD-Control, levels=design)
fit2     <- eBayes(contrasts.fit(lmFit(mVals,design), cont_mat))
top_cpgs <- topTable(fit2, num=Inf, coef=1, sort.by="p")
write.csv(top_cpgs, file=file.path(outdir,"DMCpGs_ASD_vs_Control_allRegions.csv"))
message("  Sig CpGs (FDR<0.05): ", sum(top_cpgs$adj.P.Val<0.05,na.rm=TRUE))

message("Step 6b: Per-region analysis...")
for (reg in unique(regions)) {
  idx_reg <- regions == reg
  mV_r    <- mVals[,idx_reg]
  diag_r  <- diagnosis[idx_reg]
  age_r   <- age[idx_reg]; sex_r <- factor(sex[idx_reg])
  des_r   <- model.matrix(~0 + diag_r + age_r + sex_r)
  colnames(des_r)[1:2] <- c("Control","ASD")
  cm_r    <- makeContrasts(ASD-Control, levels=des_r)
  fit_r   <- eBayes(contrasts.fit(lmFit(mV_r,des_r), cm_r))
  top_r   <- topTable(fit_r, num=Inf, coef=1, sort.by="p")
  write.csv(top_r, file=file.path(outdir,paste0("DMCpGs_ASD_vs_Control_",reg,".csv")))
  message("  ", reg, " sig CpGs: ", sum(top_r$adj.P.Val<0.05,na.rm=TRUE))
}

message("Step 7: DMRcate...")
cpg_anno  <- cpg.annotate(object=mVals,datatype="array",what="M",
  analysis.type="differential",design=design,contrasts=TRUE,
  cont.matrix=cont_mat,coef="ASD - Control",arraytype="450K")
dmr_ranges <- extractRanges(dmrcate(cpg_anno,lambda=1000,C=2))
write.csv(as.data.frame(dmr_ranges),
  file=file.path(outdir,"DMRs_ASD_vs_Control_allRegions.csv"),row.names=FALSE)
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
top500   <- head(rownames(top_cpgs[order(top_cpgs$P.Value),]),500)
anno_col <- data.frame(Diagnosis=diagnosis,Region=tissue,Sex=sex,
                       row.names=colnames(betaVals))
anno_colors <- list(
  Diagnosis=c(Control="steelblue",ASD="tomato"),
  Region=c(PFC="purple",TC="darkorange",CB="darkgreen"),
  Sex=c(M="gold",F="plum"))
pdf(file.path(plotdir,"heatmap_top500.pdf"),width=18,height=10)
pheatmap(betaVals[top500,],annotation_col=anno_col,
  annotation_colors=anno_colors,
  show_rownames=FALSE,show_colnames=FALSE,
  clustering_distance_rows="correlation",
  clustering_distance_cols="correlation",
  color=colorRampPalette(c("navy","white","firebrick3"))(100),
  main="Top 500 Variable CpGs — GSE80017 (ASD vs Control, 3 brain regions)")
dev.off()

message("Step 10: Saving...")
save(mSet.norm,mVals,betaVals,top_cpgs,dmr_ranges,
     file=file.path(outdir,"GSE80017_full_analysis.RData"))
message("=== GSE80017 COMPLETE ===")
