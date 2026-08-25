## Donor-correct cell-type proportions for integration (76 frontal donors).
## Subsets the existing per-sample deconvolution to unique donors and collapses
## GSE80655 split runs. No re-deconvolution needed.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
logf <- function(...) cat(sprintf(...),"\n")
CS  <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq/counts_scz"
NP  <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/RNA_seq"
RT  <- "/rds/projects/e/elhamsak-dissertation-autism-2/MDD/transcriptomics/GSE80655/metadata/PRJNA319583_SraRunTable.csv"
DEC <- file.path(CS,"Deconvolution")
OUT <- file.path(CS,"INTEGRATION"); dir.create(OUT, showWarnings=FALSE)

d <- read.csv(file.path(DEC,"cell_type_proportions_long.csv"), stringsAsFactors=FALSE)
logf("loaded: %d rows, %d samples", nrow(d), length(unique(d$srr)))

## --- FRONTAL ONLY: GSE80655 DLPFC + GSE78936 BA11 ---
d <- d[d$dataset %in% c("GSE80655_DLPFC","GSE78936"),]

## GSE78936 -> BA11 only (drop the duplicated BA24/BA9 donors)
d <- d[!(d$dataset=="GSE78936" & d$brain_region!="Brodmann's Area 11"),]
logf("after frontal + BA11 filter: %d samples", length(unique(d$srr)))

## GSE80655 -> collapse split runs to BioSample (average proportions)
rt <- read.csv(RT, check.names=FALSE); bio <- setNames(rt$BioSample, rt$Run)
d$donor <- ifelse(d$dataset=="GSE80655_DLPFC" & d$srr %in% names(bio), bio[d$srr], d$srr)

## average proportion per donor x signature x cell_type (collapses the split-run pairs)
agg <- aggregate(proportion ~ donor + dataset + signature + cell_type + diagnosis,
                 data=d, FUN=mean)
logf("after donor collapse: %d unique donors", length(unique(agg$donor)))
print(table(unique(agg[,c("donor","dataset","diagnosis")])$dataset,
            unique(agg[,c("donor","dataset","diagnosis")])$diagnosis))

## --- wide matrix for integration: donors x (signature_celltype) ---
agg$feature <- paste(agg$signature, agg$cell_type, sep="_")
wide <- reshape(agg[,c("donor","diagnosis","dataset","feature","proportion")],
                idvar=c("donor","diagnosis","dataset"), timevar="feature",
                direction="wide")
colnames(wide) <- sub("^proportion\\.","",colnames(wide))
logf("\nwide matrix: %d donors x %d features", nrow(wide), ncol(wide)-3)
write.csv(wide, file.path(OUT,"SCZ_composition_donorcorrect_wide.csv"), row.names=FALSE)
write.csv(agg,  file.path(OUT,"SCZ_composition_donorcorrect_long.csv"), row.names=FALSE)

## quick sanity: the key cell types by diagnosis (CA signature)
for (ct in c("Inhibitory","Excitatory","Astrocytes")) {
  col <- paste0("CA_",ct)
  if (col %in% names(wide)) {
    m <- tapply(wide[[col]], wide$diagnosis, median, na.rm=TRUE)
    p <- wilcox.test(wide[[col]] ~ wide$diagnosis)$p.value
    logf("  CA %-12s Control=%.3f SCZ=%.3f  p=%.3f", ct, m["Control"], m["SCZ"], p)
  }
}
logf("\nOutputs -> %s", OUT)
logf("  SCZ_composition_donorcorrect_wide.csv  <- USE THIS FOR INTEGRATION")
