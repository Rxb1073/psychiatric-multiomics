D   <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final"
OUT <- file.path(D,"INTEGRATION_inputs","07_chd8_axis")
suppressMessages(library(GenomicRanges))

gt <- read.delim("/rds/projects/e/elhamsak-dissertation-autism/RNA_seq/reference/genome/Homo_sapiens.GRCh38.109.gtf",
                 header=FALSE, comment.char="#", quote="",
                 col.names=c("chr","src","feature","start","end","score","strand","frame","attr"),
                 colClasses=c("character","NULL","character","integer","integer","NULL","character","NULL","character"))
gt <- gt[gt$feature=="gene",]
gt$ens <- sub('.*gene_id "([^"]+)".*',"\\1",gt$attr)
gt$bt  <- sub('.*gene_biotype "([^"]+)".*',"\\1",gt$attr)
gt <- gt[gt$bt %in% c("protein_coding","lncRNA"),]
gt$tss <- ifelse(gt$strand=="+", gt$start, gt$end)
tss <- GRanges(paste0("chr",gt$chr), IRanges(gt$tss,gt$tss), ens=gt$ens)

bed <- read.delim("/rds/projects/e/elhamsak-dissertation-autism/Epigenomics/epigenomics-CHD8/GSE61492/peaks_corrected/CHD8_specific_consensus.bed", header=FALSE)[,1:3]
pk  <- GRanges(bed[,1], IRanges(bed[,2]+1, bed[,3]))
nr  <- distanceToNearest(pk, tss)
ok  <- !is.na(mcols(nr)$distance) & mcols(nr)$distance <= 2000

df <- data.frame(peak=queryHits(nr)[ok], ens=tss$ens[subjectHits(nr)[ok]],
                 dist=mcols(nr)$distance[ok])
write.csv(df, file.path(OUT,"CHD8_peak_to_gene_prededup.csv"), row.names=FALSE)
cat(sprintf("%d of %d sites within 2 kb (%.1f%%) -> %d unique genes\n",
            nrow(df), length(pk), 100*nrow(df)/length(pk), length(unique(df$ens))))
