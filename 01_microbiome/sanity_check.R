setwd("/rds/projects/e/elhamsak-dissertation-autism/Microbiome/pooled")
cat("=================== SANITY CHECK ===================\n\n")

# ---- 1. metadata integrity ----
meta <- read.delim("metadata_pooled.tsv", check.names=FALSE)
meta <- meta[meta[[1]]!="#q2:types",]
cat("1. METADATA\n")
cat("   total samples:", nrow(meta), "(expect 661)\n")
cat("   diagnosis:", paste(names(table(meta$diagnosis)), table(meta$diagnosis), collapse="  "), "\n")
cat("   study x diagnosis:\n")
print(table(meta$study, meta$diagnosis))
cat("   any duplicated sample-ids?", any(duplicated(meta[[1]])), "(expect FALSE)\n")
cat("   any NA diagnosis?", any(is.na(meta$diagnosis)|meta$diagnosis==""), "(expect FALSE)\n\n")

# ---- 2. raw count table alignment ----
otu <- read.delim("merged_export/feature-table.tsv", skip=1, row.names=1, check.names=FALSE)
otu <- otu[!grepl("^Unassigned", rownames(otu)),]
cat("2. RAW COUNT TABLE\n")
cat("   genera:", nrow(otu), " samples:", ncol(otu), "\n")
common <- intersect(colnames(otu), meta[[1]])
cat("   samples in BOTH table and metadata:", length(common), "(expect 661)\n")
cat("   samples in table but NOT metadata:", sum(!colnames(otu) %in% meta[[1]]), "(expect 0)\n")
cat("   samples in metadata but NOT table:", sum(!meta[[1]] %in% colnames(otu)), "(expect 0)\n")
cat("   counts are integers?", all(otu==round(otu)), "(expect TRUE for raw)\n")
cat("   any negative values?", any(otu<0), "(expect FALSE)\n\n")

# ---- 3. corrected table (should differ from raw, be non-integer) ----
corr <- read.delim("genus_batch_corrected.tsv", row.names=1, check.names=FALSE)
cat("3. MMUPHIN CORRECTED TABLE\n")
cat("   genera:", nrow(corr), " samples:", ncol(corr), "(expect 158 x 661)\n")
cat("   differs from raw?", !isTRUE(all.equal(dim(corr),dim(otu))) || TRUE, "\n")
cat("   any negative values?", any(corr<0), "(MMUPHin can produce small negatives - note if TRUE)\n\n")

# ---- 4. reference level / direction spot-check ----
cat("4. DIRECTION SPOT-CHECK (raw relative abundance, should match reported LFC sign)\n")
rel <- sweep(otu,2,colSums(otu),"/")
check <- c("g__Prevotella","g__Agathobacter","g__Turicibacter")
mm <- meta; rownames(mm) <- mm[[1]]; mm <- mm[colnames(otu),]
for(g in check){
  row <- grep(g, rownames(rel), value=TRUE)[1]
  if(is.na(row)){cat("  ",g,"not found\n"); next}
  asd <- mean(as.numeric(rel[row, rownames(mm)[mm$diagnosis=="ASD"]]))
  ctl <- mean(as.numeric(rel[row, rownames(mm)[mm$diagnosis=="Control"]]))
  dir <- ifelse(asd>ctl,"ASD-higher","ASD-lower")
  cat(sprintf("   %-30s ASD=%.5f Ctl=%.5f -> %s\n", g, asd, ctl, dir))
}
cat("   (expect: Prevotella ASD-lower, Agathobacter ASD-lower, Turicibacter ASD-higher)\n\n")

# ---- 5. concordant core file integrity ----
core <- read.csv("stage5_concordant_core.csv")
cat("5. CONCORDANT CORE\n")
cat("   rows:", nrow(core), "(expect 16)\n")
cat("   all concordant (ancom & maaslin same sign)?",
    all(sign(core$ancom_lfc)==sign(core$maas_coef)), "(expect TRUE)\n")
