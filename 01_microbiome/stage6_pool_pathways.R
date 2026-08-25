setwd("/rds/projects/e/elhamsak-dissertation-autism/Microbiome/pooled")
studies <- c(PRJNA453621="../PRJNA453621", PRJNA529598="../PRJNA529598",
             PRJNA769228="../PRJNA769228", PRJNA1022510="../PRJNA1022510")

tabs <- list()
for(st in names(studies)){
  f <- file.path(studies[st], "picrust2_out/pathways_out/path_abun_unstrat.tsv.gz")
  df <- read.delim(gzfile(f), row.names=1, check.names=FALSE)
  cat(sprintf("%s: %d samples, %d pathways\n", st, ncol(df), nrow(df)))
  tabs[[st]] <- df
}

# union of all pathway IDs
all_paths <- sort(unique(unlist(lapply(tabs, rownames))))
# build merged matrix (pathways x samples), missing = 0
merged <- do.call(cbind, lapply(tabs, function(df){
  m <- matrix(0, nrow=length(all_paths), ncol=ncol(df),
              dimnames=list(all_paths, colnames(df)))
  m[rownames(df), ] <- as.matrix(df)
  m
}))
cat(sprintf("MERGED: %d samples, %d pathways\n", ncol(merged), nrow(merged)))

# align to pooled metadata
meta <- read.delim("metadata_pooled.tsv", check.names=FALSE)
meta <- meta[meta[[1]] != "#q2:types", ]
rownames(meta) <- meta[[1]]
common <- intersect(colnames(merged), rownames(meta))
cat(sprintf("samples matching metadata: %d of %d\n", length(common), ncol(merged)))

merged <- merged[, common]
write.table(merged, "pooled_pathways.tsv", sep="\t", quote=FALSE, col.names=NA)
write.table(meta[common, ], "pooled_pathways_meta.tsv", sep="\t", quote=FALSE, row.names=FALSE)
cat("wrote pooled_pathways.tsv + pooled_pathways_meta.tsv\n")
