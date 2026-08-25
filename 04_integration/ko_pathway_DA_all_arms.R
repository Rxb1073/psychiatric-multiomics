## KEGG-orthologue pathway differential abundance, three arms.
##
## Same PICRUSt2 runs, same MaAsLin2 models, same normalisation and prevalence
## filter as each arm's MetaCyc analysis -- only the input table changes, from
## pathways_out/path_abun_unstrat to KO_metagenome_out/pred_metagenome_unstrat.
## KOs are summed within KEGG pathway BEFORE testing, so the unit of inference
## is the pathway, matching the MetaCyc analysis and the host KEGG results.
##
## Models, unchanged from the MetaCyc scripts:
##   ASD  ~ study + diagnosis   (4 projects, reference PRJNA453621 / Control)
##   SCZ  ~ diagnosis           (single cohort PRJNA880407, reference HC)
##   MDD  ~ study + disease     (2 projects, reference Control)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(library(Maaslin2))

BASE_A <- "/rds/projects/e/elhamsak-dissertation-autism/Microbiome"
BASE_S <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA880407/picrust2"
BASE_M <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome"
HARM   <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/harmonization"
AX     <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final/INTEGRATION_inputs/06_kegg_axis"
dir.create(AX, recursive=TRUE, showWarnings=FALSE)

MAP <- read.csv(file.path(AX,"ko_pathway_map.csv"), stringsAsFactors=FALSE)
cat(sprintf("KO map: %d links, %d KOs, %d pathways\n\n",
            nrow(MAP), length(unique(MAP$KO)), length(unique(MAP$pathway))))

read_ko <- function(f){
  stopifnot(file.exists(f))
  d <- read.delim(gzfile(f), row.names=1, check.names=FALSE)
  d[!is.na(rownames(d)) & rownames(d)!="", , drop=FALSE]
}

## sum KO abundances within pathway; a KO in several pathways contributes to each
collapse_to_pathway <- function(ko_mat, label){
  m <- MAP[MAP$KO %in% rownames(ko_mat), c("KO","pathway")]
  hit <- length(unique(m$KO))
  cat(sprintf("  %s: %d KOs, %d (%.1f%%) map to >=1 pathway\n",
              label, nrow(ko_mat), hit, 100*hit/nrow(ko_mat)))
  agg <- rowsum(as.matrix(ko_mat)[m$KO, , drop=FALSE], group=m$pathway)
  cat(sprintf("  %s: collapsed to %d pathways x %d samples\n", label, nrow(agg), ncol(agg)))
  agg
}

run_maaslin <- function(feat_pathways, meta, fixed, ref, arm, contrast_value=NULL){
  out <- file.path(AX, sprintf("maaslin2_KO_%s", arm))
  fit <- Maaslin2(input_data      = as.data.frame(t(feat_pathways)),
                  input_metadata  = meta,
                  output          = out,
                  fixed_effects   = fixed,
                  reference       = ref,
                  normalization   = "TSS", transform = "LOG",
                  min_prevalence  = 0.10, standardize = FALSE,
                  plot_heatmap    = FALSE, plot_scatter = FALSE)
  dxvar <- tail(fixed, 1)
  res <- fit$results[fit$results$metadata == dxvar, ]
  if(!is.null(contrast_value)) res <- res[res$value == contrast_value, ]
  res <- res[order(res$qval), ]
  res$suffix <- sub("^ko", "", res$feature)
  res <- merge(res, unique(MAP[,c("pathway","pathway_name","is_global")]),
               by.x="feature", by.y="pathway", all.x=TRUE)
  res <- res[order(res$qval), ]
  f <- file.path(AX, sprintf("KO_pathway_DA_%s.tsv", arm))
  write.table(res, f, sep="\t", quote=FALSE, row.names=FALSE)
  n <- sum(res$qval < 0.05, na.rm=TRUE)
  cat(sprintf("\n  === %s: %d of %d pathways at q<0.05 (%d down, %d up) ===\n",
              arm, n, nrow(res),
              sum(res$qval<0.05 & res$coef<0, na.rm=TRUE),
              sum(res$qval<0.05 & res$coef>0, na.rm=TRUE)))
  if(n) print(head(res[res$qval<0.05, c("feature","pathway_name","coef","qval","is_global")], 12), row.names=FALSE)
  cat(sprintf("  wrote %s\n", f))
  invisible(res)
}

## ---------------------------------------------------------------- ASD -------
cat("\n================ ASD ================\n")
studies_a <- c("PRJNA453621","PRJNA529598","PRJNA769228","PRJNA1022510")
tabs <- lapply(studies_a, function(st)
  read_ko(file.path(BASE_A, st, "picrust2_out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz")))
names(tabs) <- studies_a
for(st in studies_a) cat(sprintf("  %s: %d KOs x %d samples\n", st, nrow(tabs[[st]]), ncol(tabs[[st]])))

all_ko <- sort(unique(unlist(lapply(tabs, rownames))))
merged_a <- do.call(cbind, lapply(tabs, function(df){
  m <- matrix(0, nrow=length(all_ko), ncol=ncol(df), dimnames=list(all_ko, colnames(df)))
  m[rownames(df), ] <- as.matrix(df); m }))
cat(sprintf("  merged: %d KOs x %d samples\n", nrow(merged_a), ncol(merged_a)))

meta_a <- read.delim(file.path(BASE_A,"pooled/metadata_pooled.tsv"), check.names=FALSE)
meta_a <- meta_a[meta_a[[1]] != "#q2:types", ]
rownames(meta_a) <- meta_a[[1]]
common <- intersect(colnames(merged_a), rownames(meta_a))
cat(sprintf("  samples matching pooled metadata: %d of %d\n", length(common), ncol(merged_a)))
merged_a <- merged_a[, common]; meta_a <- meta_a[common, , drop=FALSE]
meta_a$diagnosis <- factor(meta_a$diagnosis, levels=c("Control","ASD"))
meta_a$study     <- factor(meta_a$study)
print(table(meta_a$study, meta_a$diagnosis))
pw_a <- collapse_to_pathway(merged_a, "ASD")
res_a <- run_maaslin(pw_a, meta_a, c("study","diagnosis"),
                     "study,PRJNA453621;diagnosis,Control", "ASD")

## ---------------------------------------------------------------- SCZ -------
cat("\n================ SCZ ================\n")
ko_s <- read_ko(file.path(BASE_S,"picrust2_output/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"))
cat(sprintf("  PRJNA880407: %d KOs x %d samples\n", nrow(ko_s), ncol(ko_s)))
meta_s <- read.delim(file.path(dirname(BASE_S),"metadata/metadata_qiime2.tsv"), check.names=FALSE)
rownames(meta_s) <- meta_s$`sample-id`
common <- intersect(colnames(ko_s), rownames(meta_s))
cat(sprintf("  samples matching metadata: %d of %d\n", length(common), ncol(ko_s)))
ko_s <- ko_s[, common]; meta_s <- meta_s[common, , drop=FALSE]
meta_s$diagnosis <- factor(meta_s$diagnosis, levels=c("HC","SZ"))
print(table(meta_s$diagnosis))
pw_s <- collapse_to_pathway(ko_s, "SCZ")
res_s <- run_maaslin(pw_s, meta_s, c("diagnosis"), c("diagnosis,HC"), "SCZ", contrast_value="SZ")

## ---------------------------------------------------------------- MDD -------
cat("\n================ MDD ================\n")
read_harm <- function(name){
  h <- read.delim(file.path(HARM, paste0(name,"_harmonized.tsv")), check.names=FALSE)
  h <- h[!grepl("^#", h[[1]]), ]
  setNames(h$disease, h$`sample-id`) }
p1 <- read_ko(file.path(BASE_M,"picrust2/PRJNA591924/picrust2_output/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"))
p2 <- read_ko(file.path(BASE_M,"picrust2/PRJNA776170/picrust2_output/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz"))
cat(sprintf("  PRJNA591924: %d KOs x %d samples\n  PRJNA776170: %d KOs x %d samples (pre-filter)\n",
            nrow(p1), ncol(p1), nrow(p2), ncol(p2)))
lab1 <- read_harm("PRJNA591924"); lab2 <- read_harm("PRJNA776170")
k1 <- intersect(colnames(p1), names(lab1)); p1 <- p1[, k1]
k2 <- intersect(colnames(p2), names(lab2)); p2 <- p2[, k2]
cat(sprintf("  after harmonisation filter: %d and %d samples\n", ncol(p1), ncol(p2)))
common_ko <- intersect(rownames(p1), rownames(p2))
cat(sprintf("  KOs common to both projects: %d\n", length(common_ko)))
merged_m <- cbind(p1[common_ko, ], p2[common_ko, ])
meta_m <- data.frame(sample  = colnames(merged_m),
                     disease = c(lab1[k1], lab2[k2]),
                     study   = c(rep("PRJNA591924", ncol(p1)), rep("PRJNA776170", ncol(p2))),
                     row.names = colnames(merged_m), stringsAsFactors=FALSE)
meta_m$disease <- factor(meta_m$disease, levels=c("Control","MDD"))
print(table(meta_m$study, meta_m$disease))
pw_m <- collapse_to_pathway(merged_m, "MDD")
res_m <- run_maaslin(pw_m, meta_m, c("study","disease"), c("disease,Control"), "MDD")

cat("\n================ summary ================\n")
for(a in c("ASD","SCZ","MDD")){
  r <- get(paste0("res_", tolower(a)))
  cat(sprintf("  %-4s %4d pathways tested | %3d at q<0.05 | %3d excluding global/overview maps\n",
              a, nrow(r), sum(r$qval<0.05, na.rm=TRUE),
              sum(r$qval<0.05 & !r$is_global, na.rm=TRUE)))
}
cat(sprintf("\noutputs in %s\n", AX))
quit(save="no")
