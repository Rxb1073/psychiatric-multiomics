suppressMessages(library(vegan))
D <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final"
W <- file.path(D, "_fillwork")
me <- read.delim(file.path(D,"Microbiome_MDD/sheets/phylo_metadata.tsv"), check.names=FALSE)
me <- me[me[[1]] != "#q2:types", ]; rownames(me) <- me[[1]]
gcol <- grep("^disease$|^diagnosis$", names(me), value=TRUE)[1]
run <- function(dir, label) {
  f <- file.path(dir, "distance-matrix.tsv")
  if (!file.exists(f)) { cat("  ", label, ": not exported\n"); return(invisible()) }
  dm <- as.matrix(read.delim(f, row.names=1, check.names=FALSE))
  i <- intersect(rownames(dm), rownames(me)); dm <- dm[i,i]
  g <- factor(me[i, gcol])
  bd <- betadisper(as.dist(dm), g); pd <- permutest(bd, permutations=999)
  pa <- adonis2(as.dist(dm) ~ g, permutations=999)
  cat(sprintf("  %-22s n=%d  PERMANOVA R2=%.4f p=%.3f | PERMDISP p=%.3f\n", label, length(i), pa$R2[1], pa$`Pr(>F)`[1], pd$tab$`Pr(>F)`[1]))
}
run(file.path(W,"dmU"), "unweighted UniFrac")
run(file.path(W,"dmW"), "weighted UniFrac")
run(file.path(W,"dmJ"), "Jaccard")
