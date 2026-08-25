suppressMessages(library(vegan))
MG <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric2/MDD/microbiome/merged_genus"

raw <- read.delim(file.path(MG,"merged_genus_raw.tsv"), skip=1, row.names=1, check.names=FALSE)
mat <- t(as.matrix(raw))                      # samples x genera
mode(mat) <- "numeric"
meta <- read.delim(file.path(MG,"merged_genus_metadata.tsv"), check.names=FALSE)
rownames(meta) <- meta[[1]]
k <- intersect(rownames(mat), rownames(meta))
mat <- mat[k, ]; meta <- meta[k, ]
cat("samples:", nrow(mat), " genera:", ncol(mat), "\n")

d <- rowSums(mat)
cat(sprintf("depth: min=%.0f  median=%.0f  max=%.0f\n", min(d), median(d), max(d)))
DEPTH <- min(d)
cat("rarefying all samples to", DEPTH, "(the minimum -> no samples lost)\n\n")

set.seed(42); NIT <- 100
O <- S <- P <- matrix(0, nrow(mat), NIT)
for (i in 1:NIT) {
  r <- rrarefy(mat, DEPTH)
  O[,i] <- specnumber(r); S[,i] <- diversity(r,"shannon"); P[,i] <- diversity(r,"simpson")
}
alpha <- data.frame(sample=rownames(mat),
                    observed=rowMeans(O), shannon=rowMeans(S), simpson=rowMeans(P),
                    disease=factor(meta$disease, levels=c("Control","MDD")),
                    study=factor(meta$study))
write.table(alpha, file.path(MG,"alpha_diversity_rarefied.tsv"), sep="\t",
            quote=FALSE, row.names=FALSE)

cat("===== RAREFIED alpha: disease effect controlling for study =====\n")
cat("(", NIT, "rarefaction iterations, averaged )\n")
for (m in c("observed","shannon","simpson")) {
  cat("\n---", m, "---\n")
  f <- as.formula(paste(m, "~ study + disease"))
  print(summary(aov(f, data=alpha))[[1]][, c("F value","Pr(>F)")])
  for (st in levels(alpha$study)) {
    sub <- alpha[alpha$study==st, ]
    w <- suppressWarnings(wilcox.test(sub[[m]] ~ sub$disease))
    cat(sprintf("   %s within %s: W p = %.4f   (Ctrl %.2f vs MDD %.2f)\n",
        m, st, w$p.value,
        median(sub[[m]][sub$disease=="Control"]), median(sub[[m]][sub$disease=="MDD"])))
  }
}
cat("\n===== COMPARISON TO THE NON-RAREFIED RESULT =====\n")
cat("  non-rarefied:  observed p=0.050  shannon p=0.012  simpson p=0.037\n")
cat("  if the rarefied p-values are similar, depth is definitively not the driver.\n")
