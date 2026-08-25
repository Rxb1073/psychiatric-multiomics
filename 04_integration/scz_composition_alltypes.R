## SCZ cell-type composition: one documented model, all cell types, saved output.
##
## WHY. The chapter reports reduced inhibitory (P = 0.008) and total neuronal
## (P = 0.040) fractions, but the analysis that produced those values left no
## output table on disk. Microglial values are therefore not recoverable from the
## same model, and three different microglial p-values exist from three different
## tests. This script recomputes every cell type under one stated model and writes
## the table, so all values in the chapter come from a single documented analysis.
##
## The model matches the one described in the Methods: cell-type proportion as the
## response, dataset as a covariate, diagnosis as the term of interest, on the
## donor-corrected 76-donor table. Both reference signatures are run.
##
## A Wilcoxon test is reported alongside for every cell type, because an earlier
## unrecorded analysis appears to have used one; this makes the difference between
## the two visible rather than leaving it to be discovered.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

D   <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final"
CP  <- file.path(D, "Transcriptomics_SCZ/integration/SCZ_composition_donorcorrect_wide.csv")
OUT <- file.path(D, "Transcriptomics_SCZ/tables")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cp <- read.csv(CP, check.names = FALSE, stringsAsFactors = FALSE)
cat(sprintf("donors: %d | %s\n", nrow(cp),
            paste(names(table(cp$diagnosis)), table(cp$diagnosis), sep = "=", collapse = " ")))
cat(sprintf("datasets: %s\n\n", paste(names(table(cp$dataset)), table(cp$dataset),
                                      sep = "=", collapse = " ")))

## every proportion column, both signatures
cols <- grep("^(CA|MB)_", names(cp), value = TRUE)

## derived totals, as used in the chapter
cp$CA_Neuronal_total <- cp$CA_Excitatory + cp$CA_Inhibitory
cols <- c(cols, "CA_Neuronal_total")

res <- do.call(rbind, lapply(cols, function(v) {
  x  <- cp[[v]]
  ct <- x[cp$diagnosis == "Control"]
  sz <- x[cp$diagnosis == "SCZ"]
  ## primary: linear model adjusting for dataset
  m  <- summary(lm(x ~ dataset + diagnosis, cp))$coefficients
  r  <- m[grep("^diagnosis", rownames(m)), ]
  ## comparator: unadjusted Wilcoxon
  w  <- suppressWarnings(wilcox.test(sz, ct))
  data.frame(
    signature   = ifelse(grepl("^CA_", v), "Cell Atlas", "MultiBrain"),
    cell_type   = sub("^(CA|MB)_", "", v),
    mean_control= mean(ct), mean_SCZ = mean(sz),
    beta_lm     = unname(r[1]),
    p_lm        = unname(r[4]),
    p_wilcoxon  = unname(w$p.value),
    stringsAsFactors = FALSE)
}))

## correction within each signature, across cell types
res$q_lm <- ave(res$p_lm, res$signature, FUN = function(p) p.adjust(p, "BH"))
res <- res[order(res$signature, res$p_lm), ]

cat("=== all cell types, linear model adjusting for dataset ===\n")
print(format(res[, c("signature","cell_type","mean_control","mean_SCZ","beta_lm","p_lm","q_lm","p_wilcoxon")],
             digits = 3), row.names = FALSE)

f <- file.path(OUT, "SCZ_composition_tests_alltypes.csv")
write.csv(res, f, row.names = FALSE)
cat(sprintf("\nwritten: %s\n", f))

cat("\n=== the three values the chapter needs ===\n")
key <- res[res$cell_type %in% c("Inhibitory","Neuronal_total","Neurons","Microglia"), ]
print(format(key[, c("signature","cell_type","p_lm","q_lm","p_wilcoxon")], digits = 3), row.names = FALSE)

cat("\nCompare the inhibitory and neuronal rows against the chapter's P = 0.008 and\n")
cat("0.040. If they match, quote the microglial rows from this table alongside them.\n")
cat("If they do not, the chapter's values came from a different model and the two\n")
cat("sets must not be mixed in one sentence.\n")
