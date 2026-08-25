D   <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final"
CP  <- file.path(D, "Transcriptomics_SCZ/integration/SCZ_composition_donorcorrect_wide.csv")
OUT <- file.path(D, "Transcriptomics_SCZ/tables")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cp <- read.csv(CP, check.names = FALSE, stringsAsFactors = FALSE)
cat(sprintf("donors: %d | %s\n", nrow(cp),
            paste(names(table(cp$diagnosis)), table(cp$diagnosis), sep="=", collapse=" ")))
cat(sprintf("datasets: %s\n\n", paste(names(table(cp$dataset)), table(cp$dataset), sep="=", collapse=" ")))

cols <- grep("^(CA|MB)_", names(cp), value = TRUE)
cp$CA_Neuronal_total <- cp$CA_Excitatory + cp$CA_Inhibitory
cols <- c(cols, "CA_Neuronal_total")

res <- do.call(rbind, lapply(cols, function(v) {
  x  <- cp[[v]]
  ct <- x[cp$diagnosis == "Control"]; sz <- x[cp$diagnosis == "SCZ"]
  m  <- summary(lm(x ~ dataset + diagnosis, cp))$coefficients
  r  <- m[grep("^diagnosis", rownames(m)), ]
  w  <- suppressWarnings(wilcox.test(sz, ct))
  data.frame(signature = ifelse(grepl("^CA_", v), "Cell Atlas", "MultiBrain"),
             cell_type = sub("^(CA|MB)_", "", v),
             mean_control = mean(ct), mean_SCZ = mean(sz),
             beta_lm = unname(r[1]), p_lm = unname(r[4]),
             p_wilcoxon = unname(w$p.value), stringsAsFactors = FALSE)
}))
res$q_lm <- ave(res$p_lm, res$signature, FUN = function(p) p.adjust(p, "BH"))
res <- res[order(res$signature, res$p_lm), ]

cat("=== all cell types, linear model adjusting for dataset ===\n")
print(format(res[, c("signature","cell_type","mean_control","mean_SCZ","beta_lm","p_lm","q_lm","p_wilcoxon")],
             digits = 3), row.names = FALSE)
f <- file.path(OUT, "SCZ_composition_tests_alltypes.csv")
write.csv(res, f, row.names = FALSE)
cat(sprintf("\nwritten: %s\n", f))
