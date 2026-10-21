suppressMessages({library(ggplot2); library(reshape2)})
D <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final"
S <- "/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA880407"
W <- file.path(D, "_figwork")
th <- theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = "grey92", colour = NA))
save2 <- function(p, dir, name, w, h) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  cat("  wrote", name, "(pdf + png)\n")
}
# last taxonomic token, preserving distinctness (torques_group stays separate)
lt <- function(x) {
  x <- as.character(x)
  x <- gsub('"', "", x)
  x <- gsub("\\.(?=[dpcofgs]__)", ";", x, perl = TRUE)   # ASD writes lineages with dots
  x <- gsub(";__", "", x)
  x <- sub(";+_*$", "", x)
  x <- sub(".*;", "", x)
  x <- sub("^g__", "", x)
  x <- gsub("^\\.|\\.$", "", x)
  x
}
pick <- function(df, pats, what, file) {
  for (p in pats) { h <- grep(p, names(df), ignore.case = TRUE, value = TRUE)
                    if (length(h)) return(h[1]) }
  stop("Could not find ", what, " in ", file, "\n  columns present: ",
       paste(names(df), collapse = ", "))
}

# ============================================================ CROSS-ARM
cat("\n=== CROSSARM ===\n")
grab <- function(file, arm) {
  if (!file.exists(file)) { cat("  MISSING:", file, "\n"); return(NULL) }
  d <- read.delim(file, sep = if (grepl("csv$", file)) "," else "\t", check.names = FALSE)
  gcol <- pick(d, c("^genus$","^taxon$","^feature$","^key$","^id$"), "genus column", basename(file))
  lcol <- pick(d, c("log2FC","ancom_lfc","ancombc_log2FC","^lfc","^coef$"), "effect size", basename(file))
  qcol <- pick(d, c("^q_value$","^qval","ancombc_q","maas_q","^q$","^padj"), "q-value", basename(file))
  cat("  ", arm, ":", basename(file), "->", gcol, "/", lcol, "/", qcol,
      "  rows:", nrow(d), "\n")
  scol <- grep("^significant$", names(d), value = TRUE)
  if (length(scol)) d <- d[which(as.logical(d[[scol]])), ]
  d <- d[is.finite(d[[qcol]]) & d[[qcol]] < 0.05, ]
  data.frame(genus = lt(d[[gcol]]), lfc = as.numeric(d[[lcol]]),
             q = as.numeric(d[[qcol]]), arm = arm, stringsAsFactors = FALSE)
}
L <- list(
  grab(file.path(D, "Microbiome_ASD/tables/stage5_concordant_core.csv"), "ASD"),
  grab(file.path(D, "Microbiome_SCZ/tables/robust_genera_both_methods_v2.tsv"), "SCZ"),
  grab(file.path(D, "Microbiome_MDD/tables/robust_genera_both_methods.tsv"), "MDD"))
X <- do.call(rbind, L[!sapply(L, is.null)])
X <- X[!grepl("^uncultured$|^Bacteria$|^$", X$genus), ]
X <- X[!duplicated(X[, c("genus","arm")]), ]
n_arms <- table(X$genus)
shared <- names(n_arms)[n_arms >= 2]
cat("  genera in >=2 arms:", length(shared), "\n")
if (length(shared)) cat("  ", paste(sort(shared), collapse = ", "), "\n")
SCFA <- c("Faecalibacterium","Fusicatenibacter","Anaerostipes","Blautia","Coprococcus",
          "Roseburia","Agathobacter","Subdoligranulum","Butyricicoccus","Marvinbryantia")
keep <- union(shared, intersect(X$genus, SCFA))
Y <- X[X$genus %in% keep, ]
Y$genus <- factor(Y$genus, levels = names(sort(tapply(Y$lfc, Y$genus, mean))))
Y$arm <- factor(Y$arm, levels = c("ASD","SCZ","MDD"))
Y$lab <- ifelse(Y$genus %in% SCFA, "SCFA producer", "other")
p <- ggplot(Y, aes(arm, genus)) +
  geom_point(aes(fill = lfc), size = 4.2, shape = 21, colour = "grey25", stroke = .3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = expression(log[2]*"FC (case vs control)")) +
  facet_grid(lab ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL) + th
save2(p, file.path(D, "CROSS_ARM_figures"), "FigX_crossarm_shared_genera",
      6.5, max(4, 0.22 * length(unique(Y$genus)) + 2))
write.csv(Y[order(Y$genus, Y$arm), ], file.path(D, "CROSS_ARM_figures", "FigX_source_data.csv"),
          row.names = FALSE)

# ============================================================ SCZ Fig1 alpha
cat("\n=== SCZ Fig1 alpha ===\n")
meta <- read.delim(file.path(S, "metadata/metadata_qiime2.tsv"), check.names = FALSE)
meta <- meta[meta[[1]] != "#q2:types", ]
dcol <- pick(meta, c("^diagnosis$","^group$"), "diagnosis", "metadata_qiime2.tsv")
rownames(meta) <- meta[[1]]
nice <- c(shannon = "Shannon", observed_features = "Observed features",
          evenness = "Pielou evenness", faith_pd = "Faith's PD")
A <- do.call(rbind, lapply(names(nice), function(v) {
  f <- file.path(W, paste0("a_", v), "alpha-diversity.tsv")
  if (!file.exists(f)) { cat("  MISSING", v, "\n"); return(NULL) }
  d <- read.delim(f, row.names = 1); ids <- intersect(rownames(d), rownames(meta))
  cat("  ", v, ":", length(ids), "samples\n")
  data.frame(id = ids, value = d[ids, 1], metric = nice[[v]],
             diagnosis = meta[ids, dcol], stringsAsFactors = FALSE)
}))
if (!is.null(A)) {
  pv <- sapply(split(A, A$metric), function(s)
          tryCatch(kruskal.test(value ~ factor(diagnosis), data = s)$p.value, error = function(e) NA))
  A$metric <- factor(A$metric, levels = nice)
  lab <- data.frame(metric = factor(names(pv), levels = nice),
                    txt = paste0("Kruskal-Wallis p = ", signif(pv, 3)))
  p <- ggplot(A, aes(diagnosis, value, fill = diagnosis)) +
    geom_boxplot(outlier.size = .6, width = .6, alpha = .85) +
    geom_jitter(width = .12, size = .5, alpha = .35) +
    facet_wrap(~ metric, scales = "free_y", nrow = 2) +
    geom_text(data = lab, aes(x = 1.5, y = Inf, label = txt), vjust = 1.5,
              size = 3, inherit.aes = FALSE) +
    scale_fill_manual(values = c("#4477AA", "#CC6677")) +
    labs(x = NULL, y = NULL, title = "Alpha diversity is reduced in schizophrenia",
         subtitle = "SEPP tree; rarefied to 7,243 reads; n = 50 per group") +
    th + theme(legend.position = "none")
  save2(p, file.path(D, "Microbiome_SCZ/figures"), "Fig1_alpha_diversity", 7, 5.5)
}

# ============================================================ SCZ Fig2 PCoA
cat("\n=== SCZ Fig2 PCoA ===\n")
dmf <- file.path(W, "uu_dm", "distance-matrix.tsv")
if (file.exists(dmf)) {
  dm <- as.matrix(read.delim(dmf, row.names = 1, check.names = FALSE))
  ids <- intersect(rownames(dm), rownames(meta)); dm <- dm[ids, ids]
  cm <- cmdscale(as.dist(dm), k = 2, eig = TRUE)
  ve <- round(100 * cm$eig[1:2] / sum(cm$eig[cm$eig > 0]), 1)
  P <- data.frame(PC1 = cm$points[,1], PC2 = cm$points[,2], diagnosis = meta[ids, dcol])
  cat("  ", length(ids), "samples; PC1", ve[1], "% PC2", ve[2], "%\n")
  p <- ggplot(P, aes(PC1, PC2, colour = diagnosis, fill = diagnosis)) +
    stat_ellipse(geom = "polygon", alpha = .12, colour = NA, level = .95) +
    geom_point(size = 2, alpha = .85) +
    scale_colour_manual(values = c("#4477AA","#CC6677")) +
    scale_fill_manual(values = c("#4477AA","#CC6677")) +
    labs(x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"),
         title = "Community composition differs by diagnosis",
         subtitle = "Unweighted UniFrac (SEPP). PERMANOVA pseudo-F = 8.27, p = 0.001; PERMDISP p = 0.93") + th
  save2(p, file.path(D, "Microbiome_SCZ/figures"), "Fig2_beta_pcoa_unweighted_unifrac", 6.5, 5)
} else cat("  distance matrix not exported - skipped\n")

# ============================================================ SCZ Fig3 + Fig4
cat("\n=== SCZ Fig3/Fig4 (26 genera, filter fixed) ===\n")
R <- read.delim(file.path(D, "Microbiome_SCZ/tables/robust_genera_both_methods_v2.tsv"),
                check.names = FALSE)
R$g <- lt(R[[pick(R, c("^key$","^genus$"), "genus", "robust_v2")]])
# the fix: drop unresolved AND rank-prefixed keys, or d__Bacteria returns as the largest bar
R <- R[!grepl("uncultured|^Bacteria$|^$|^(d|p|c|o|f)__", R$g), ]
cat("  genera plotted:", nrow(R), "(expect 26)\n")
R$dir <- ifelse(R$log2FC_SZ_vs_HC < 0, "depleted in SZ", "enriched in SZ")
R$scfa <- R$g %in% SCFA
p <- ggplot(R, aes(log2FC_SZ_vs_HC, maaslin_coef)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(aes(colour = dir, shape = scfa), size = 2.6) +
  scale_colour_manual(values = c("#2166AC","#B2182B"), name = NULL) +
  scale_shape_manual(values = c(16, 17), labels = c("other","SCFA producer"), name = NULL) +
  labs(x = "ANCOM-BC2 log2FC (SZ vs HC)", y = "MaAsLin2 coefficient",
       title = paste0("Two-method concordance: all ", nrow(R), " genera agree in direction"),
       subtitle = "ANCOM-BC2 (BH) and MaAsLin2, both q < 0.05") + th
save2(p, file.path(D, "Microbiome_SCZ/figures"), "Fig3_concordance_scatter", 6.5, 5)
R$g <- factor(R$g, levels = R$g[order(R$log2FC_SZ_vs_HC)])
p <- ggplot(R, aes(log2FC_SZ_vs_HC, g, fill = dir)) +
  geom_col(width = .7) +
  geom_vline(xintercept = 0, colour = "grey30") +
  scale_fill_manual(values = c("#2166AC","#B2182B"), name = NULL) +
  labs(x = "ANCOM-BC2 log2FC (SZ vs HC)", y = NULL,
       title = paste0("Genera concordantly altered in schizophrenia (n = ", nrow(R), ")"),
       subtitle = "Both methods q < 0.05; ranked by effect size") +
  th + theme(axis.text.y = element_text(face = ifelse(
        levels(R$g) %in% SCFA, "bold", "plain")))
save2(p, file.path(D, "Microbiome_SCZ/figures"), "Fig4_genus_effectsizes", 6.5, 6.5)

# ============================================================ ASD Fig7 + Fig8
cat("\n=== ASD Fig7 pathway volcano / Fig8 bridge ===\n")
f <- file.path(D, "Microbiome_ASD/tables/stage6_pathway_results.csv")
if (file.exists(f)) {
  P6 <- read.delim(f, sep = ",", check.names = FALSE)
  cc <- pick(P6, c("^coef$"), "coef", "stage6"); qq <- pick(P6, c("^qval","^q_value$"), "q", "stage6")
  P6 <- P6[is.finite(P6[[cc]]) & is.finite(P6[[qq]]), ]
  P6$sig <- ifelse(P6[[qq]] < 0.05, ifelse(P6[[cc]] < 0, "ASD-down", "ASD-up"), "ns")
  cat("  pathways:", nrow(P6), " significant:", sum(P6$sig != "ns"), "\n")
  p <- ggplot(P6, aes(.data[[cc]], -log10(.data[[qq]]), colour = sig)) +
    geom_point(size = 1.5, alpha = .8) +
    geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey40") +
    scale_colour_manual(values = c("ASD-down" = "#2166AC", "ASD-up" = "#B2182B", "ns" = "grey75"),
                        name = NULL) +
    labs(x = "MaAsLin2 coefficient", y = expression(-log[10]*q),
         title = "Predicted pathway differential abundance (PICRUSt2, MetaCyc)",
         subtitle = "Diffuse, low-effect-size spread; weighted NSTI 0.078-0.203, two cohorts above 0.15") + th
  save2(p, file.path(D, "Microbiome_ASD/figures"), "Fig7_pathway_volcano", 6.5, 5)
} else cat("  stage6_pathway_results.csv missing\n")

f <- file.path(D, "Microbiome_ASD/tables/taxonomy_function_bridge_asd.tsv")
if (file.exists(f)) {
  B <- read.delim(f, check.names = FALSE)
  cat("  bridge columns:", paste(names(B), collapse = ", "), "\n")
  gcol <- pick(B, c("^genus$"), "genus", "bridge"); scol <- pick(B, c("share","pct"), "share", "bridge")
  kcol <- grep("^kind$|block", names(B), value = TRUE)
  B[[scol]] <- as.numeric(B[[scol]])
  if (length(kcol)) {
    B <- do.call(rbind, lapply(split(B, B[[kcol]]), function(s)
           head(s[order(-s[[scol]]), ], 12)))
  } else B <- head(B[order(-B[[scol]]), ], 12)
  B$scfa <- lt(B[[gcol]]) %in% SCFA
  B[[gcol]] <- factor(B[[gcol]], levels = unique(B[[gcol]][order(B[[scol]])]))
  p <- ggplot(B, aes(.data[[scol]], .data[[gcol]], fill = scfa)) +
    geom_col(width = .7) +
    scale_fill_manual(values = c("grey70","#2166AC"),
                      labels = c("other","SCFA producer"), name = NULL) +
    labs(x = "% of contribution", y = NULL,
         title = "Taxonomy-function bridge: who builds the depleted pathways",
         subtitle = "Descriptive, not a test - contributions summed across all samples, ASD not compared to control") + th
    if (length(kcol)) p <- p + facet_wrap(as.formula(paste("~", kcol)), scales = "free_y")
  save2(p, file.path(D, "Microbiome_ASD/figures"), "Fig8_taxonomy_function_bridge", 7, 6)
} else cat("  taxonomy_function_bridge_asd.tsv missing\n")

cat("\n=== DONE ===\n")
