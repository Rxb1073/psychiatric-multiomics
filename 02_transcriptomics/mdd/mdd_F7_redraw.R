## F7 mode-of-death panel, standalone redraw.
## Horizontal bars; annotations in a right-hand column so they cannot reach
## the axis labels; full GEO accessions rather than truncated numbers.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
if (isTRUE(capabilities("cairo"))) options(bitmapType = "cairo")
logf <- function(...) cat(sprintf(...), "\n")

M   <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final/Transcriptomics_MDD"
TAB <- file.path(M, "tables"); OUT <- file.path(M, "figures")
MDDC <- "#C0392B"; C1 <- "#2E86C1"

sd_ <- read.delim(file.path(TAB, "16_suicide_decomposition.tsv"), check.names = FALSE)
i <- sd_[sd_$set == "IEG", ]
logf("%d contrasts: %s", nrow(i), paste(unique(i$test), collapse = " "))

pretty <- c(
  T1_modeofdeath_in_controls = "suicide vs natural death,\nwithin controls",
  T2_MDD_vs_suicide_controls = "MDD vs suicide controls",
  MDD_vs_natural_controls    = "MDD vs natural-death controls",
  depression_no_suicide      = "depression without suicide",
  suicide_and_MDD            = "depression with suicide",
  suicide_within_MDD         = "suicide vs non-suicide,\nwithin depression")
i$lab <- ifelse(i$contrast %in% names(pretty), pretty[i$contrast], i$contrast)
i$lab <- paste0(i$lab, "\n", i$test)

o <- order(i$ieg_median_logFC); v <- i$ieg_median_logFC[o]
cols <- ifelse(v > 0, C1, MDDC); rng <- range(c(v, 0))
ann_x <- rng[2] + diff(rng) * 0.10
xlim_ <- c(rng[1] - diff(rng) * 0.04, rng[2] + diff(rng) * 0.85)

draw <- function() {
  par(mar = c(4.6, 16, 3.2, 1.5), mgp = c(2.7, .7, 0))
  bp <- barplot(v, horiz = TRUE, names.arg = i$lab[o], las = 1, cex.names = .72,
                col = cols, border = NA,
                xlab = expression(median ~ log[2] ~ FC ~ (IEG ~ set)),
                main = "Immediate-early module by contrast", xlim = xlim_)
  abline(v = 0)
  text(ann_x, bp, sprintf("p = %.3f    %d of 16 reduced", i$p_corr[o], i$ieg_negative[o]),
       pos = 4, cex = .68, xpd = NA)
  mtext("p-values from CAMERA with inter-gene correlation estimated from the data",
        side = 1, line = 3.3, cex = .66, col = "grey35")
}
pdf(file.path(OUT, "F7_mode_of_death.pdf"), width = 10, height = 5.6); draw(); dev.off()
png(file.path(OUT, "F7_mode_of_death.png"), width = 10, height = 5.6, units = "in", res = 150); draw(); dev.off()
logf("wrote F7_mode_of_death .pdf/.png to %s", OUT)
