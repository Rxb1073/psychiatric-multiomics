## Panel d for the CHD8 figure: dose-response of transcriptomic displacement
## against binding strength, replacing the four-row table.
##
## Plots the observed median log2 fold change at each nested set size against the
## expression-matched null, so the gradient is visible rather than read off rows.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressMessages(library(ggplot2))

D   <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final"
OUT <- file.path(D, "CROSS_ARM_figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## Values as computed by chd8_dose_response_both_arms.R (ASD arm).
## sd_null is the spread of the 10,000 matched null medians, used for the band.
d <- data.frame(
  sites   = c(1000, 2000, 3000, 4351),
  label   = c("strongest\n1,000", "strongest\n2,000", "strongest\n3,000", "all\n4,351"),
  genes   = c(896, 1758, 2621, 3820),
  median  = c(0.02772747, 0.02526009, 0.02182977, 0.01935687),
  null    = c(-0.001880227, -0.003994195, -0.005463142, -0.007336792),
  sd_null = c(0.0057, 0.0040, 0.0033, 0.0026)
)
d$x <- factor(d$label, levels = d$label)

p <- ggplot(d, aes(x = x, group = 1)) +
  ## null band: mean +/- 2 SD of the matched null distribution
  geom_ribbon(aes(ymin = null - 2*sd_null, ymax = null + 2*sd_null),
              fill = "grey85", alpha = 0.9) +
  geom_line(aes(y = null), colour = "grey45", linetype = "22", linewidth = 0.5) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  ## observed
  geom_line(aes(y = median), colour = "#B2182B", linewidth = 0.7) +
  geom_point(aes(y = median), colour = "#B2182B", size = 2.6) +
  geom_text(aes(y = median, label = sprintf("%+.4f", median)),
            vjust = -1.3, size = 2.7, colour = "#B2182B") +
  geom_text(aes(y = null, label = sprintf("%+.4f", null)),
            vjust = 2.0, size = 2.5, colour = "grey35") +
  annotate("text", x = 1, y = 0.0345, hjust = 0, size = 2.7,
           label = "CHD8 targets") +
  annotate("text", x = 1, y = -0.0125, hjust = 0, size = 2.5, colour = "grey35",
           label = "expression-matched null (band: \u00b12 SD)") +
  scale_y_continuous(limits = c(-0.016, 0.037),
                     breaks = seq(-0.01, 0.03, 0.01)) +
  labs(x = "CHD8 sites, ranked by minimum fold enrichment across antibody preparations",
       y = expression("median log"[2]*" fold change (ASD vs control)")) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title = element_text(size = 8))

ggsave(file.path(OUT, "FigX_CHD8_dose_response.pdf"), p, width = 4.6, height = 3.4)
ggsave(file.path(OUT, "FigX_CHD8_dose_response.png"), p, width = 4.6, height = 3.4, dpi = 300)
cat("wrote FigX_CHD8_dose_response (pdf + png) to", OUT, "\n")
cat("\nAll four points are p < 0.0001 against the matched null;\n")
cat("state that in the legend rather than on the plot.\n")
