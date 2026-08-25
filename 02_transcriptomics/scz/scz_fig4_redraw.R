.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
if (isTRUE(capabilities("cairo"))) options(bitmapType = "cairo")
suppressMessages(library(ggplot2))
OUT <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final/Transcriptomics_SCZ/figures"

dl <- data.frame(
  subset = factor(c("all 78 samples","CA1+CA3 (52)","CA1 (26)","CA3 (26)","DG (26)"),
           levels=c("all 78 samples","CA1+CA3 (52)","CA1 (26)","CA3 (26)","DG (26)")),
  donors = c("26 donors\n3 samples each","26 donors\n2 samples each",
             "26 donors\n1 sample each","26 donors\n1 sample each","26 donors\n1 sample each"),
  genes  = c(146, 50, 6, 0, 0))

p <- ggplot(dl, aes(subset, genes)) +
  geom_col(fill="#922B21", width=.65) +
  geom_text(aes(label=genes), vjust=-.5, size=3.8) +
  geom_text(aes(label=donors, y=-14), size=2.7, colour="grey35", lineheight=.95) +
  coord_cartesian(ylim=c(-22, 160), clip="off") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=15, hjust=1),
        plot.margin=margin(6, 10, 6, 10)) +
  labs(x=NULL, y="DE genes (FDR < 0.05)",
    title="GSE138082: apparent signal scales with samples per donor, not with donors",
    subtitle="Identical model and the same 26 donors (13 control, 13 case) in every analysis;\nonly the number of samples contributed per donor differs")

ggsave(file.path(OUT,"Fig_4_subfield_pseudoreplication.pdf"), p, width=8.5, height=5.8)
ggsave(file.path(OUT,"Fig_4_subfield_pseudoreplication.png"), p, width=8.5, height=5.8, dpi=150)
cat("done\n")
