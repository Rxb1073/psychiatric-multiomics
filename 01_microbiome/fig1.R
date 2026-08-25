library(ggplot2)
d <- data.frame(level=factor(c("ASV","Genus"),levels=c("ASV","Genus")), shared=c(0,166))
ggplot(d, aes(level, shared, fill=level)) + geom_col(width=0.6) +
  geom_text(aes(label=shared), vjust=-0.5, size=5) +
  scale_fill_manual(values=c("ASV"="#c0392b","Genus"="#27ae60")) +
  labs(y="Features shared by all 4 cohorts", x=NULL) +
  theme_minimal(base_size=13) + theme(legend.position="none")
ggsave("fig1_overlap.png", width=4, height=4, dpi=300)
