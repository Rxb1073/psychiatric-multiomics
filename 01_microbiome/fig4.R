library(ggplot2)
d <- data.frame(metric=factor(c("Bray-Curtis","Jaccard","Unweighted UniFrac","Weighted UniFrac"),
  levels=c("Bray-Curtis","Jaccard","Weighted UniFrac","Unweighted UniFrac")), pct=c(4.0,9.3,12.4,11.2))
ggplot(d, aes(metric,pct,fill=metric)) + geom_col(width=0.65) +
  geom_text(aes(label=sprintf("+%.0f%%",pct)),vjust=-0.4,size=4) +
  labs(y="ASD excess dispersion (%)", x=NULL) + coord_cartesian(ylim=c(0,15)) +
  theme_minimal(base_size=12) + theme(legend.position="none",axis.text.x=element_text(angle=20,hjust=1))
ggsave("fig4_dispersion.png", width=6, height=4, dpi=300)
