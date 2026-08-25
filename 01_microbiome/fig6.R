library(ggplot2)
lo <- read.csv("stage5_lodo_survival.csv")
lo$genus <- sub(".*;g__","",lo$taxon); lo$genus <- sub(";s__.*","",lo$genus)
lo <- lo[order(lo$survives_of_4),]; lo$genus <- factor(lo$genus,levels=lo$genus)
lo$broken_by[is.na(lo$broken_by)] <- ""
ggplot(lo, aes(survives_of_4, genus, fill=factor(survives_of_4))) + geom_col() +
  geom_text(aes(label=ifelse(broken_by=="","",paste("drop:",broken_by))), hjust=-0.05, size=2.8) +
  scale_fill_brewer(palette="RdYlGn", name="rounds\nsurvived") +
  coord_cartesian(xlim=c(0,5)) + labs(x="Leave-one-out rounds concordant (of 4)", y=NULL) +
  theme_minimal(base_size=11)
ggsave("fig6_lodo.png", width=7, height=5, dpi=300)
