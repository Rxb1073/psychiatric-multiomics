library(ggplot2)
core <- read.csv("stage5_concordant_core.csv")
core$genus <- sub(".*;g__","",core$taxon.x); core$genus <- sub(";s__.*","",core$genus)
core$dir <- ifelse(core$ancom_lfc>0,"Enriched in ASD","Depleted in ASD")
robust <- c("[Ruminococcus]_torques_group","Agathobacter","CAG-352")
core$robust <- ifelse(core$genus %in% robust,"*","")
core <- core[order(core$ancom_lfc),]; core$genus <- factor(core$genus,levels=core$genus)
ggplot(core, aes(ancom_lfc, genus, fill=dir)) + geom_col() +
  geom_text(aes(label=robust), hjust=ifelse(core$ancom_lfc>0,-0.3,1.3), size=5) +
  scale_fill_manual(values=c("Depleted in ASD"="#2980b9","Enriched in ASD"="#c0392b")) +
  labs(x="ANCOM-BC2 log-fold change (ASD vs Control)", y=NULL, fill=NULL,
       caption="* robust core: directionally unanimous across all 4 cohorts") +
  theme_minimal(base_size=11)
ggsave("fig5_diffabund.png", width=7, height=5, dpi=300)
