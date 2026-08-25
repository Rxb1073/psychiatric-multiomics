library(ggplot2)
a <- read.csv("stage4_alpha_values.csv", row.names=1)

mk <- function(col, ylab, title){
  ggplot(a, aes(diagnosis, .data[[col]], fill=diagnosis)) +
    geom_boxplot(outlier.size=0.5) +
    labs(title=title, x=NULL, y=ylab) +
    theme_minimal(base_size=12) + theme(legend.position="none")
}
p1 <- mk("observed","Genera","Observed genera")
p2 <- mk("shannon","H'","Shannon")

# Faith's PD from the rarefied export
f <- read.delim("unifrac/faith_export/alpha-diversity.tsv", row.names=1)
fcol <- names(f)[1]                      # whatever the column is called
f$diagnosis <- a[rownames(f),"diagnosis"]
f <- f[!is.na(f$diagnosis),]
p3 <- ggplot(f, aes(diagnosis, .data[[fcol]], fill=diagnosis)) +
  geom_boxplot(outlier.size=0.5) +
  labs(title="Faith's PD", x=NULL, y="PD") +
  theme_minimal(base_size=12) + theme(legend.position="none")

# arrange 3 panels without patchwork, via a grid device
png("fig3_alpha.png", width=9, height=3.5, units="in", res=300)
gridExtra::grid.arrange(p1, p2, p3, nrow=1)
dev.off()
cat("wrote fig3_alpha.png\n")
