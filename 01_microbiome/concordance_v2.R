setwd("/rds/projects/e/elhamsak-dissertation-neuropsychiatric/Microbiome/PRJNA880407")
last_tok <- function(x) sub(".*;", "", sub(";+_*$", "", gsub(";__", "", x)))
drop_unresolved <- function(k) grepl("uncultured|^(d|p|c|o|f)__|^Bacteria$|^$|^f__|^o__|^c__|^p__", k)
a <- read.delim("ancombc2_SZ_vs_HC.tsv", check.names=FALSE); a <- a[which(a$significant), ]
a$key <- last_tok(a$genus)
m <- read.delim("maaslin2/maaslin2_SZ_vs_HC.tsv", check.names=FALSE); m <- m[m[[8]] < 0.05, ]
m$key <- last_tok(gsub("\\.", ";", m[[1]]))
cat("ANCOM-BC2 sig:", nrow(a), " (unresolved dropped:", sum(drop_unresolved(a$key)), ")\n")
cat("MaAsLin2   sig:", nrow(m), " (unresolved dropped:", sum(drop_unresolved(m$key)), ")\n")
a <- a[!drop_unresolved(a$key), ]; m <- m[!drop_unresolved(m$key), ]
da <- a$key[duplicated(a$key)]; dm <- m$key[duplicated(m$key)]
if (length(da) || length(dm)) { cat("!! duplicate keys remain:\n"); print(unique(c(da,dm))) }
a <- a[!duplicated(a$key), ]; m <- m[!duplicated(m$key), ]
j <- merge(a[, c("key","genus","log2FC_SZ_vs_HC","q_value")],
           m[, c("key", names(m)[4], names(m)[8])], by="key")
names(j)[5:6] <- c("maaslin_coef","maaslin_q")
j$concordant <- sign(j$log2FC_SZ_vs_HC) == sign(j$maaslin_coef)
j <- j[order(j$q_value), ]
cat("\n=== named genera: intersect", nrow(j), " direction-concordant", sum(j$concordant), "===\n")
write.table(j, "robust_genera_both_methods_v2.tsv", sep="\t", quote=FALSE, row.names=FALSE)
print(j[j$concordant, c("key","log2FC_SZ_vs_HC","maaslin_coef","q_value")], row.names=FALSE)
