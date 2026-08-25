suppressMessages(library(Maaslin2))
setwd("/rds/projects/e/elhamsak-dissertation-autism/Microbiome/pooled")
path <- read.delim("pooled_pathways.tsv", row.names=1, check.names=FALSE)
meta <- read.delim("pooled_pathways_meta.tsv", row.names=1, check.names=FALSE)
meta <- meta[colnames(path), , drop=FALSE]
meta$diagnosis <- factor(meta$diagnosis, levels=c("Control","ASD"))
meta$study <- factor(meta$study)
cat("pathways:", nrow(path), " samples:", ncol(path), "\n")
fit <- Maaslin2(input_data=as.data.frame(t(path)), input_metadata=meta,
                output="stage6_maaslin_pathways",
                fixed_effects=c("study","diagnosis"),
                reference="study,PRJNA453621;diagnosis,Control",
                normalization="TSS", transform="LOG",
                min_prevalence=0.10, plot_heatmap=FALSE, plot_scatter=FALSE)
res <- fit$results; dx <- res[res$metadata=="diagnosis", ]
sig <- dx[dx$qval<0.05, ]
cat("\nSignificant pathways (q<0.05):", nrow(sig), "of", nrow(dx), "tested\n")
cat("  direction:", sum(sig$coef>0),"ASD-up,",sum(sig$coef<0),"ASD-down\n\n")
scfa <- c("PWY-5022","CENTFERM-PWY","PWY-5088","P163-PWY","PWY-5676","PWY-5677",
          "butanoate","butyrate","propanoate","PWY-6590","PWY0-1312","PWY-6588")
cat("===== SCFA / butyrate-related pathways =====\n")
hits <- dx[grepl(paste(scfa,collapse="|"), dx$feature, ignore.case=TRUE), ]
if(nrow(hits)>0){ hits <- hits[order(hits$qval), ]
  for(i in seq_len(nrow(hits)))
    cat(sprintf("  %-16s coef=%+.3f q=%.3g %s\n", hits$feature[i], hits$coef[i], hits$qval[i],
        ifelse(hits$qval[i]<0.05, ifelse(hits$coef[i]<0,"** ASD-DOWN **","** ASD-UP **"),"(ns)")))
} else cat("  (no SCFA-named pathways matched)\n")
write.csv(dx, "stage6_pathway_results.csv", row.names=FALSE)
cat("\n[NSTI CAVEAT] LOW CONFIDENCE: mean NSTI 0.08-0.20 (529598=0.20, 1022510=0.17 exceed 0.15).\n")
