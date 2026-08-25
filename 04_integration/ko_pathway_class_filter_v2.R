## Restrict the KO pathway results to KEGG classes that exist in bacteria.
##
## KEGG's Human Diseases and Organismal Systems maps contain a handful of KOs
## shared with prokaryotes (chaperones, ATP synthase, ribosomal proteins), so
## PICRUSt2 assigns them non-zero abundance and they can reach significance on
## a few housekeeping genes. They are not microbial pathways and must not enter
## a cross-omics comparison. Global/overview maps are dropped for the same
## reason: they aggregate most of metabolism into a single feature.
##
## Classes come from the KEGG BRITE hierarchy br:ko00001, fetched as one flat
## file and parsed -- KEGGREST's keggGet() parser errors on some entries.
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

AX <- "/rds/projects/e/elhamsak-dissertation-autism-2/RNA_seq/DE_analysis/DISSERTATION_final/INTEGRATION_inputs/06_kegg_axis"
ARMS <- c("ASD","SCZ","MDD")

res <- lapply(ARMS, function(a)
  read.delim(file.path(AX, sprintf("KO_pathway_DA_%s.tsv", a)), check.names=FALSE))
names(res) <- ARMS
for(a in ARMS) cat(sprintf("%s: %d pathways tested\n", a, nrow(res[[a]])))

## ---- BRITE hierarchy, cached ----------------------------------------------
CLS  <- file.path(AX, "kegg_pathway_class.csv")
HTXT <- file.path(AX, "br_ko00001.keg")

if(file.exists(CLS)){
  cls <- read.csv(CLS, stringsAsFactors=FALSE)
  cat(sprintf("using cached class table (%d pathways)\n", nrow(cls)))
} else {
  if(!file.exists(HTXT)){
    cat("fetching KEGG BRITE hierarchy br:ko00001 ...\n")
    ok <- tryCatch({ download.file("https://rest.kegg.jp/get/br:ko00001",
                                   HTXT, quiet=TRUE); TRUE },
                   error=function(e){ message("download failed: ", conditionMessage(e)); FALSE })
    if(!ok) stop("could not fetch the hierarchy -- run this on a node with internet, ",
                 "or download https://rest.kegg.jp/get/br:ko00001 manually to ", HTXT)
  }
  ln <- readLines(HTXT, warn=FALSE)
  ## A/B/C prefixed lines: A = top class, B = subclass, C = pathway
  A <- B <- NA_character_
  rows <- list()
  for(l in ln){
    if(grepl("^A", l))      A <- trimws(sub("^A0?[0-9]*\\s*", "", l))
    else if(grepl("^B", l)) B <- trimws(sub("^B\\s*[0-9]*\\s*", "", l))
    else if(grepl("^C", l)){
      id <- sub("^C\\s*([0-9]{5}).*$", "\\1", l)
      if(grepl("^[0-9]{5}$", id)){
        nm <- trimws(sub("^C\\s*[0-9]{5}\\s*", "", l))
        nm <- trimws(sub("\\s*\\[PATH:.*\\]$", "", nm))
        rows[[length(rows)+1]] <- data.frame(pathway=paste0("ko", id), name=nm,
                                             top=A, sub=B, stringsAsFactors=FALSE)
      }
    }
  }
  cls <- unique(do.call(rbind, rows))
  cls <- cls[!duplicated(cls$pathway), ]
  write.csv(cls, CLS, row.names=FALSE)
  cat(sprintf("parsed %d pathways from the hierarchy -> %s\n", nrow(cls), CLS))
}

cat("\n=== all KEGG pathways by top-level class ===\n")
print(sort(table(cls$top), decreasing=TRUE))

## ---- the filter ------------------------------------------------------------
KEEP_CLASS <- c("Metabolism",
                "Genetic Information Processing",
                "Environmental Information Processing",
                "Cellular Processes")
DROP_SUB   <- c("Immune system","Nervous system","Endocrine system","Circulatory system",
                "Digestive system","Excretory system","Sensory system","Aging",
                "Development and regeneration","Environmental adaptation",
                "Global and overview maps")

cls$prokaryote_relevant <- cls$top %in% KEEP_CLASS & !(cls$sub %in% DROP_SUB)
cat(sprintf("\nprokaryote-relevant: %d of %d KEGG pathways\n",
            sum(cls$prokaryote_relevant), nrow(cls)))

## ---- apply to each arm -----------------------------------------------------
summ <- list()
for(a in ARMS){
  r <- merge(res[[a]], cls[,c("pathway","top","sub","prokaryote_relevant")],
             by.x="feature", by.y="pathway", all.x=TRUE)
  r$prokaryote_relevant[is.na(r$prokaryote_relevant)] <- FALSE
  f <- r[r$prokaryote_relevant & !r$is_global, ]
  f <- f[order(f$qval), ]
  ## BH recomputed on the retained set: the original correction included
  ## pathways now known not to be microbial
  f$qval_refiltered <- p.adjust(f$pval, method="BH")
  write.table(f, file.path(AX, sprintf("KO_pathway_DA_%s_filtered.tsv", a)),
              sep="\t", quote=FALSE, row.names=FALSE)
  summ[[a]] <- data.frame(arm=a,
    tested_all      = nrow(r),
    sig_all         = sum(r$qval < 0.05, na.rm=TRUE),
    tested_filtered = nrow(f),
    sig_original_q  = sum(f$qval < 0.05, na.rm=TRUE),
    sig_refiltered  = sum(f$qval_refiltered < 0.05, na.rm=TRUE),
    down            = sum(f$qval_refiltered < 0.05 & f$coef < 0, na.rm=TRUE),
    up              = sum(f$qval_refiltered < 0.05 & f$coef > 0, na.rm=TRUE))
  cat(sprintf("\n--- %s: top 10 retained ---\n", a))
  print(head(f[, c("feature","pathway_name","sub","coef","qval_refiltered")], 10), row.names=FALSE)
  cat(sprintf("  dropped as non-microbial: %d\n", sum(!r$prokaryote_relevant)))
}

cat("\n================ summary ================\n")
S <- do.call(rbind, summ); print(S, row.names=FALSE)
write.csv(S, file.path(AX,"KO_pathway_filter_summary.csv"), row.names=FALSE)
cat(sprintf("\nfiltered tables written to %s\n", AX))
